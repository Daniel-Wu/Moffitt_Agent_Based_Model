%Makes a whole lot of simulations
%For macs, make the comp stop turning off
clearvars; clc;
%system('caffeinate -ims &')

simulStartNum = 1;
simulEndNum = 10;
folderName = 'C:\Users\User\Desktop\simulDatas';
%If the folder doesn't exist, make it
if ~exist('folderName','dir')
    mkdir(folderName);
end
fileName = 'Naive';

%% Start simul iterations
parfor simulCount = simulStartNum:simulEndNum
    
    %Create save space for pictures
    mkdir([folderName '\' fileName num2str(simulCount) 'Pictures']);
    
    %% Parameters
    param = struct();
    param.fieldSize = 400;
    param.cancerStartSize = 10;

    %This is how many iterations we are going through.
    time = 5000;%1500;

    %Treatment Schedules
    AbiSchedule = [];
    testoSchedule = [];
    doceSchedule = [];

    %This is the chance a cancer cell will die, [t+ t- tp]]
    param.deathChance = [1/60, 1/60, 1/60];

    %This is the base chance a cell will fill a nearby empty one.
    param.growthChanceConstants = [0, 1/50, 0];

    %This is the chance a cell will move into a random nearby empty spot
    param.moveChance = [1/4, 1/4, 1/4];
    % PSA
    param.PSAProd = [0.01 0.006, 0.006];
    param.PSADecay = 0.3;
    
    decayRate = 0.02;
    
    %testosterone
    param.origTestoProduction = 1;        %units of testosterone per tp per step
    testoProduction = param.origTestoProduction;            %changed by Abi
    param.testoDecay = decayRate;
    testoHelpMat = zeros(param.fieldSize);
    param.testoHelpMatSpreadRadius = 3;               %How "much" diffusion occurs, how round the surface is
    param.testoDose = 20;
    param.initTesto = 5;
    %In this model we have the assumption that Tp benefits less from the
    %testo than t+ due to the cost of production
    param.TpTestoDisad = 3/4;     %How much testo benefits Tp compared to t+

    %Abi drug function
    %%
    param.Dose = 50000;%50;%1500;         %Some sort of arbitrary param.drugDose.
    %%
    param.drugDecay = decayRate;      %Some sort of decay, param.drugDose/step.
    param.initDrug = 0;
    param.drugEfficacy = 0.005;     %How much testo production per unit of drug that is shut down
    param.minAbi = 0;     %Minimum abi level - good for modelling continuous abi
    drugLevel = param.initDrug;      %Current amount of drugs in system.

    %DOCETAXEL  
    %%
    param.doceDose = 0;%1000;
    %%
    param.doceDecay = decayRate;
    param.doceEfficacy = 0.001;
    param.initDoce = 0;
    doceLevel = param.initDoce;
    doceDeathChance = doceLevel*param.doceEfficacy;
    
    %% Internal Initialization

    %Here's the playing field!
    field = zeros(param.fieldSize, 'int8');
    %And the testosterone field, with some inital testo
    testoField = zeros(param.fieldSize) + param.initTesto;
    %Indexes of cells that aren't the center of a 3x3 matrix
    edgeInd = [1 2 3 4 6 7 8 9];

    %Initialize growth chances
    origGrowthChance = zeros(param.fieldSize,param.fieldSize,3);
    origGrowthChance(:,:,1) = param.growthChanceConstants(1);
    origGrowthChance(:,:,2) = param.growthChanceConstants(2);
    origGrowthChance(:,:,3) = param.growthChanceConstants(3);
    growthChance = origGrowthChance;

    %Recording populations and drug/testo/PSA levels
    TMPop = nan(1,time);
    TIPop = nan(1,time);
    TPPop = nan(1,time);
    totalPop = nan(1,time);
    testoTracker = nan(1,time);
    drugTracker = nan(1,time);
    doceTracker = nan(1,time);
    PSATracker = nan(1,time);
    PSA = 0;

    %Adaptive treatment
    wantToDoAdaptTreatment = false;
    offsetTimer = [300 300 300];    %[abi, testo, doce]
    %%
    offsetTimerResetVal = [10 300 300]; %Min time btw adaptive therapy
    %%
    PSAStartAbi = 8;
    TIPorpStartDoce = 0.33;
    
    %Here's a nice group of random starting cells.
    middleLoc = round(size(field)/2);
    startField = field(middleLoc:middleLoc+param.cancerStartSize-1,middleLoc:middleLoc+param.cancerStartSize-1);

    %Now we randomly give those starting cells a phenotype.

    for i = 1:length(startField)
       for j = 1:length(startField) 
            startField(i,j) = ceil(rand()*3);
       end
    end
    %And then we put it back into the field.
    field(middleLoc:middleLoc+param.cancerStartSize-1,middleLoc:middleLoc+param.cancerStartSize-1) = startField;

    %Now we build a wall on the border of the field.
    %I wish I could just use padarray here.
    field(1,:) = -1;
    field(param.fieldSize,:) = -1;
    field(:,1) = -1;
    field(:,param.fieldSize) = -1;
    %padarray(field, [1 1], -1);

    %% Simulation

    for i = 1:time                          %For the entire simulation
       %Make a list of cells
       [cellX, cellY] = find(field>0);        %Tells me where all the crap is
       randomOrder = randperm(length(cellX));        %A random order
       cellX = cellX(randomOrder);                  %Shuffle the order of cells
       cellY = cellY(randomOrder);
       %% For Every Cell
       for j=1:length(cellX)
            
          cellType = field(cellX(j),cellY(j));              %Tell me what kind of cell it is. 
          PSA = PSA + param.PSAProd(cellType) * testoHelpMat(cellX(j),cellY(j));      %Make PSA
          
          %% Movement
          if rand()<param.moveChance(cellType)  %If you move
              neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around the cell
              edgeInd = edgeInd(randperm(length(edgeInd)));              %Scramble checking order
              for iMove = 1:8               %For the spots in that order
                  if neighborhood(edgeInd(iMove))==0            %If it's empty
                      neighborhood(edgeInd(iMove)) = cellType;  %Move there
                      neighborhood(5) = 0;
                      break;
                  end
              end
              field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in

          end
          
          
          %% Growth
          growthFactor = growthChance(cellX(j),cellY(j),cellType);

          if rand()<growthFactor                              %If you are lucky, you try to divide
             if rand()>doceDeathChance                       %If you are more lucky, you survive Dox and divide
                 neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around the cell
                 edgeInd = edgeInd(randperm(length(edgeInd)));              %Scramble checking order
                 for c=1:8                                                  %For this order
                     if neighborhood(edgeInd(c)) ==0                         %If there isn't anything there
                         neighborhood(edgeInd(c)) = cellType;                %Divide there
                         break;                                             %And stop checking
                     end
                 end
                 field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in
             
             else                           %If you aren't that lucky, you die from dox.
                 field(cellX(j), cellY(j)) = 0;
                 
             end
          end
          %% Testosterone Production
          %If it's tp, it'll make testosterone
           if cellType == 3
               testoField(cellX(j),cellY(j)) = testoField(cellX(j),cellY(j)) + testoProduction;
           end
          %% Death
          %But the cell might die
           if rand()<param.deathChance(cellType)
              field(cellX(j),cellY(j)) = 0; 
           end
       end
%% End of For Every Cell
       %RecalcDrugs
       testoProduction = max(param.origTestoProduction-drugLevel*param.drugEfficacy,0);     %Shuts off production of testosterone
       if ~isempty(AbiSchedule) && i>AbiSchedule(1)
            drugLevel = max(drugLevel*(1-param.drugDecay),param.minAbi);                          %Decay Drugs bounded at minAbi
       else
            drugLevel = max(drugLevel*(1-param.drugDecay), 0);                          %Decay Drugs bounded at 0
       end
       drugTracker(i) = drugLevel;                                      %Record Drug Level

       %docetaxel
       doceLevel = max(doceLevel*(1-param.doceDecay),0);    %Decay Docetaxel
       doceDeathChance = doceLevel*param.doceEfficacy;
       doceTracker(i) = doceLevel;      %Record Doce Level
       
       %Updating testosterone levels. 
       %testoField = diffuseTesto(testoField);  Really slow masstransfer formula
       testoField = spreadTesto(testoField,param.testoDecay);   %Faster (maybe accurate?) version made by me
       testoField = spreadTesto(testoField,0);  %Spread it twice?
       %testoField = spreadTesto(testoField,0);  %Spread it thrice?
       testoTracker(i) = sum(sum(testoField)); 

       for a = (param.testoHelpMatSpreadRadius+1):(param.fieldSize-param.testoHelpMatSpreadRadius)
           for b = (param.testoHelpMatSpreadRadius+1):(param.fieldSize-param.testoHelpMatSpreadRadius)
               testoHelpMat(a,b) = sum(sum(testoField((a-param.testoHelpMatSpreadRadius):(a+param.testoHelpMatSpreadRadius),(b-param.testoHelpMatSpreadRadius):(b+param.testoHelpMatSpreadRadius))));       %Cells are affected in a 5x5 square?
           end
       end
       testoHelpMat = testoHelpMat./5000;        %BSed constant of division for t+
       
       %Logarithmic adjustment to testo's effects.
       testoHelpMat = calcTestoHelpMat(testoHelpMat);
       
       %Reassign growth chance based on testosterone
       growthChance(:,:,1) = origGrowthChance(:,:,1) + testoHelpMat;
       growthChance(:,:,3) = origGrowthChance(:,:,3) + testoHelpMat*param.TpTestoDisad;

       %Alter PSA
       PSA = PSA*(1-param.PSADecay);
       PSATracker(i) = PSA;
       
       %% Export data to .png
       fig = figure('Visible', 'off');
       set(fig,'units','normalized','outerposition',[0 0 1 1])
       
       %Field Display
       timeDisplay = uicontrol('style','text','Position',[30 120 100 20]);
       set(timeDisplay, 'string', strcat('Time: ', num2str(i)));
       subplot(2,3,[1 2 4 5])
       hold on
       set(gca, 'Box', 'on', 'YTick', [], 'XTick', []);
       spy(field == 1, 'k')
       spy(field == 2, 'r')
       spy(field == 3, 'g')
       xlabel('')
       title('Field')
       hold off
       
       %Androgen Display
       subplot(2,3,3)
       %surf(testoHelpMat)                       %Surface display
       pcolor(flipud(testoHelpMat))              %Flat Display
       axis off;
       title('Androgen Field')
       shading interp;
       drawnow;
       
       %Population Display
       subplot(2, 3, 6)
       hold on;
       plot(TIPop, 'r', 'lineWidth', 2)
       plot(TPPop, 'g', 'lineWidth', 2)
       plot(TMPop, 'k', 'lineWidth', 2)
       lgd = legend('Independent', 'Producing', 'Consuming');
       set(lgd, 'Position', [0.5, 0.1, 0.1, 0.05])
       title('      Cell Population')
       xlabel('Simulation Time')
       ylabel('Population')
       
       saveas(fig, [folderName '\' fileName num2str(simulCount) 'Pictures\'  num2str(i) '.png']);


       
       % record cell populations
       TMPop(i) = length(find(field==1));
       TIPop(i) = length(find(field==2));
       TPPop(i) = length(find(field==3));
       totalPop(i) = TMPop(i) + TIPop(i) + TPPop(i);
       
       %Apply Scheduled Treatment
       if any(AbiSchedule==i)              %Abi first
           drugLevel = drugLevel + param.Dose;
       end

       if any(testoSchedule==i)             %Testo next
           testoField = testoField + param.testoDose;
       end
       
       if any(doceSchedule==i)
           doceLevel = doceLevel + param.doceDose;
       end

       %Apply Adaptive Treatment
       if wantToDoAdaptTreatment
           offsetTimer = offsetTimer - 1;  %Decrement timer
           if PSA>PSAStartAbi && offsetTimer(1)<=0
               offsetTimer(1) = offsetTimerResetVal(1);    %Reset timer
               drugLevel = drugLevel + param.Dose;         %Give drugs
               AbiSchedule = [AbiSchedule i];              %Record administration
           end
           
           %Do I need adaptive testo here?
           
           if TIPop(i)/totalPop(i) > TIPorpStartDoce && offsetTimer(3)<=0
               offsetTimer(3) = offsetTimerResetVal(3);     %Reset Timer
               drugLevel = drugLevel + 1000;      %Give doce
               AbiSchedule = [AbiSchedule i];
               doceSchedule = [doceSchedule, i+50];             %Record admin
               testoSchedule = [testoSchedule, i+100];       %Give helpful testo
           end
       end
    end %End of for each timestep


    %SAVE DATA
    %save([folderName '/' fileName num2str(simulCount)],'fieldData','testoFieldData','TMPop','TIPop','TPPop','testoTracker','drugTracker','doceTracker','PSATracker', 'AbiSchedule', 'testoSchedule','doceSchedule','param','-v7.3')
    cheatSave(folderName, fileName, simulCount, TMPop, TIPop, TPPop, testoTracker, drugTracker, doceTracker, PSATracker, AbiSchedule, testoSchedule, doceSchedule, param);
    disp(['Simulation ' num2str(simulCount) ' done!'])
    beep;
end

%system('killall caffeinate &')