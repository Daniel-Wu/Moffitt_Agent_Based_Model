%Makes a whole lot of simulations
%For macs, make the comp stop turning off
clearvars;
system('caffeinate -ims &')

simulNum = 1;
folderName = 'simulData';
mkdir(folderName);

%%Start simul iterations
for simulCount = 13:13%1:simulNum
    %% Parameters
    clc;

    param.fieldSize = 100;
    param.cancerStartSize = 10;

    %This is how many iterations we are going through.
    time = 1500;

    %Treatment Schedules
    drugSchedule = [300 600 900 1200];
    testoSchedule = [];%[500 700 900 1100 1300];

    %This is the chance a cancer cell will die, [t+ t- tp]
    param.deathChance = [0.02 0.02 0.02];
    %This is the chance a cell will fill a nearby empty one.
    param.growthChanceConstants = [0, 0.05, 0];
    
    % PSA
    param.PSAProd = [0.01 0.006, 0.006];
    param.PSADecay = 0.3;
    
    param.origTestoProduction = 1;        %units of testosterone per tp per step
    testoProduction = param.origTestoProduction;            %changed by Abi
    param.testoDecay = 0.02;
    testoHelpMat = zeros(param.fieldSize);
    param.testoHelpMatSpreadRadius = 2;               %How "much" diffusion occurs, how round the surface is
    param.testoDose = 5;
    param.initTesto = 2;
    %In this model we have the assumption that Tp benefits less from the
    %testo than t+ due to the cost of production
    param.TpTestoDisad = 2/3;     %How much testo benefits Tp compared to t+

    %Abi drug function
    param.drugDose = 500;         %Some sort of arbitrary param.drugDose.
    param.drugDecay = 0.02;      %Some sort of decay, param.drugDose/step.
    param.initDrug = 0;
    param.drugEfficacy = 0.005;     %How much testo production per unit of drug that is shut down
    drugLevel = param.initDrug;      %Current amount of drugs in system.

    %% Internal Initialization

    %Here's the playing field!
    field = zeros(param.fieldSize);
    fieldData = zeros(param.fieldSize, param.fieldSize, time);
    %And the testosterone field, with some inital testo
    testoField = zeros(param.fieldSize) + param.initTesto;
    testoFieldData = zeros(param.fieldSize,param.fieldSize, time);

    %Indexes of cells that aren't the center of a 3x3 matrix
    edgeInd = [1 2 3 4 6 7 8 9];

    %Initialize growth chances
    origGrowthChance = zeros(param.fieldSize,param.fieldSize,3);
    origGrowthChance(:,:,1) = param.growthChanceConstants(1);
    origGrowthChance(:,:,2) = param.growthChanceConstants(2);
    origGrowthChance(:,:,3) = param.growthChanceConstants(3);
    growthChance = origGrowthChance;

    %Recording populations and drug/testo/PSA levels
    TMPop = zeros(1,time);
    TIPop = zeros(1,time);
    TPPop = zeros(1,time);
    testoTracker = zeros(1,time);
    drugTracker = zeros(1,time);
    PSATracker = zeros(1,time);
    PSA = 0;

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


    %% Simulation

    for i = 1:time                          %For the entire simulation
       [cellX, cellY] = find(field>0);        %Tells me where all the crap is
       randomOrder = randperm(length(cellX));        %A random order
       cellX = cellX(randomOrder);                  %Shuffle the order of cells
       cellY = cellY(randomOrder);
       for j=1:length(cellX)                %For every cell
            
          cellType = field(cellX(j),cellY(j));              %Tell me what kind of cell it is. 
          PSA = PSA + param.PSAProd(cellType) * testoHelpMat(cellX(j),cellY(j));      %Make PSA
          growthFactor = growthChance(cellX(j),cellY(j),cellType);
          luck = rand();

          if luck<growthFactor                              %If you are lucky, you divide
             neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around the cell
             edgeInd = edgeInd(randperm(length(edgeInd)));              %Scramble checking order
             for c=1:8                                                  %For this order
                 if neighborhood(edgeInd(c)) ==0                         %If there isn't anything there
                     neighborhood(edgeInd(c)) = cellType;                %Divide there
                     break;                                             %And stop checking
                 end
             end
             field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in
          end

           %If it's tp, it'll make testosterone
           if cellType == 3
               testoField(cellX(j),cellY(j)) = testoField(cellX(j),cellY(j)) + testoProduction;
           end

           %But the cell might die
           if rand()<param.deathChance(cellType)
              field(cellX(j),cellY(j)) = 0; 
           end
       end

       %RecalcDrugs
       testoProduction = max(param.origTestoProduction-drugLevel*param.drugEfficacy,0);     %Shuts off production of testosterone
       drugLevel = max(drugLevel*(1-param.drugDecay),0);                          %Decay Drugs
       drugTracker(i) = drugLevel;                                      %Record Drug Level

       %Updating testosterone levels. 
       %testoField = diffuseTesto(testoField);  Really slow masstransfer formula
       testoField = spreadTesto(testoField,param.testoDecay);   %Faster version made by me
       testoField = spreadTesto(testoField,0);  %Spread it twice?
       %testoField = spreadTesto(testoField,0);  %Spread it thrice?
       testoTracker(i) = sum(sum(testoField)); 

       for a = (param.testoHelpMatSpreadRadius+1):(param.fieldSize-param.testoHelpMatSpreadRadius)
           for b = (param.testoHelpMatSpreadRadius+1):(param.fieldSize-param.testoHelpMatSpreadRadius)
               testoHelpMat(a,b) = sum(sum(testoField((a-param.testoHelpMatSpreadRadius):(a+param.testoHelpMatSpreadRadius),(b-param.testoHelpMatSpreadRadius):(b+param.testoHelpMatSpreadRadius))));       %Cells are affected in a 5x5 square?
           end
       end
       testoHelpMat = testoHelpMat./1000;        %BSed constant of division for t+
       %Reassign growth chance based on testosterone
       growthChance(:,:,1) = origGrowthChance(:,:,1) + testoHelpMat;
       growthChance(:,:,3) = origGrowthChance(:,:,3) + testoHelpMat*param.TpTestoDisad;

       %Alter PSA
       PSA = PSA*(1-param.PSADecay);
       PSATracker(i) = PSA;
       
       %Apply Treatment
       if any(drugSchedule==i)              %Abi first
           drugLevel = drugLevel + param.drugDose;
       end

       if any(testoSchedule==i)             %Testo next
           testoField = testoField + param.testoDose;
       end

       %SAVE DATA
       fieldData(:,:,i) = field;
       testoFieldData(:,:,i) = testoHelpMat;
       % record cell populations
       TMPop(i) = length(find(field==1));
       TIPop(i) = length(find(field==2));
       TPPop(i) = length(find(field==3));
    end


    %SAVE DATA
    save([folderName '/simulation' num2str(simulCount)],'fieldData','testoFieldData','TMPop','TIPop','TPPop','testoTracker','drugTracker','PSATracker', 'drugSchedule', 'testoSchedule','param')
    disp(['Simulation ' num2str(simulCount) ' done!'])
    beep
end

system('killall caffeinate &')