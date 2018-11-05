%My try at some sort of really simple cellular automata model. Two rules:
%the cell spreads to it's Moore neighborhood or dies. Except now there is
%testosterone. Daniel Wu 6/17/16

%% Parameters
clearvars;clf;clc;
%filename = 'ProstateSuicide.gif';

fieldSize = 100;
cancerStartSize = 10;

%This is how many iterations we are going through.
time = 1500;

%Treatment Schedules
drugSchedule = [];%[300 600 900 1200];
testoSchedule = [];%[500 700 900 1100 1300];
%This is the chance a cancer cell will die, [t+ t- tp]
deathChance = [0.05 0.05 0.05];
%This is the chance a cell will fill a nearby empty one.
growthChanceConstants = [0.005, 0.02, 0.005];


origTestoProduction = 1;        %units of testosterone per tp per step
testoProduction = origTestoProduction;            %changed by Abi
testoDecay = 0.02;
testoHelpMat = zeros(fieldSize);
testoDose = 10;

%Abi drug function
dose = 400;         %Some sort of arbitrary dose.
drugDecay = 0.02;      %Some sort of decay, porportionDose/step.
drugLevel = 0;      %Current amount of drugs in system.

%% Internal Initialization

%Here's the playing field!
field = zeros(fieldSize);
%And the testosterone field
testoField = zeros(fieldSize);

%Initialize growth chances
origGrowthChance = zeros(fieldSize,fieldSize,3);
origGrowthChance(:,:,1) = growthChanceConstants(1);
origGrowthChance(:,:,2) = growthChanceConstants(2);
origGrowthChance(:,:,3) = growthChanceConstants(3);
growthChance = origGrowthChance;

%Recording populations and drug levels
TMPop = zeros(1,time);
TIPop = zeros(1,time);
TPPop = zeros(1,time);
testoTracker = zeros(1,time);
drugTracker = zeros(1,time);

%Here's a nice 4x4 group of random starting cells.
middleLoc = round(size(field)/2);
startField = field(middleLoc:middleLoc+cancerStartSize-1,middleLoc:middleLoc+cancerStartSize-1);

%Now we randomly give those starting cells a phenotype.

for i = 1:length(startField)
   for j = 1:length(startField) 
        startField(i,j) = ceil(rand()*3);
   end
end
%And then we put it back into the field.
field(middleLoc:middleLoc+cancerStartSize-1,middleLoc:middleLoc+cancerStartSize-1) = startField;

%Now we build a wall on the border of the field.
%I wish I could just use padarray here.
field(1,:) = -1;
field(fieldSize,:) = -1;
field(:,1) = -1;
field(:,fieldSize) = -1;

%% Now we make a pretty interface
figure(1)
%Drug button
button = uicontrol('string','Hit them with Abi!','Callback','drugLevel = drugLevel + dose;');
button.Position(2) = 10;                %Make the button stop covering the plot
button.Position(3) = 150;               %Make the button longer

%Testosterone button
testoButton = uicontrol('string','Throw in testosterone!', 'Callback', 'testoField = testoField + testoDose;');
testoButton.Position(1) = 200;
testoButton.Position(2) = 10;                %Make the button stop covering the plot
testoButton.Position(3) = 150;               %Make the button longer

%Drug quantity display
drugDisplay = uicontrol('style','text','Position',[20 60 80 20]);

%Timer
timeCount = uicontrol('style','text','Position',[20 100 80 20]);

%Figure
hold on
subplot(2,3,[1 2 4 5])

%Legend
legendTM = uicontrol('style','text','Position',[20 140 80 20],'ForegroundColor','b');
legendTI = uicontrol('style','text','Position',[20 160 80 20],'ForegroundColor','r');
legendTP = uicontrol('style','text','Position',[20 180 80 20],'ForegroundColor','g');
set(legendTM, 'string', 'T+ Cheaters');
set(legendTI, 'string', 'T- Independents');
set(legendTP, 'string', 'Tp Producers');

%Record parameters maybe
%parameterSummary = [strcat('Drug Dosage Times: ', num2str(drugSchedule)); strcat(); strcat(); strcat; strcat(); strcat()];
%% Simulation

for i = 1:time                          %For the entire simulation
   [cellX, cellY] = find(field>0);        %Tells me where all the crap is
   randomOrder = randperm(length(cellX));        %A random order
   cellX = cellX(randomOrder);                  %Shuffle the order of cells
   cellY = cellY(randomOrder);
   for j=1:length(cellX)                %For every cell
      neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around it
      cellType = field(cellX(j),cellY(j));              %And tell me what kind of cell it is. 
      luck = rand(3);
      growthFactor = growthChance(cellX(j),cellY(j),cellType);
      for xTest = 1:3
           for yTest = 1:3              %And for that stuff
                        
               if neighborhood(xTest,yTest) == 0        %If there's nothing living there
                   if luck(xTest,yTest)<growthFactor               %And you're lucky
                      neighborhood(xTest,yTest) = cellType;    %Now there's something there. 
                   end
               end
            
           end
       end
       
       field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in
       
       %If it's tp, it'll make testosterone
       if cellType == 3
           testoField(cellX(j),cellY(j)) = testoField(cellX(j),cellY(j)) + testoProduction;
       end
       
       %But the cell might die
       if rand()<deathChance(cellType)
          field(cellX(j),cellY(j)) = 0; 
       end
   end
    
   %RecalcDrugs
   testoProduction = max(origTestoProduction-drugLevel/200,0);     %Shuts off production of testosterone
   drugLevel = max((1-drugDecay)*drugLevel,0);
   drugTracker(i) = drugLevel;
   set(drugDisplay, 'string',strcat('Abi Level: ', num2str(drugLevel)));
   set(timeCount, 'string', strcat('Time: ', num2str(i)));
   
   %Updating testosterone levels. 
   %testoField = diffuseTesto(testoField);  Really slow masstransfer formula
   testoField = spreadTesto(testoField,testoDecay);   %Faster version made by me
   testoField = spreadTesto(testoField,0);  %Spread it twice?
   testoField = spreadTesto(testoField,0);  %Spread it thrice?
   testoTracker(i) = sum(sum(testoField)); 
   
   for a = 2:(fieldSize-1)
       for b = 2:(fieldSize-1)
           testoHelpMat(a,b) = sum(sum(testoField((a-1):(a+1),(b-1):(b+1))));
       end
   end
   testoHelpMat = testoHelpMat./500;        %BSed constant of division for t+
   %Reassign growth chance based on testosterone
   growthChance(:,:,1) = origGrowthChance(:,:,1) + testoHelpMat;
   growthChance(:,:,3) = origGrowthChance(:,:,3) + testoHelpMat./2; %BSed constant for Tp
   
   %Apply Treatment
   if any(drugSchedule==i)              %Abi first
       drugLevel = drugLevel + dose;
   end
   
   if any(testoSchedule==i)             %Testo dose
       testoField = testoField + testoDose;
   end
   
   %plot results nicely
   spy(field==1,25)
   hold on
   spy(field==2,'r',25)
   spy(field==3,'g',25)
   pause(0.01)
   subplot(2,3,3)
   surf(testoHelpMat)
   subplot(2,3,[1 2 4 5])
   
   % record cell populations
   TMPop(i) = length(find(field==1));
   TIPop(i) = length(find(field==2));
   TPPop(i) = length(find(field==3));
   
   
   % Do I want to .gif this simulation? If so, uncomment these lines and
   % the one at the top that says filename.
%    frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if i == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
   hold off
end

%% Results Figure
figure(2);
hold on;
totalPop = TMPop + TIPop + TPPop;
plot(TMPop);
plot(TIPop,'r');
plot(TPPop,'g');
plot(totalPop,'m');
plot(testoTracker/20,'c')   %Scale to fit.
plot(drugTracker)
title('Cancer Population')
legend('Consumers','Independents','Producers','Total Cancer Population','Testosterone','Abi Level')
