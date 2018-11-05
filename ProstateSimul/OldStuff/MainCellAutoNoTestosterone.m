%My try at some sort of really simple cellular automata model. Two rules:
%the cell spreads to it's Moore neighborhood or dies. 
%Daniel Wu
%6/17/16

clearvars;clf;clc
%filename = 'testAutomata1.gif';

fieldSize = 1000;
cancerStartSize = 10;

%Here's the playing field!
field = zeros(fieldSize);
%And the testosterone field
testoField = zeros(fieldSize);

%This is how many iterations we are going through.
time = 1000;

%This is the chance a cancer cell will die, [t+ t- tp]
OrigDeathChance = [0.3 0.3 0.3];
deathChance = OrigDeathChance;              %This is changed by drugs
%This is the chance a cell will fill a nearby empty one.
growthChance = [0.1, 0.07, 0.09];
testoProduction = 1;            %Units of testosterone per tp per step
testoDecay = 0.02;

%Now we do the same for T-


%Drug function
dose = 300;         %Some sort of arbitrary dose.
drugDecay = 5;      %Some sort of decay, dose/step.
drugLevel = 0;      %Current amount of drugs in system.
drugDeath = [0.001, 0, 0.001];      %How deadly the drugs are to [T+, T-, Tp]

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


%Now we make a pretty interface
figure(1)
button = uicontrol('string','Hit them with the drugs!','Callback','drugLevel = drugLevel + dose;');
button.Position(2) = 10;                %Make the button stop covering the plot
button.Position(3) = 150;               %Make the button longer

drugDisplay = uicontrol('style','text');
drugDisplay.Position(1) = 10;
drugDisplay.Position(2) = 60;
drugDisplay.Position(3) = 80;

timeCount = uicontrol('style','text','Position',[10 100 80 20]);

for i = 1:time                          %For the entire simulation
   [cellX, cellY] = find(field);        %Tells me where all the crap is
   randomOrder = randperm(length(cellX));        %A random order
   cellX = cellX(randomOrder);                  %Shuffle the order of cells
   cellY = cellY(randomOrder);
   for j=1:length(cellX)                %For every cell
      neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around it
      cellType = field(cellX(j),cellY(j));              %And tell me what kind of cell it is. 
      luck = rand(3);
      for xTest = 1:3
           for yTest = 1:3              %And for that stuff
            
               if neighborhood(xTest,yTest) == 0        %If there's nothing living there
                   if luck(xTest,yTest)<growthChance(cellType)               %And you're lucky
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
   deathChance = OrigDeathChance + drugLevel.*drugDeath;
   drugLevel = max(drugLevel-drugDecay,0);
   set(drugDisplay, 'string',strcat('Drug Level: ', num2str(drugLevel)));
   set(timeCount, 'string', strcat('Time: ', num2str(i)));
   
   %Hard part: updating testosterone levels. Still broken for now.
   %testoField = diffuseTesto(testoField);
   
   spy(field==1)
   hold on
   spy(field==2,'r')
   spy(field==3,'g')
   pause(0.01)
   
   % Do I want to .gif this simulation? If so, uncomment these lines and
   % the one at the top that says filename.
   %frame = getframe(1);
   %   im = frame2im(frame);
   %   [imind,cm] = rgb2ind(im,256);
   %   if i == 1;
   %       imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
   %   else
   %       imwrite(imind,cm,filename,'gif','WriteMode','append');
   %   end
   hold off
end