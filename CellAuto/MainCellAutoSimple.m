%My try at some sort of really simple cellular automata model. Two rules:
%the cell spreads to it's Moore neighborhood or dies.

fieldSize = 1000;
cancerStartSize = 2;

%Here's the playing field!
field = zeros(fieldSize);

%This is how many iterations we are going through.
time = 2000;

%This is the chance a cell will die
OrigDeathChance = 0.4;
deathChance = OrigDeathChance;              %This is changed by drugs
%This is the chance a cell will fill a nearby empty one.
growthChance = 0.1;

%Drug function
dose = 300;         %Some sort of arbitrary dose.
drugDecay = 5;      %Some sort of decay, dose/step.
drugLevel = 0;      %Current amount of drugs in system.

%Here's a nice 2x2 group of starting cells.
middleLoc = round(size(field)/2);
field(middleLoc:middleLoc+cancerStartSize-1,middleLoc:middleLoc+cancerStartSize-1) = 1;

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
   
   for j=1:length(cellX)                %For every cell
      neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around it
       
       for xTest = 1:3
           for yTest = 1:3              %And for that stuff
            
               if neighborhood(xTest,yTest) == 0        %If there's nothing living there
                   if rand()<growthChance               %And you're lucky
                      neighborhood(xTest,yTest) = 1;    %Now there's something there. 
                   end
               end
            
           end
       end
       
       field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in
       
       %But the cell might die
       if rand()<deathChance
          field(cellX(j),cellY(j)) = 0; 
       end
   end
    
   %RecalcDrugs
   deathChance = OrigDeathChance + drugLevel/1000;
   drugLevel = max(drugLevel-drugDecay,0);
   set(drugDisplay, 'string',strcat('Drug Level: ', num2str(drugLevel)));
   set(timeCount, 'string', strcat('Time: ', num2str(i)));
   
   spy(field)
   pause(0.01)
end