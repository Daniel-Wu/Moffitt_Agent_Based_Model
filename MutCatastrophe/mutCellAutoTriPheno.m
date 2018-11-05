%This program will try to model evolutionary dynamics of a constantly
%mutating population. The phenotypes are simply a 2 d gradiant of base
%reproduction rates and drug resistance arranged in a triangle.
%Daniel Wu
%6/19/16
%Qs - how get color of spy dot to reflect phenotype?
clc;clf;clearvars;hold on;
%filename = 'mutCellSimul.gif';


%Constants, all in terms of chance per step.
mutChance = 0.1; %0.001
maxGrowthRate = 0.15;
minGrowthRate = 0.05;
maxDrugResist = 1;      %Porportion drug resisted
fieldSize = 600;
time = 500;
cancerStartSize = 4;
deviation = 5;

%Set up field
field = zeros(fieldSize);

%Set up phenotypes
phenotypeSize = [100 100];
emptyPhenotypes = ones(phenotypeSize(1),phenotypeSize(2));
validPhenotypes = ~rot90(quarterCircleOnes(phenotypeSize(1)),2); %rot90(triu(emptyPhenotypes)); %Circle/triangle
validPhenotypeInd = find(validPhenotypes);
phenotypeCount = emptyPhenotypes-1;
phenotypeNum = length(validPhenotypeInd);

%Here's some vectors that contain probabilities for the ith and jth
%phenotypes
phenotypeGrowth = linspace(minGrowthRate,maxGrowthRate, phenotypeSize(1));
phenotypeDrugResist = linspace(0,maxDrugResist, phenotypeSize(2));
[origGrowthChanceMat, drugResistMat] = meshgrid(phenotypeGrowth, phenotypeDrugResist);
growthChanceMat = origGrowthChanceMat;
origDeathChanceMat = zeros(phenotypeSize(1),phenotypeSize(2)) + 0.3;
deathChanceMat = origDeathChanceMat;        %Modified in code

%Drug function
dose = 500;         %Some sort of arbitrary dose.
drugDecay = 5;      %Some sort of decay, dose/step.
drugLevel = 0;%500;      %Current amount of drugs in system.

%Let's track population of cells over time.
population = NaN(1,ceil(time/5));

%Here's a nice 4x4 group of random starting cells.
middleLoc = round(size(field)/2);
startField = field(middleLoc:middleLoc+cancerStartSize-1,middleLoc:middleLoc+cancerStartSize-1);

%Now we randomly give those starting cells a phenotype and count them.

for i = 1:length(startField)
   for j = 1:length(startField) 
        startField(i,j) = validPhenotypeInd(ceil(rand()*phenotypeNum));
        phenotypeCount(startField(i,j)) = phenotypeCount(startField(i,j))+1;
   end
end
%And then we put it back into the field.
field(middleLoc:middleLoc+cancerStartSize-1,middleLoc:middleLoc+cancerStartSize-1) = startField;

%Now we make a pretty interface
figure(1)
button = uicontrol('string','Hit them with the drugs!','Callback','drugLevel = drugLevel + dose;');
button.Position(2) = 10;                %Make the button stop covering the plot
button.Position(3) = 150;               %Make the button longer

drugDisplay = uicontrol('style','text','Position',[10 60 80 40]);

timeCount = uicontrol('style','text','Position',[10 100 80 20]);

%%
%Now the crazy buggy simulation part
for i = 1:time                          %For the entire simulation
   [cellX, cellY] = find(field);        %Tells me where all the crap is
   randomOrder = randperm(length(cellX));        %A random order
   cellX = cellX(randomOrder);                  %Shuffle the order of cells
   cellY = cellY(randomOrder);
   for j=1:length(cellX)                %For every cell
      neighborhood = field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1));    %Show me stuff around it
      cellType = field(cellX(j),cellY(j));              %And tell me what kind of cell it is. 
       for xTest = 1:3
           for yTest = 1:3              %And for that stuff
            luck = rand(3);
               if neighborhood(xTest,yTest) == 0        %If there's nothing living there
                   if luck(xTest,yTest)<growthChanceMat(cellType)               %And you're lucky
                      neighborhood(xTest,yTest) = cellType;    %Now there's something there.
                      phenotypeCount(cellType) = phenotypeCount(cellType)+1;  %And we count it

                   end
               end
            
           end
       end
       
       field((cellX(j)-1):(cellX(j)+1),(cellY(j)-1):(cellY(j)+1))=neighborhood;     %And we put that new stuff back in
       
       %But the cell might die
       if rand()<deathChanceMat(cellType)
          field(cellX(j),cellY(j)) = 0;
          phenotypeCount(cellType) = phenotypeCount(cellType)-1;  %And we count it
       %And it also might mutate
       else
           if rand()<mutChance
               phenotypeCount(cellType) = phenotypeCount(cellType)-1;  %We count the loss of this cell
               field(cellX(j),cellY(j)) = stdDevMutTri(cellType,phenotypeSize,validPhenotypeInd,deviation); %Mutate the cell
               if field(cellX(j),cellY(j))~=0;                                                  %If the cell didn't die
                    phenotypeCount(field(cellX(j),cellY(j))) = phenotypeCount(field(cellX(j),cellY(j))) + 1; %And count it
               end
           end
       end
   end
   %RecalcDrugs - back to per step
   %Make a matrix of additional drug death rates after resetting
   currDrugDeath = zeros(phenotypeSize) + drugLevel/1000;      %Affects Death rate
   deathChanceMat = origDeathChanceMat + (1-drugResistMat).*currDrugDeath;
   growthChanceMat = origGrowthChanceMat - (1-drugResistMat).*currDrugDeath/50;
   
   drugLevel = max(drugLevel-drugDecay,0);
   set(drugDisplay, 'string',strcat('Drug Level: ', num2str(drugLevel)));
   set(timeCount, 'string', strcat('Time: ', num2str(i)));
   
   %Now we graph some nice stuff
      spy(field);        %Not color Coded
     
     %ColorCoded
%     hold on
%    modField = mod(field,phenotypeSize(2));  %Focuses us on drug resist
%    spy((field~=0)&(modField<phenotypeSize(2)/5),'r');
%    spy((modField<phenotypeSize(2)*2/5)&(modField>phenotypeSize(2)/5),'m');
%    spy((modField<phenotypeSize(2)*3/5)&(modField>phenotypeSize(2)*2/5),'c');
%    spy((modField<phenotypeSize(2)*4/5)&(modField>phenotypeSize(2)*3/5),'b');
%    spy(modField>phenotypeSize(2)*4/5,'k');
%    
     title('I am a tumor');     %Get rid of this for speed
     pause(0.0001);
%   hold off;
   if mod(i,5) ==0      %Let's plot population over time every 5 steps.
       subplot(2,3,6);
       population(i/5) = length(find(field));
       plot(population);
       title('Number of tumor cells');

       
       
        if mod(i,10)==0          %Let's make nice bar graph every 10 steps
            subplot(2,3,3);
            bar3(phenotypeCount);
            title('Current Phenotypes');
            ylabel('Drug Resistance');
            xlabel('Base Growth Rate');
            zlabel('Number of cells');
            subplot(2,3,[1 2 4 5]);
        else
            subplot(2,3,[1 2 4 5]);

        end
   end
   
   % Do I want to .gif this simulation? If so, uncomment these lines and
   % the one at the top that says filename.
%    frame = getframe(1);
%      im = frame2im(frame);
%      [imind,cm] = rgb2ind(im,256);
%      if i == 1;
%          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%      else
%          imwrite(imind,cm,filename,'gif','WriteMode','append');
%      end
   cla;
end

%Nice Histogram
figure(2);
subplot(1,2,1);
bar3(phenotypeCount);
title('Current Phenotypes');
ylabel('Drug Resistance');
xlabel('Base Growth Rate');
zlabel('Number of cells');
%Nice population Graph
subplot(1,2,2)
plot(5:5:time,population);
title('Number of tumor cells');