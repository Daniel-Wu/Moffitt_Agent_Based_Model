function [] = simulAggregate(path, simulNum)
%simulAggregate Aggregates and displays data from lots of simulations
%   Lots of data will be combined and displayed. 
%   Args must be a String path(w/o file number but with file name)
%   and a matrix of fileNumbers

%% Init
clf; close all; clc;

simulCategory = '';      %Will hold the category of the simul

%Simul catagories
simulCatNames = {'totalExtinctionSimuls', 'earlyTIDeathSimuls', 'allOtherSimuls', 'TIDominantSimuls'};
totalExtinctionSimuls.num = [];

earlyTIDeathSimuls.num = [];

allOtherSimuls.num = [];

TIDominantSimuls.num = [];

%Other Data
TIEndPorpData = zeros(1,length(simulNum));
TIPorpData = zeros(1, 1500);
%% Aggregation
for iSimul = 1:length(simulNum)           %For each simul
    
    load([path num2str(simulNum(iSimul)) '.mat'])    %Get data

    TIEndPorp = TIPop(end)/(TIPop(end) + TPPop(end) + TMPop(end)); %Porportion of T- at end
    
    %Categorize Simuls
    if TIPop(899) < 1               %If TIPop dies kinda early
        simulCategory = 'earlyTIDeathSimuls';
        
    else if (TIEndPorp) > 0.3   %If TI is more than 30% of pop at end
            simulCategory = 'TIDominantSimuls';
            
        else if (TIPop(end) + TPPop(end) + TMPop(end)) == 0     %If everyone is dead
                simulCategory = 'totalExtinctionSimuls';
            else
                simulCategory = 'allOtherSimuls';                 
            end 
        end
    end
    
    %Remember simul num
    eval([simulCategory '.num = [' simulCategory '.num simulNum(iSimul)];'])
    
    %Record other data
    TIEndPorpData(iSimul) = TIEndPorp;
    TIPorpData = TIPorpData + TIPop./(TIPop + TMPop + TPPop);
end                                     %End of for each simul

%Do stuff to the data
TIPorpData = TIPorpData/length(simulNum);
eval(['clearvars -except path simulCatNames simulNum TIEndPorpData TIPorpData ' strjoin(simulCatNames)]) %Get rid of extra vars

save([path 'Info.mat'])   %Save companion file
%% Results display
for iDisp = 1:length(simulCatNames)   %Disp num of simuls in each cat
    disp(simulCatNames(iDisp))
    eval([ 'disp(length(' cell2mat(simulCatNames(iDisp)) '.num))' ])
end

figure(2)
title('T- Final Porportions')
xlabel('Porportion of population')
histogram(TIEndPorpData,10)
end