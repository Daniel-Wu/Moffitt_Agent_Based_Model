%This is my shot at making a prostate cancer model, with certain
%assumptions listed in the accompanying text file. Three populations,
%Producers, moochers, and independants interact, and we look at how
%medication hits them.
clearvars;clc;clf;

%Simulation stuffs:
step = .001;           %Step time (days)
time = 3%65;         %Simulation time (days)

stepNum = round(time/step);     %integer number of steps.  
t = 0:step:time;    %time vector

%Coefficients
r_tP = 0.6;         %growth rates (for tP)
r_tM = 1;
r_tI = 0.5;
testosteroneHelp = 0.00005;       %How much testosterone impacts tM (really same as r_tM)

drugMurder = 1;     %Drug deathliness constant
drugDecayUnadjusted  = 1;     %How much drug concentration decreases per day (1 is 100%).
drugDecay = step * drugDecayUnadjusted;       %drug decay per step.

testosteroneProduction = 2;                 %Units Per tP
testosteroneDecayUnadjusted = 1;            %1 is 100%
testosteroneDecay = step * testosteroneDecayUnadjusted;     %testosterone decay per step

K = 10000;          %Carrying capacity (cells)

%Initialize cells
tP0 = 200;           %Number of inital Producer Cells
tM0 = 1000;          %Number of Moochers
tI0 = 20;            %Number of Independants


%Actual program here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialize population vectors
tP = zeros(1,stepNum);
tM = zeros(1,stepNum);
tI = zeros(1,stepNum);
totalCells = zeros(1,stepNum);

tP(1) = tP0;
tM(1) = tM0;
tI(1) = tI0;
totalCells(1) = tP0 + tM0 + tI0;

drug = zeros(1,stepNum);                            %DrugConcentrationVector
testosterone = zeros(1,stepNum);                           %Testosterone Concentration
testosterone(1) = 1000;

%pop change for tP is: 
%r_tP*tP(i-1)*log(K/totalCells(i-1))                     %Gompertz growth
% - drugMurder*tP*drug(i-1)                             %Drug murder
%

%I included a testosterone factor for the moochers based on porportion of
%producers

for i = 2:stepNum+1
    tP(i) = max(r_tP*tP(i-1)*log(K/totalCells(i-1)),0);
    tM(i) = max(r_tM*tM(i-1)*log(K/totalCells(i-1))*testosterone(i-1)*testosteroneHelp,0);
    tI(i) = max(r_tI*tI(i-1)*log(K/totalCells(i-1)),0);
    totalCells(i) = tP(i)+tM(i)+tI(i);             %Gives me total cancer pop

    drug(i) = max(drug(i-1) - drugDecay,0);        %Decays me some drugs
    testosterone(i) = max((1-testosteroneDecay)*testosterone(i-1) + testosteroneProduction*tP(i-1),0);
end


%Plot totalCells, tP, tM, tI here, maybe drug
hold on
plot(t,tP,'b')
plot(t,tM,'k')
plot(t,tI,'r')
legend('Tp','T+','T-')