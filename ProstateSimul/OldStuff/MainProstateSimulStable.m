%This is my shot at making a prostate cancer model, with certain
%assumptions listed in the accompanying text file. Three populations,
%Producers, moochers, and independants interact, and we look at how
%medication hits them.

%Current issues: changing step changes behavior in ways it shouldn't: not
%more detail. I use step to adjust all values so that they are rates in
%terms of step rather than day (cells/day vs cells/step)

clear all;clc;clf;

%Simulation stuffs:
step = 0.5;           %Step time (days)
timeApprox = 365;         %Simulation time (days)
time = (step - mod(timeApprox,step)) + timeApprox;    %Make time divisible by step.

stepNum = time/step;     %integer number of steps.  
t = 0:step:time;    %time vector

%Coefficients
r_tPUnadj = 0.03;         %absolute growth rates (for tP)
r_tMUnadj =  .1;
r_tIUnadj = 0.01;
testosteroneHelpUnadj = 0.0008;       %How much testosterone impacts tM (really same as r_tM)

r_tP = r_tPUnadj*step;
r_tM = r_tMUnadj*step;
r_tI = r_tIUnadj*step;
testosteroneHelp = testosteroneHelpUnadj*step;

drugMurderUnadj = 1;     %Drug deathliness constant (cells/time)
drugMurder = drugMurderUnadj*step;
drugDecayUnadjusted  = 1;     %How much drug concentration decreases per day (1 is 100%).
drugDecay = step * drugDecayUnadjusted;       %drug decay per step.

testosteroneProductionUnadj = 1;                 %Units Per tP
testosteroneProduction = step*testosteroneProductionUnadj;
testosteroneDecayUnadjusted = 3.5;            %1 is 100% %Test/day
testosteroneDecay = step * testosteroneDecayUnadjusted;     %testosterone decay per step
testosteroneConsumptionUnadj = 0.05;             %Is this a thing? UnitsTest/cell/day
testosteroneConsumption = step*testosteroneConsumptionUnadj;

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
testosterone = zeros(1,stepNum);                           %Testosterone Units
testosterone(1) = 100;

%pop change for tP is: 
%r_tP*tP(i-1)*log(K/totalCells(i-1))                     %Gompertz growth
% - drugMurder*tP*drug(i-1)                             %Drug murder
%

%I included a testosterone factor for the moochers based on porportion of
%producers

for i = 2:stepNum+1
    tP(i) = max(    tP(i-1) + step * (r_tP * tP(i-1) * log(K/totalCells(i-1)))     ,0);
    tM(i) = max(    tM(i-1) + step * (r_tM * tM(i-1) * log(K/totalCells(i-1)) * testosterone(i-1) * testosteroneHelp)     ,0);
    tI(i) = max(    tI(i-1) + step * (r_tI * tI(i-1) * log(K/totalCells(i-1)))     ,0);
    totalCells(i) = tP(i)+tM(i)+tI(i);             %Gives me total cancer pop

    drug(i) = max(  drug(i-1) + step * (- drugDecay)   ,0);        %Decays me some drugs
    testosterone(i) = max(  testosterone(i-1) + step * ((1 - testosteroneDecay) * testosterone(i - 1) + testosteroneProduction * tP(i - 1) - testosteroneConsumption * (tM(i-1) + tP(i-1)))   ,0);
end


%Plot totalCells, tP, tM, tI here, maybe drug
hold on
plot(t,tP,'b','linewidth',2)
plot(t,tM,'k','linewidth',2)
plot(t,tI,'r','linewidth',2)
plot(t,totalCells,'y')
plot(t,testosterone,'g')
plot(t,drug,'m')
legend('Tp','T+','T-','Total Cells','testosterone')