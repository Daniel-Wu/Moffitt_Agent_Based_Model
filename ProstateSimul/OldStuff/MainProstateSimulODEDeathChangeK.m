%MainProstateSimulODE

%This is my shot at making a prostate cancer model, with certain
%assumptions listed in the accompanying text file. Three populations,
%Producers, moochers, and independants interact, and we look at how
%medication hits them. This is the ODE version.



clear all;clc;clf;

%Simulation stuffs:
time= 20000;           %Steps
doseSchedule = [1000:1000:time-1000];

%Coefficients
r_tP = 0.06;         %absolute growth rates (for tP)
r_tM = 0.1;
r_tI = 0.05;
testosteroneHelp = .001;       %How much testosterone impacts tM (really same as r_tM)

d_tP = 0.02;         %Natural/normal death rates (for tP)
d_tM =  0.02;
d_tI = 0.02;

drugMurder = .01;     %Drug deathliness constant (cells/time)
drugDecay  = 0.0015;     %How much drug concentration decreases per day (1 is 100%).
dose = 1000;            %DrugDose

testosteroneProduction = .1;                 %Units Per tP
testosteroneDecay = .2;            %1 is 100% %Test/day
testosteroneConsumption = 0.005;             %Is this a thing? UnitsTest/cell/day

K = 10000;          %Carrying capacity (cells)
Kgrowth = .5;       %Change in k per step

%Initialize cells
tP0 = 500;           %Number of inital Producer Cells
tM0 = 1000;          %Number of Moochers
tI0 = 20;            %Number of Independants

%Actual program here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testosteroneInit = 500;
drugInit = 0;

%pop change for tP is: 
%r_tP*tP(i-1)*log(K/totalCells(i-1))                     %Gompertz growth
% - drugMurder*tP*drug(i-1)                             %Drug murder
%

%I included a testosterone factor for the moochers based on porportion of
%producers

%x is my vector:
%[Tp, T+, T-, total cells, testosterone, drugs,K]

func = @(t,x)[
    r_tP * x(1) * log(x(6)/(x(1) + x(2) + x(3))) - d_tP*x(1)- drugMurder*x(5);
    r_tM * x(2) * log(x(6)/(x(1) + x(2) + x(3))) * x(4) * testosteroneHelp - d_tM*x(2)-drugMurder*x(5);
    r_tI * x(3) * log(x(6)/(x(1) + x(2) + x(3))) - d_tI*x(3);
     - testosteroneDecay * x(4) + testosteroneProduction * x(1) - testosteroneConsumption * (x(1) + x(2));
    -drugDecay*x(5);
    Kgrowth;
    
];

options = odeset('NonNegative',1:5);        %Makes sure all quantities are non-negative
%[t, data] = ode45(func,[0 time],[tP0, tM0, tI0, testosteroneInit, drugInit],options);

[t, data] = ode45(func,[0 doseSchedule(1)],[tP0, tM0, tI0, testosteroneInit, drugInit,K],options);
doseSchedule = [doseSchedule,time];
for i = 2:length(doseSchedule)      %Generates data
data(end,5) = data(end,5)+dose;
[newT,newData] = ode45(func,[doseSchedule(i-1), doseSchedule(i)],data(end,:),options);

t = [t;newT];
data = [data;newData];
end

%for i = 2:stepNum+1
%    tP(i) = max(    tP(i-1) + step * (r_tP * tP(i-1) * log(K/totalCells(i-1)))     ,0);
%    tM(i) = max(    tM(i-1) + step * (r_tM * tM(i-1) * log(K/totalCells(i-1)) * testosterone(i-1) * testosteroneHelp)     ,0);
%    tI(i) = max(    tI(i-1) + step * (r_tI * tI(i-1) * log(K/totalCells(i-1)))     ,0);
%    totalCells(i) = tP(i)+tM(i)+tI(i);             %Gives me total cancer
% 
%
%    drug(i) = max(  drug(i-1) + step * (- drugDecay)   ,0);        %Decays me some drugs
%    testosterone(i) = max(  testosterone(i-1) + step * ((1 - testosteroneDecay) * testosterone(i - 1) + testosteroneProduction * tP(i - 1) - testosteroneConsumption * (tM(i-1) + tP(i-1)))   ,0);
%end


%Plot totalCells, tP, tM, tI here, maybe drug
hold on
plot(t,data(:,1),'b','linewidth',2)
plot(t,data(:,2),'k','linewidth',2)
plot(t,data(:,3),'r','linewidth',2)
plot(t,data(:,1)+data(:,2) + data(:,3),'c')
plot(t,data(:,4),'g')
plot(t,data(:,5),'m')
legend('Tp','T+','T-','Total Cells','testosterone','drug')