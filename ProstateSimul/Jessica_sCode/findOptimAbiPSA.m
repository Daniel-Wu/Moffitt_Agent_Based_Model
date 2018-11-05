%Finds some nice drug doses

%Threshold for quality of life, in PSA. qual for 100 = qual for 200, but
%quality linearly decays for all PSA's after, down to 0 at 300. Random
%algorithm.


testValues = linspace(0.1,0.9999,9000);
data = zeros(1,length(testValues));
PSAdata = zeros(1,length(testValues));
for iSimuls = 1:length(testValues)
    StopAbiFrac = testValues(iSimuls);
    July5_MainNoOutput;
    data(iSimuls) = abiCycleInfo(end,1);
    PSAdata(iSimuls) = mean(all_PSA);
    disp(StopAbiFrac)
    clearvars -except iSimuls data testValues PSAdata;
    
end
figure(1)
plot(testValues,data);
xlabel('Abiraterone Stop PSA (As a fraction of Abiraterone Start PSA)')
ylabel('Simulation time until progression')
title('Impact of Therapy Strategies on Progression Time')

figure(2)
%Calc quality of life

lifeQualThresh = 200;
qualLife = data.*(1 + max(-1, max(0, PSAdata-lifeQualThresh)*-0.01));
plot(testValues, qualLife)
xlabel('Abiraterone Stop PSA (As a fraction of Abiraterone Start PSA)')
ylabel('Quality of life measure')
title('Impact of Therapy Strategies on Total Quality of Life')
save('bruteForceOptimResults.mat')