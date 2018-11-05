%My shot at the discrete function - models prostate dynamics system
%I tried making an explicit abiraterone quantity
clc; clear all;

%% Var Init
%Simul Params
endTime = 20000;
time = 0;

%Set coefficients here
a = 0.5;
b = 0.99;
c = 0.01;
d = 0.04;
e = 0.03;
f = 0.02;
preCoeffMat = [0 a b; c 0 d; e f 0];   %Nice matrix of coefficients

%Post Luperon
%Jessica's worst case
a = 0.91;
b = 0.01;
c = 0.96;
d = 0.02;
e = 0.95;
f = 0.9;
postCoeffMat = [0 a b; c 0 d; e f 0];   %Nice matrix of coefficients

%Values of R and K
rPreLup = [0.005, 0.002, 0.001];
kInitPreLup = [600 200 100];
kPreLup = kInitPreLup;
rPreAbi  = [0.005 0.002 0.001];
rPostAbi = [0.001 0.001 0.001];
kInitPostLup = [0 200 200];
kPostLup = kInitPostLup;

%Population init
testo = 0;
pop = [0.1 0.1 0.1];           %[T+ Tp T-]
PSA = 0;
testoPSAImpactPreLup = [0.1 0.8 0.8];
testoPSAImpactPostLup = [0.01 0.006, 0.006];
basePSAProdPreLup = [0 0 0];
basePSAProdPostLup = [0 0.5 0.5];
PSADecay = 0.3;

%Drug variables
abiLevel = 0;
abiDose = 9.9;
doceLevel = 0;
doceDose = 300000;
drugDecay = 0.05;

%Data recording
popData = nan(endTime,3);
popData(1,:) = pop;

PSAData = nan(endTime,1);
cycleCount = 0;

popPorp = pop./sum(pop);
popPorpData = nan(endTime,3);
popPorpData(1,:) = popPorp;

abiTracker = nan(endTime,1);
doceTracker = nan(endTime,1);

%PSA Brightlines
PSA_GiveADT = 80;
PSA_GiveAbi = 300;
PSA_StopAbi = PSA_GiveAbi*0.999;


%% Simul
%Carcinogenesis
    while PSA<PSA_GiveADT           %While we are off Luperon
        
        time = time + 1;            %Increment time
        
        kPreLup(1) = pop(2)*500;          %Change K of T+ based on Tp pop
        
        [pop, popPorp] = growCancer(pop, rPreLup, kPreLup, preCoeffMat);
        PSA = PSA*(1-PSADecay) + sum(pop.*testoPSAImpactPreLup);          %Not t+ dependant nvm calcPSA(pop, popPorp, PSA, PSADecay, testoPSAImpactPreLup, basePSAProdPreLup);

        %Decay drugs
        doceLevel = max(doceLevel*(1-drugDecay), 0);
        abiLevel = max(abiLevel*(1-drugDecay), 0);
        
        %Record Data
        popData(time, :) = pop;
        popPorpData(time, :) = popPorp;
        PSAData(time) = PSA;
        doceTracker(time) = doceLevel;
        abiTracker(time) = abiLevel;
        
    end     


while true
    
    if time>=endTime        %Check to see if we should stop the simul
        break
    end
    
    disp(['Cycle ' num2str(cycleCount) ' complete.'])
    cycleCount = cycleCount + 1; %Count cycles
    
    
    while PSA<PSA_GiveAbi           %While we are off abi
        
        time = time + 1;            %Increment time
        
        if time>=endTime        %Check to see if we should stop the simul
            break
        end
        
        kPostLup(1) = popPorp(2)*500;          %Change K of T+ based on Tp pop
        kPostLup(2) = kInitPostLup(2) - abiLevel;
        if any(kPostLup<0)     %Catch K in case it becomes nonnegative
            kPostLup(kPostLup<0) = 0;
            disp(['k went negative at time ' num2str(time)])
        end
        
        [pop, popPorp] = growCancer(pop, rPreAbi, kPostLup, postCoeffMat);
        PSA = calcPSA(pop, popPorp, PSA, PSADecay, testoPSAImpactPostLup, basePSAProdPostLup);
        
        %Decay drugs
        doceLevel = max(doceLevel*(1-drugDecay), 0);
        abiLevel = max(abiLevel*(1-drugDecay), 0);
        
        %Record Data
        popData(time, :) = pop;
        popPorpData(time, :) = popPorp;
        PSAData(time) = PSA;
        doceTracker(time) = doceLevel;
        abiTracker(time) = abiLevel;
        %% Debug
        kTracker(time,:) = kPostLup;
    end     %End of vacation
    
    
    
    while PSA>PSA_StopAbi           %Until we can reach that brightline, we treat
        
        time = time + 1;            %Increment time
        
        if time>=endTime        %Check to see if we should stop the simul
            break
        end
        
        %Continuously add in abi
        abiLevel = abiLevel + abiDose;
        
        %I don't agree with this K - we implement abi 3 ways? new r, new k,
        %new X? 
        %TI population does the same thing regardless of T+ and Tp it
        %seems
        kPostLup(1) = popPorp(2)*(500-abiLevel); %Constant X?
        kPostLup(2) = kInitPostLup(2) - abiLevel;  %Linear relationship? Graph is actually exponential
        if any(kPostLup<0)     %Catch K in case it becomes nonnegative
            kPostLup(kPostLup<0) = 0;
            disp(['k went negative at time ' num2str(time)])
        end
        
        [pop, popPorp] = growCancer(pop, rPostAbi, kPostLup, postCoeffMat);
        PSA = calcPSA(pop, popPorp, PSA, PSADecay, testoPSAImpactPostLup, basePSAProdPostLup);
        
        %Decay drugs
        doceLevel = max(doceLevel*(1-drugDecay), 0);
        abiLevel = max(abiLevel*(1-drugDecay), 0);
        
        %Record Data
        popData(time, :) = pop;
        popPorpData(time, :) = popPorp;
        PSAData(time) = PSA;
        doceTracker(time) = doceLevel;
        abiTracker(time) = abiLevel;
        %% Debug
        kTracker(time,:) = kPostLup;
    end     %end of treatment
    
end

%% Plot Results

figure
subplot(3,1,1) %Populations
hold on;
xlabel('Time')
ylabel('Population')
plot(popData,'lineWidth',3)

subplot(3,1,2) %Porportions
hold on;
xlabel('Time')
ylabel('Population porportions')
plot(popPorpData,'lineWidth',3)

subplot(3,1,3) %PSA
hold on;
xlabel('Time')
ylabel('PSA Level')
plot(PSAData)

