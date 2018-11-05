

%% These are all 22 Post-ADT matrix orders
% Currently unused.

% testOrders = {...
%     'f','e','b','d','c','a';...
%     'f','e','b','d','a','c';...
%     'f','e','b','a','d','c';...
%     'f','b','d','e','c','a';...
%     'f','b','d','e','a','c';...
%     'f','b','d','a','e','c';...
%     'f','b','e','d','c','a';...
%     'f','b','e','d','a','c';...
%     'f','b','e','a','d','c';...
%     'f','b','a','e','d','c';...
%     'f','b','a','d','e','c';...
%     'b','d','f','e','c','a';...
%     'b','d','f','e','a','c';...
%     'b','d','f','a','e','c';...
%     'b','f','d','e','c','a';...
%     'b','f','d','e','a','c';...
%     'b','f','d','a','e','c';...
%     'b','f','e','d','c','a';...
%     'b','f','e','d','a','c';...
%     'b','f','e','a','d','c';...
%     'b','f','a','e','d','c';...
%     'b','f','a','d','e','c';...
%     };
% params = testOrders(21,:);
%clear all;
% Set PSA decay rate
sigmaPSA = 0.3;

tempTimer = 0;

%0 = Constant Abi. 1 = Cycling.
constantAbi = 0;

if ~exist('StopAbiFrac','var')
    StopAbiFrac = 1/3; %Fraction at which to stop Abi
end

doceDisplacement = 0; %How many time steps after abi to give docetaxel (>=0)

if ~exist('doceDose','var')
    doceDose = 10;
end
doceLevel = 0;
doceDecay = 0.02;

% Naive Pre-Lupron
% values = [0.0100    0.2060    0.4020    0.5980    0.7940    0.99];

% Best Case 21 Post
testCase_values(1, :) = [0.0100    0.02    0.03    0.4    0.95    0.99];
testCase_params(1, :) = {'b','f','a','e','d','c'};


% Worst Case 14 Post
testCase_values(2,:) = [0.0100    0.02    0.9    0.91    0.95    0.96];
testCase_params(2,:) = {'b','d','f','a','e','c'};


% Middle ground 1 equal Post
testCase_values(3,:) = [0.994    0.995   0.996   0.997    0.998   0.999];
testCase_params(3,:) = {'b','f','d','a','e','c'};




for TestCaseIndex = 2:2
    
    disp(TestCaseIndex)
    values = testCase_values(TestCaseIndex, :);
    params = testCase_params(TestCaseIndex, :);
    
    
    %% ADT ONLY
    %     PSA_GiveADT = 80;
    %     PSA_GiveABI = 600;
    %
    %     NAIVE;
    %
    %     ADT_ONLY;
    %
    %     SimulationEndFunction;
    
    
    
    if (constantAbi == 1)
        
        % Constant ABI
        PSA_GiveADT = 80;
        PSA_GiveABI = 300;
        PSA_StopABI = 0;
        
        NAIVE;
        
        ADT_ONLY;
        ABI_ADT;
        
        SimulationEndFunction;
        
        
    else
        % Cycle ABI
        PSA_GiveADT = 80;
        PSA_GiveABI = 300;
        PSA_StopABI = PSA_GiveABI*StopAbiFrac;
        
        NAIVE;
        
        % Always give one cycle of ABI
        cycleNum = 0;
        nextCycleFlag = 1;
        
        while (nextCycleFlag == 1)
            cycleNum = cycleNum + 1;
            %fprintf('Cycle %d\n', cycleNum)
            ADT_ONLY;
            ABI_ADT;
        end
        
        SimulationEndFunction;
        ExtractAbiCycleDynamics;
        
        
    end
    
end


