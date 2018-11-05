%% Pre_ADT
% All three species will be creating their highest level of PSA.

% Set growth rates
r = [0.005, 0.002, 0.001];

% Set carrying capacities
k_max = [600, 200, 100];

% Set Pre-ADT matrix parameters (This is from T+ maximization analysis)
c = 0.01;
f = 0.02;
e = 0.03;
d = 0.04;
a = 0.5;
b = 0.99;

%Create matrix
payoffMatrix = [0 a b; c 0 d; e f 0];

% Set all cells low to simulate carcinogensis
x = [0.01 0.01, 0.01];

% Set initial pre-cancer PSA level.
PSA = 0;

time = 0;

while (PSA < PSA_GiveADT)
    time = time + 1;
    
    % TreatmentIndex = 1 for Naive
    treatmentIndex(time) = 1;
    
    
    p = x./sum(x);
    
    % Get dynamics
    for i = 1:1:3
        
        E(i) = payoffMatrix(i,:) * p';
        
        G(i) = r(i) * (k_max(i) - (1-E(i)) * sum(x))/k_max(i);
        
        if (G(i) == -inf)
            G(i) = -1;
        end
        
        delta_x(i) = x(i) * G(i);
        
    end
    
    x = x + delta_x;
    
    if(x(1) <= 0.1)
        x(1) = 0.1;
    end
    
    if(x(2) <= 0.1)
        x(2) = 0.1;
    end
    
    if(x(3) <= 0.1)
        x(3) = 0.1;
    end
    
    delta_PSA = 0.1 * x(1) + .8 * x(2) + .8 * x(3) - sigmaPSA*PSA;
    PSA = PSA + delta_PSA;
    all_x(time,:) = x;
    all_p(time, :) = p;
    all_PSA(time) = PSA;
    
    
end