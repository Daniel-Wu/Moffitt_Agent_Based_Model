%% Post-ADT

% Set growth rates
r = [0.005, 0.002, 0.001];

% Set carrying capacities
k_max = [0, 200, 200];

% Set max symbiotic level
X = 500;

% Assign values for current order index
for i = 1:1:6
    assignin('base', cell2mat(params(1,i)), values(i));
end

%Create matrix
payoffMatrix = [0 a b; c 0 d; e f 0];

startFlag = 1;
while (PSA < PSA_GiveABI || startFlag == 1)
    time = time + 1;
    
    if (PSA < PSA_GiveABI)
        startFlag = 0;
    end
    
    p = x./sum(x);
    
    k_max(1) = p(2) * X;
    
    % Get dynamics
    for i = 1:1:3
        
        E(i) = payoffMatrix(i,:) * p';
        
        G(i) = r(i) * (k_max(i) - (1-E(i)) * sum(x))/k_max(i);
        
        if (G(i) == -inf)
            G(i) = -1;
        end
        if(G(i)==inf)
            disp('AHHHHH Why is G infinity?')
            G(-1) = -1; %Crash the prog
        end
        
        delta_x(i) = x(i) * G(i);
        
    end
    
    doceLevel = max(doceLevel*(1-doceDecay), 0);    %Decay Doce

    x = x + delta_x;
    
    if(x(1) <= 0.1)
        x(1) = 0.1;
        %disp(['T+ Saved during ADT at time :' num2str(time)])
    end
    
    if(x(2) <= 0.1)
        x(2) = 0.1;
        %disp(['Tp Saved during ADT at time :' num2str(time)])
    end
    
    if(x(3) <= 0.1)
        x(3) = 0.1;
        disp(['T- Saved during ADT at time :' num2str(time)])

    end
    
    delta_PSA = getDeltaPSA(x, p);
    PSA = PSA + delta_PSA - sigmaPSA*PSA;

    all_x(time,:) = x;
    all_p(time, :) = p;
    all_PSA(time) = PSA;
    all_k(time,:) = k_max;
    
    
end