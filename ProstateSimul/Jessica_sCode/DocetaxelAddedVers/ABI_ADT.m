%% Abiraterone Treatment

% Set growth rates
r = [0.001, 0.001, 0.001];

% Set carrying capacities
k_max = [0, 200, 200];

% Abiraterone affects the level that T+ can benefit from TP

% Set max symbiotic level
X = 500;

doceTime = time + doceDisplacement;

AbiLength = 100;
createDownSlope = X*(exp(0.05 * (1-(1:AbiLength))));

X(time:time+50000) = 0;
X(time:time + AbiLength - 1) = createDownSlope';

createDownSlope = 200*(exp(0.05 * (1-(1:AbiLength))));
k_TP(time:time + AbiLength - 1) = createDownSlope';
k_TP(time + AbiLength:time + 50000) = createDownSlope(end);

% Assign values for current order index
for i = 1:1:6
    assignin('base', cell2mat(params(1,i)), values(i));
end

%Create matrix
payoffMatrix = [0 a b; c 0 d; e f 0];

IncomingPSA = PSA;

while (PSA > PSA_StopABI)
    if doceTime == time
       doceLevel = doceLevel + doceDose; 
    end
    
    time = time + 1;
    
    % TreatmentIndex = 3 for Abi
    treatmentIndex(time) = 3;
    
    p = x./sum(x);
    
    %Don't we hit T+ twice as hard here?
    k_max(1) = p(2) * X(time);
    k_max(2) = k_TP(time);
    
    % Get dynamics
    for i = 1:1:3
        
        E(i) = payoffMatrix(i,:) * p';
        
        G(i) = r(i) * (k_max(i) - (1-E(i)) * sum(x))/k_max(i);
        
        if (G(i) == -inf)
            G(i) = -1;
        end
        G(i) = doceCalcTemp(doceLevel,G(i));    %Apply Doce
        delta_x(i) = x(i) * G(i);
        
    end
    
    doceLevel = max(doceLevel*(1-doceDecay), 0);    %Decay Doce
    
    x = x + delta_x;
    
    if(x(1) <= 0.1)
        x(1) = 0.1;
        %disp(['T+ Saved during Abi at time :' num2str(time)])

    end
    
    if(x(2) <= 0.1)
        x(2) = 0.1;
        %disp(['Tp Saved during Abi at time :' num2str(time)])

    end
    
    if(x(3) <= 0.1)
        x(3) = 0.1;
        disp(['T- Saved during Abi at time :' num2str(time)])

    end
    
    delta_PSA = getDeltaPSA(x, p);
    PSA = PSA + delta_PSA - sigmaPSA*PSA;
    
    
    all_x(time,:) = x;
    all_p(time, :) = p;
    all_PSA(time) = PSA;
    all_k(time,:) = k_max;

    % Enter abiraterone failure catch. If we are in constant ABI then just stop at 140.
    % If we are in cycles and don't see a decrease in PSA from ABI then break out.
    
    
    if (PSA >= IncomingPSA * 1.1)
        disp('Abiraterone Failure 10% greater than incoming')
        nextCycleFlag = 0;
        break;
    else
        nextCycleFlag = 1;
    end
    
    
end
