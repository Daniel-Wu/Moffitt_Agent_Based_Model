
abiON(treatmentIndex == 3) = 1;

abiCycle = 0;

for i = 1:1:length(treatmentIndex) - 1;
%     for i = 1:1:100000
    
    if abiON(i + 1) - abiON(i) == 1;
        
        abiCycle = abiCycle + 1;
        
        abiCycleInfo(abiCycle, 1) = i;
        
    elseif abiON(i + 1) - abiON(i) == -1;
        
        abiCycleInfo(abiCycle, 2) = i - abiCycleInfo(abiCycle, 1);
        
    end
    
end
disp (abiCycleInfo)




