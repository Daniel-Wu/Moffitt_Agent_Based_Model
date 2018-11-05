function [newPSA] = calcPSA(pop, popPorp, PSA, PSADecay, testoPSAImpact, basePSAProd)
%calcPSA Calculates and returns PSA
%   Uses made up numbers and values to figure out how much the current PSA
%   is.

newPSA = PSA*(1-PSADecay) + sum(pop.*(testoPSAImpact.*popPorp(2) + basePSAProd));

end

