function [testoHelpMat] = calcTestoHelpMat(testoMat)
%calcTestoHelpMat Converts adjusted testo to testo's effect on growth
%   Uses a hyperbola style equation to adjust testo-growth-rate - capping
%   it at 0.1.
%   Takes testoMat and returns testoHelpMat.

testoHelpMat = 0.2./(1 + exp(-25.*testoMat)) - 0.1;  %Log Function

end