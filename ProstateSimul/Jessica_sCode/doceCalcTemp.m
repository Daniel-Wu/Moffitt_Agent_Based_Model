function [newG] = doceCalcTemp(drugLevel, G)
%doceCalcTemp Returns the new G value
%   Decreases the current G by a function of the current G and the log of 
%   the drug level.
%   Takes double G value
%   Takes double non-neg drug level

effectiveDrug = log((drugLevel/100)+1);
    
newG = G - (effectiveDrug.*(exp(G)))/10;

end

