function [mutPhenotype] = stdDevMut(phenotype, phenMatDimensions, phenotypeNum, deviation)
%stdDevMut Calculates and returns a phenotype value for mutated cell. 
%   With input arguments int phenotype, the matrix size of phenotype
%   matrix, and the mutation deviation, gives mutated phenotype.

%Does rounding this screw up the whole stdDev thing?
xChange = round(deviation*randn());
yChange = round(deviation*randn());


mutPhenotype = phenotype + xChange*phenMatDimensions(1) + yChange;

%Bad solution to wraparound issue - mutate randomly way off
if((mutPhenotype>phenotypeNum)||(mutPhenotype<0))
    mutPhenotype = mod(mutPhenotype,phenotypeNum);
else if mutPhenotype==0
    mutPhenotype = 1;
end
%Right now there is a crash if the mut phenotype lands exactly on 0.

%Another bad solution - assume that going out means that it's dead.


%Need to add check if out of phenotype matrix. Perhaps do while, or reverse
%mut direction? I should also change this to accept matrix of input
%phenotypes to make crazy fast.
 
end