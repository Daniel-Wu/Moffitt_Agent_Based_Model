function [mutPhenotype] = stdDevMutTri(phenotype, phenMatDimensions, validPhenotypes, deviation)
%stdDevMut Calculates and returns a phenotype value for mutated cell. 
%   With input arguments int phenotype, the matrix size of phenotype
%   matrix, and the mutation deviation, gives mutated phenotype, in a
%   triangular matrix of phenotypes.

%Does rounding this screw up the whole stdDev thing?
xChange = round(deviation*randn());
yChange = round(deviation*randn());


mutPhenotype = phenotype + xChange*phenMatDimensions(1) + yChange;

%Bad solution: if new phenotype not on triangle, die.
if ~any(mutPhenotype==validPhenotypes)
    mutPhenotype = 0;
end

end