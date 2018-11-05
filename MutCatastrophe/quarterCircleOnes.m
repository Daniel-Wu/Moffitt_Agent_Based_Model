function [finalMat] = quarterCircleOnes(radius)
%quarterCircleOnes Creates a quarter circle of ones in a matrix of size [radius, radius].
%   The circle is centered on the upper-left corner, and all squares within
%   the circle are set equal to one, and all else are equal to 0.
tic
finalMat = zeros(radius);
for i = 1:radius
    for j = 1:radius
        if (i-1)*(i-1) + (j-1)*(j-1) <= radius*radius;
            finalMat(i,j) = 1;
        end
        
    end
end
toc
end

