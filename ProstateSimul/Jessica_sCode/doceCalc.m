function [newGX] = doceCalc(GX, doceLevel)
%doceCalc Adjusts delta X by killing cells in accordance with docetaxel
%level

%Decided doce impact function for now
%Ideally would like a much smoother function
%when x<-1, y=x-doseadj. When x>-1, y=ln(x+2)-1-doseadj
%Dose adj function is also weird
%Doseadj = ln(dose+1)

GX = linspace(-1,1,100);
docelevel = 0:99;

[XTemp, Y] = meshgrid(GX, docelevel);
Z = zeros(100);

X = XTemp*10;

adjY = log(Y+1)*3;
ind1 = find(X<-1);
ind2 = find(X>=-1);

for i = 1:length(ind1)
    ind = ind1(i);
    Z(ind) = (X(ind) - adjY(ind));
end

for i = 1:length(ind2)
    ind = ind2(i);
    Z(ind) = (1-(adjY(ind)/1000)).*(X(ind) - adjY(ind)) + (adjY(ind)/1000) .* (log(X(ind) + 2) - 1); %- adjY(ind));
end

surf(X, Y, Z);
newGX = 1;
xlabel('Original G')
ylabel('Dox Level')
zlabel('Altered G')
title('Dox Effectiveness')
end

