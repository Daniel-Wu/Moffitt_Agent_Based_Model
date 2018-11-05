function [finalTesto] = diffuseTesto(initTesto)
%diffuseTesto diffuses the testosterone
%   Uses some sort of fancy diffusion equation I don't understand
%testoLoc = find(initialTesto>0);        %Locations of testosterone
%testoLoc = testoLoc(randperm(length(testoLoc)));        %Now in a random order
[height,width] = size(initTesto);
finalTesto = zeros(height,width);
d=1;
dt=0.5;
h=2;

for i = 1:height
    for j = 1:width
        try
        finalTesto(i,j) = initTesto(i,j) + (d*dt)/(h*h)*(sum(sum(initTesto((i-1):(i+1),(j-1):(j+1))))-9*initTesto(i,j));
        catch
        end
    end
end

