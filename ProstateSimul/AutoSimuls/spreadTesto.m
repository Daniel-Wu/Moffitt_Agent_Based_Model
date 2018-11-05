function [finalTesto] = spreadTesto(initialTesto,testoDecay)
%spreadTesto Spreads some of the testosterone around
%   Takes a given input matrix landscape of testosterone, and spreads it
%   around.
neighborhoodDiameter = 5;       
radius=floor(neighborhoodDiameter/2);       %Have to adjust this for even diameters

testoLoc = find(initialTesto>0);        %Locations of testosterone
testoLoc = testoLoc(randperm(length(testoLoc)));        %Now in a random order
[height,width] = size(initialTesto);

finalTesto = zeros(size(initialTesto));     %Initialize answer
initialTesto=initialTesto.*(1-testoDecay);      %And decay testosterone

for i = 1:length(testoLoc)              %For each place with testosterone
    location = testoLoc(i);
    [locationI,locationJ] = ind2sub([height, width],location);  %Make them subscript
    onEdge = false;
    try                                 %Try making a neighborhood
        neighborhood = initialTesto((locationI-radius):(locationI+radius),(locationJ-radius):(locationJ+radius));
    catch
        onEdge = true;
    end
    
    if ~onEdge
        spreadLoc = neighborhood<neighborhood(radius+1,radius+1);
        spreadAmt = mean([neighborhood(spreadLoc); neighborhood(radius+1,radius+1)]);
        neighborhood(spreadLoc) = spreadAmt; 
        neighborhood(radius+1,radius+1) = spreadAmt;
        finalTesto((locationI-radius):(locationI+radius),(locationJ-radius):(locationJ+radius)) = neighborhood;
    end
    
end

