%Finds some nice drug doses

testValues = linspace(0,26,1301);
data = zeros(1,length(testValues));

for iSimuls = 1:length(testValues)
    doceDose = testValues(iSimuls);
    July5_MainNoOutput;
    data(iSimuls) = sum(abiCycleInfo(:,1));
    disp(doceDose)
    clearvars -except iSimuls data testValues;
    
end

plot(testValues,data);