%Views simul data - Import data fieldData, testoFieldData, TIPop, TMPop,
%TPPop, totalPop, testoTracker, drugTracker.
close all;
time = length(TPPop);

%Overview graph
figure(1);
hold on;
plot(TMPop);
plot(TIPop,'r');
plot(TPPop,'g');
plot((TMPop+TIPop+TPPop),'m');
plot(testoTracker/20,'c')   %Scale to fit.
plot(drugTracker)
title('Cancer Population')
legend('Consumers','Independents','Producers','Total Cancer Population','Testosterone','Abi Level')

pause(5)

%Field simulation
figure(2);
for i=1:time
   subplot(2,3,[1 2 4 5])
   cla;
   hold on
   spy(fieldData(:,:,i) ==1, 25)
   spy(fieldData(:,:,i) ==2,'r',25)
   spy(fieldData(:,:,i) ==3,'g',25)
   hold off
   subplot(2,3,3)
   surf(testoFieldData(:,:,i))
   pause(0.005)
end