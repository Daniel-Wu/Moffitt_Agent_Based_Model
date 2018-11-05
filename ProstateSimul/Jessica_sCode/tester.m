%Tester
%figure(3)
colormap(bone)
drugLevel = 0:999;
G = linspace(-1,1);
[X, Y] = meshgrid(drugLevel, G);
Z = doceCalcTemp(X, Y);
surf(X,Y,Z)
title('Docetaxel effectiveness 2')
xlabel('Docetaxel Level')
ylabel('Original Growth Rate')
zlabel('Adjusted Growth Rate')