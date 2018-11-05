%Views simul data - Import data fieldData, testoFieldData, TIPop, TMPop,
%TPPop, totalPop, testoTracker, drugTracker.
close all;


showSimul = 0;
makeVideo = 0;
vidName = 'NoTreat3';

%Greyscale Version
set(groot,'DefaultAxesColorOrder',[[1 0 0];[1 0.2 0.2];[1 0.4 0.4];[1 0.6 0.6];[1 0.8 0.8];[1 0.9 0.9]])
%Color Version
set(groot,'DefaultAxesColorOrder','remove')



%Make simul pretty
if(size(fieldData,1)==100)
    isSmallField = true;
    dotSize = 25;
else
    isSmallField = false;
end

time = length(TPPop);

totalPop = TMPop + TIPop + TPPop;
%Overview graph
if makeVideo
    overview = figure(1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end

hold on;
plot(TMPop,'--','lineWidth',5);
plot(TIPop,'-.','lineWidth',5);

if isSmallField
    plot(testoTracker/20,'lineWidth',2)   %Scale to fit.
else
    plot(200:time,testoTracker(200:end)/20, 'lineWidth',2) %cuts off beginning, incase graph scales strangely.
end

plot(totalPop,'lineWidth',5);
plot(TPPop,':','lineWidth',5);
%plot(drugTracker, '--','lineWidth', 1)
title('Cancer Population')
xlabel('Time')
ylabel('Cell Population')

try
    %plot(doceTracker,':', 'lineWidth', 2)
    legend('Consumers','Independents','Testosterone', 'Total Cancer Population','Producers')%,'Abi Level','Docetaxel Level')
catch
    legend('Consumers','Independents','Testosterone','Total Cancer Population','Producers')%,'Abi Level')
end


%Draw lines for Drugging times
ylimits = ylim;
if exist('AbiSchedule', 'var')
    if length(AbiSchedule) > 0
        for i = 1:length(AbiSchedule)
            linehand = line([AbiSchedule(i) AbiSchedule(i)],[0 ylimits(2)],'lineWidth',3,'color',[0.8 0.8 0.8]);
            uistack(linehand,'bottom')
        end
    end
end

if exist('doceSchedule', 'var')
    if length(doceSchedule) > 0
        for i = 1:length(doceSchedule)
            line([doceSchedule(i) doceSchedule(i)],[0 ylimits(2)],'lineWidth',3,'color','k')
        end
    end
end

if exist('testoSchedule', 'var')
    if length(testoSchedule) > 0
        for i = 1:length(testoSchedule)
            line([testoSchedule(i) testoSchedule(i)], [0 ylimits(2)], 'lineWidth', 3, 'color', 'c')
        end
    end
end

%% Field simulation
if showSimul
    pause(5)    
    myFigureThing = figure(2);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    if makeVideo
        v = VideoWriter(['C:\Users\User\Desktop\simulDatas\' vidName '.avi']);
        v.FrameRate = 150;
        open(v);
    end
    
    timeDisplay = uicontrol('style','text','Position',[20 100 80 20]);
    
    if isSmallField     %Dot size is defined
        for i=1:time
           set(timeDisplay, 'string', strcat('Time: ', num2str(i)));
           subplot(2,3,[1 2 4 5])
           cla;
           hold on
           spy(fieldData(:,:,i) == 1, 'k', dotSize-5)
           spy(fieldData(:,:,i) == 2, 'r', dotSize)
           spy(fieldData(:,:,i) == 3, 'g', dotSize)
           hold off
           subplot(2,3,3)
           
           %surf(testoFieldData(:,:,i))
           pcolor(flipud(testoFieldData(:,:,i)))  %Flat Display
           
           shading interp;
           drawnow;
           if makeVideo
               frame = getframe;
               writeVideo(v,frame);
           end
        end
        
        
    else              %Let matlab do auto dot size
        for i=1:time
           set(timeDisplay, 'string', strcat('Time: ', num2str(i)));
           subplot(2,3,[1 2 4 5])
           cla;
           hold on
           spy(fieldData(:,:,i) == 1, 'k')
           spy(fieldData(:,:,i) == 2, 'r')
           spy(fieldData(:,:,i) == 3, 'g')
           xlabel('')
           title('Field')
           hold off
           
           subplot(2,3,3)
           %surf(testoFieldData(:,:,i))
           pcolor(flipud(testoFieldData(:,:,i)))    %Flat Display
           title('Testosterone Field')
           shading interp;
           drawnow;
           
           
           
           if makeVideo
               drawnow;
               frame = getframe(myFigureThing);
               pause(0.01)
               writeVideo(v,frame);
               pause(0.01)
           end
        end
    end
end

if makeVideo
    frame = getframe(overview);
    pause(0.01)
    writeVideo(v,frame);
    pause(0.01)
    close(v);
end