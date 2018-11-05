function [] = cheatSave(folderName, fileName, simulCount, TMPop, TIPop, TPPop, testoTracker, drugTracker, doceTracker, PSATracker, AbiSchedule, testoSchedule, doceSchedule, param)
%cheatSave Saves stuff
%   Saves input args

    save([folderName '/' fileName num2str(simulCount)],'TMPop','TIPop','TPPop','testoTracker','drugTracker','doceTracker','PSATracker', 'AbiSchedule', 'testoSchedule','doceSchedule','param','-v7.3')

end

