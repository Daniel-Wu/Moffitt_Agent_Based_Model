function [newDrugTracker] = applyDoseSchedule(drugTracker, doseSchedule, decay, time)
%applyDoseSchedule Creates a drug tracker in accordance with the dose
%schedule
%applyDoseSchedule(existing mat drugTracker, mat doseSchedule, double decay, int time)
%Returns an updated dose schedule


for iTime = time:(time + length(doseSchedule)-1)
    drugTracker(iTime) = drugTracker(iTime-1)*(1-decay) + doseSchedule(iTime-time+1);
    
end

newDrugTracker = drugTracker;
end

