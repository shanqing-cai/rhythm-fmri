function n = get_trial_index(scr, phase, rn, tn)
allPhases = fields(scr);

if isempty(fsic(allPhases, phase))
    error('Cannot find phase %s in script', phase);
end

n = 1;
bFound = 0;
for i1 = 1 : numel(allPhases)
    ph = allPhases{i1};
    for i2 = 1 : scr.(ph).nReps
        rep = sprintf('rep%d', i2);
        
        for i3 = 1 : length(scr.(ph).(rep).trialOrder)
            if isequal(phase, ph) && rn == i2 && tn ==i3;
                bFound = 1;
                return;
            end
            
            n = n + 1;
        end
        
    end
    
end

if bFound == 0
    error('Cannot find phase %s, rep #%d, trial #%d in script', ...
          phase, rn, tn);
end
return