function td = init_timingDat(scr)
td = struct;
td.phase = {};
td.repN = [];
td.trialN = [];
td.trialType = [];

td.cv_ivi = [];
td.mean_ivi = [];

a_ph = fields(scr);
nph = length(a_ph);
nt = 0;

for i1 = 1 : nph
    ph = a_ph{i1};
    
    for i2 = 1 : scr.(ph).nReps
        rep = sprintf('rep%d', i2);
        t_nt = length(scr.(ph).(rep).trialOrder);
        idx_t = 1 : t_nt;
        
%         if length(ph) > 5 && (isequal(ph(1 : 5), 'pract') || isequal(ph(1 : 5), 'inter'))
%             idx_t = find(scr.(ph).(rep).trialOrder <= 2); % No baseline trials for the pract* and inter* phases
%             t_nt = length(idx_t);
%         end
        
        td.phase = [td.phase, repmat({ph}, 1, t_nt)];
        td.repN = [td.repN, repmat(i2, 1, t_nt)];
        td.trialN = [td.trialN, 1 : t_nt];
        td.trialType = [td.trialType, scr.(ph).(rep).trialOrder(idx_t)];
        
        td.cv_ivi = [td.cv_ivi, nan(1, t_nt)];
        td.mean_ivi = [td.mean_ivi, nan(1, t_nt)];
    end
end

td.trialCnt = 1;
return
