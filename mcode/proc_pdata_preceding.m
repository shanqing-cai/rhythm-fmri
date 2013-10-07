function [pdata1, idx_noPert_precNoPert, idx_noPert_precF1Up, idx_noPert_precDecel] = proc_pdata_preceding(pdata, expt)
%% proc_pdata_preceding

%% Organize expt script
a_runns = [];
a_repns = [];
a_trialns = [];
a_pertTypes = [];

runn = 1;
while (isfield(expt.script, sprintf('run%d', runn)))
    run = sprintf('run%d', runn);
    
    for i1 = 1 : expt.script.(run).nReps
        rep = sprintf('rep%d', i1);
        
        a_pertTypes = [a_pertTypes, expt.script.(run).(rep).pertType];
        a_runns = [a_runns, runn * ones(size(expt.script.(run).(rep).pertType))];
        a_repns = [a_repns, i1 * ones(size(expt.script.(run).(rep).pertType))];
        a_trialns = [a_trialns, 1 : length(expt.script.(run).(rep).pertType)];
    end
    
    runn = runn + 1;
end


%%
pdata.mainData.precPertType = nan(size(pdata.mainData.pertType));

for i1 = 1 : length(pdata.mainData.precPertType)
    phs = pdata.mainData.phases{i1};
    if ~(length(phs) > 3 && isequal(phs(1 : 3), 'run'))
        continue;
    end
    
    runn = str2double(strrep(phs, 'run', ''));
    repn = pdata.mainData.blockNums(i1);
    trialn = pdata.mainData.trialNums(i1);
    
    idx = find(a_runns == runn & a_repns == repn & a_trialns == trialn);    
    assert(length(idx) == 1);
    assert(pdata.mainData.pertType(i1) == a_pertTypes(idx));
    
    if idx == 1
        continue;
    end
    
    pdata.mainData.precPertType(i1) = a_pertTypes(idx - 1);
end

%% Return values
pdata1 = pdata;

idx_noPert_precNoPert = find(pdata.mainData.pertType == 0 & ...
                             (pdata.mainData.precPertType == 0 | pdata.mainData.precPertType == 4));

idx_noPert_precF1Up = find(pdata.mainData.pertType == 0 & ...
                           pdata.mainData.precPertType == 1);
                       
idx_noPert_precDecel = find(pdata.mainData.pertType == 0 & ...
                            pdata.mainData.precPertType == 2);
return