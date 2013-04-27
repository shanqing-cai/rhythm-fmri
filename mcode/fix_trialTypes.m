function fix_trialTypes()
load '../dacache/PILOT_ANS_F01.mat'; % gives pdata
load 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\expt.mat' % gives expt

%%
for i1 = 1 : length(pdata.mainData.trialNums)
    t_fn = pdata.mainData.rawDataFNs(i1);
    
    t_phase = pdata.mainData.phases{i1};
    t_repNum = pdata.mainData.blockNums(i1);
    t_trialNum = pdata.mainData.trialNums(i1);
    
    t_rep = sprintf('rep%d', t_repNum);
	
    t_trialType = expt.script.(t_phase).(t_rep).pertType(t_trialNum);
    pdata.mainData.pertType(i1) = t_trialType;
end

%%
save('../dacache/PILOT_ANS_F01.mat', 'pdata');
return