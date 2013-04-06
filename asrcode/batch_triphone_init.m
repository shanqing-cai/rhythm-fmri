function batch_triphone_init
%% CONFIG
% BASE_MAT = 'E:\DATA\RHYTHM-FMRI\TestExpt_audapter_1\pre\rep1\trial-1-2.mat';
BASE_MAT = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_1\pre\rep1\trial-1-2.mat';
%%

stcs = sents();

ts = nan(1, numel(stcs));
for i1 = 1 : numel(stcs)
    sent = stcs{i1};
    
    fprintf('Working on sentence %d of %d: %s\n', i1, numel(stcs), sent);
    tic;
    pa = run_julian(BASE_MAT, 'sent', sent);
    ts(i1) = toc;
end
return