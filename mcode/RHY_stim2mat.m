function RHY_stim2mat(subjID, TR, TA)
% this script has to be executed inside the subject's directory.
% read subjid from directory name
%% Constants
behavDataDir = 'G:\DATA\RHYTHM-FMRI';

%% Parameters
TR = 11500;    % ms
TA = 2470;    % ms

% DATA_dir = '/users/cais/STUT/DATA';
% funcloc_dir = '/users/cais/STUT/funcloc';

%% Load expt.mat
check_dir(behavDataDir);

subjBehavDataDir = fullfile(behavDataDir, subjID);
check_dir(subjBehavDataDir);

exptMatFN = fullfile(subjBehavDataDir, 'expt.mat');
check_file(exptMatFN);

load(exptMatFN);
assert(exist('expt', 'var') == 1);

%% Determine how many functional runs were administered
nRunsScr = 0;   % Number of runs in the script
while isfield(expt.script, sprintf('run%d', nRunsScr + 1))
    nRunsScr = nRunsScr + 1;
end

nRunsAct = 0;   % Actual number of runs
while isdir(fullfile(subjBehavDataDir, sprintf('run%d', nRunsAct + 1)))
    nRunsAct = nRunsAct + 1;
end

fprintf(1, 'INFO: Number of runs in expt.script = %d\n', nRunsScr);
fprintf(1, 'INFO: Number of runs actually delivered = %d\n', nRunsAct);

assert(nRunsAct <= nRunsScr);

%% Read the stim-schedule for the delivered runs
isInvalid = 0;
runinfo = {};
stims = {};
for i1 = 1 : nRunsAct
    runStr = sprintf('run%d', i1);
    
    stims{i1} = nan(1, expt.script.(runStr).nTrials);
%     spOnsets{i1} = nan(1, numel(runinfo{i1}.data.stim));

    trCnt = 1; % Trial count
    for i2 = 1 : expt.script.(runStr).nReps
        repStr = sprintf('rep%d', i2);
        
        for i3 = 1 : numel(expt.script.(runStr).(repStr).trialOrder)
            
            stims{i1}(trCnt) = min([3, expt.script.(runStr).(repStr).trialOrder(i3)]);
            trCnt = trCnt + 1;
        end
    end
    stims{i1} = [3, stims{i1}(1 : end - 1)];     % Padded one: the first volume doesn't capture any task-related neural responses
end

% 1 - Speech production: oral reading of disyllabic words
% 2 - Baseline
for i1 = 1 : numel(stims)
    [R1, R2] = gen_design_matrix(stims{i1}, TA, TR);
%     sess(1).R = R1(:,1:end/2);
    sess(i1).R = R1(:, 1 : end / 2);
    sess(i1).names = {};
    % sess(1).names = {'R_RHBP','R_VOICE','R_JAW','R_TONGUE','R_LIPS','R_NOISE'};%{'RHBP','Voice','Jaw','Listen'}; 
    sess(i1).onsets = {};
    sess(i1).durations = {};
    % end
    close all
end

% if iscell(sess)
%     sess=sess{1};
% end
save(fullfile(DATA_dir, subjID, 'funcloc_model.mat'), 'sess');
fprintf('%s saved.\n', fullfile(DATA_dir, subjID, 'funcloc_model.mat'));

%% Create the text files for the use of nipype
str = cell(1, numel(stims));
for i1 = 1 : numel(stims)
    stim = stims{i1};
    str{i1}.speech = sprintf('Run %d: speech: ', i1);
    str{i1}.baseline = sprintf('Run %d: baseline: ', i1);
    %if ~isempty(find(stim == 3))
    str{i1}.invalid = sprintf('Run %d: invalid: ', i1);
    %end
    
    for i2 = 1 : numel(stim)
        if stim(i2) == 1    % Valid speech trial 
            spOnset = tOK.voiceOnset{i1}(i2 - 1);
            if isnan(spOnset)
                error('NaN value found in spOnset');
            end
%             str{i1}.speech = [str{i1}.speech, ...
%                                 num2str((TA + (i2 - 1) * TR) / 1e3 + spOnset), ' '];
            str{i1}.speech = [str{i1}.speech, ...
                                num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
        elseif stim(i2) == 2    % Baseline trial
            str{i1}.baseline = [str{i1}.baseline, ...
                                num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
        elseif stim(i2) == 3    % Invalid speech trial
            str{i1}.invalid = [str{i1}.invalid, ...
                                num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
        end
    end
end

% Write to file
text_sched_fn = fullfile(DATA_dir, subjID, 'funcloc_sched.txt');
f1 = fopen(text_sched_fn, 'w');
for i1 = 1 : numel(str)
    fprintf(f1, '%s\n', str{i1}.speech);
    fprintf(f1, '%s\n', str{i1}.baseline);
    %if isfield(str{i1}, 'invalid')
        fprintf(f1, '%s\n', str{i1}.invalid);
    %end
end
fclose(f1);
fprintf('Stimulus schedlue written to file %s.\n', text_sched_fn);

return