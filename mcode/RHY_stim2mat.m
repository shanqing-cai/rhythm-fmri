function RHY_stim2mat(subjID, TR, TA, varargin)
% this script has to be executed inside the subject's directory.
% read subjid from directory name
% Inputs:
%       TR - unit: ms
%       TA - unit: ms
%
%% Paths
if ~isempty(fsic(varargin, '--host'))
    hostName = varargin{fsic(varargin, '--host') + 1};
else
    hostName = getHostName;
end

if isequal(hostName, 'ba3')
    behavDataDir = '/users/cais/RHY/BEHAV_DATA';
    dacacheDir = '/users/cais/RHY/BEHAV_DATA/_dacache';
    dataDir = '/users/cais/RHY/DATA';
    genDesMatPath = '/users/cais/RHY/scripts';
else
    behavDataDir = '/speechlab/5/scai/RHY_BEHAV_DATA';
    dacacheDir = '/speechlab/5/scai/RHY_BEHAV_DATA/_dacache';
    dataDir = '/speechlab/5/scai/RHY/DATA';
    genDesMatPath = '/speechlab/5/scai/RHY/scripts';
end

%% Parameters
% TR = 11500;    % ms
% TA = 2470;    % ms

% DATA_dir = '/users/cais/STUT/DATA';
% funcloc_dir = '/users/cais/STUT/funcloc';

%%
if ~(length(subjID) > 4 && isequal(subjID(1 : 4), 'MRI_'))
    fprintf(2, 'WARNING: the subjID for an MRI session usually begins with MRI_. \nThe subjID you inputted (%s) does not fit this expectation.\n', subjID);
end

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

%% Load dacache file
pdataFN = fullfile(dacacheDir, sprintf('%s.mat', subjID));
check_file(pdataFN);

load(pdataFN);
assert(exist('pdata', 'var') == 1);

if ~isfield(pdata, 'mainData')
    error('Cannot find mainData in pdata');
end

if ~isfield(pdata.mainData, 'rating');
    error('Cannot find mainData in pdata.mainData');
end

if ~isfield(pdata.mainData, 'fluencyCode')
    error('Cannot find fluencyCode in pdata.mainData');
end

%% Read the stim-schedule for the delivered runs
dt = TR / 1e3 / 24 / 3600;

isInvalid = 0;
runinfo = {};
stims = {};
for i1 = 1 : nRunsAct
    runStr = sprintf('run%d', i1);
    
    stims{i1} = nan(1, expt.script.(runStr).nTrials);
%     spOnsets{i1} = nan(1, numel(runinfo{i1}.data.stim));

    ts = [];

    trCnt = 1; % Trial count
    for i2 = 1 : expt.script.(runStr).nReps
        repStr = sprintf('rep%d', i2);
        
        for i3 = 1 : numel(expt.script.(runStr).(repStr).trialOrder)
            %-- Load raw data file --%
            dataFN = dir(fullfile(subjBehavDataDir, runStr, repStr, ...
                                  sprintf('trial-%d-*.mat', i3)));
            if length(dataFN) ~= 1
                error('Cannot find exactly one .mat data file for trial: %s, %s, trial #%s', ...
                      runStr, repStr, i3);
            end
            dataFN = fullfile(subjBehavDataDir, runStr, repStr, dataFN(1).name);
            
            load(dataFN);
            assert(exist('data') == 1);
            
            ts(end + 1) = datenum(data.timeStamp);
            
            %-- --%
            stims{i1}(trCnt) = min([3, expt.script.(runStr).(repStr).trialOrder(i3)]);
            
            if stims{i1}(trCnt) <= 2
                %-- Locate trial in pdata --%
                idxPhase = fsic(pdata.mainData.phases, runStr);
                idxBlock = find(pdata.mainData.blockNums == i2);
                idxTrial = find(pdata.mainData.trialNums == i3);
                idx = intersect(intersect(idxPhase, idxBlock), idxTrial);

                if numel(idx) ~= 1
                    error('Cannot find exactly one entry for trial: [%s, %s, trial-%d] in pdata', ...
                          runStr, repStr, i3);
                end

                %-- Check for fluency and speech error --%
                if ~isempty(pdata.mainData.fluencyCode{idx})
                    stims{i1}(trCnt) = 4;
                elseif pdata.mainData.rating(idx) == 0
                    stims{i1}(trCnt) = 5;
                end

            end
            
            trCnt = trCnt + 1;
        end
    end
    stims{i1} = [3, stims{i1}(1 : end - 1)];     % Padded one: the first volume doesn't capture any task-related neural responses
    
    %-- Adjust stims according to time stamps: in case trigger skips exist
    %--%
    tsteps = diff(ts);
    len0 = length(stims{i1});
    nSkips = 0;
    for j1 = 2 : len0
        if tsteps(j1 - 1) > dt * 1.5
            stims{i1} = [stims{i1}(1 : j1 - 1), 3, ...
                         stims{i1}(j1 : end)];
            nSkips = nSkips + 1;
        end
    end
    
    assert(length(stims{i1}) == len0 + nSkips);
    if nSkips > 0
        fprintf(2, 'WARNING: found %d trigger skips in %s\n', ...
                nSkips, runStr);
        stims{i1} = stims{i1}(1 : len0);
    end    
end

%% Generate the design matrices
% 1 - Speech production: normal / non-rhytmic (SN)
% 2 - Speech production: rhythmic (SR)
% 3 - Baseline: no speech: (BL)
% 4 - Disfluent
% 5 - Error (not disfluent)

tPath = which('gen_design_matrix');
if isempty(tPath)
    addpath(genDesMatPath);
end

tPath = which('gen_design_matrix');
assert(~isempty(tPath));

for i1 = 1 : numel(stims)
    %--- The following lines are necessary because gen_design_matrix
    %    expects continuous trial type numbers with no skips ---%
    if ~isempty(find(stims{i1} == 5)) && isempty(find(stims{i1} == 4))
        stims{i1}(stims{i1} == 5) = 4;
    end
    
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

%% Save design mat file
check_dir(dataDir);

subjDataDir = fullfile(dataDir, strrep(subjID, 'MRI_', ''));
check_dir(subjDataDir, '-create');

modelFN = fullfile(subjDataDir, 'fmri_model.mat');
save(modelFN, 'sess');
check_file(modelFN);
fprintf(1, 'Design matrix saved to file: %s\n', modelFN);

%% Create the text files for the use of nipype
% str = cell(1, numel(stims));
% for i1 = 1 : numel(stims)
%     stim = stims{i1};
%     str{i1}.speech = sprintf('Run %d: speech: ', i1);
%     str{i1}.baseline = sprintf('Run %d: baseline: ', i1);
%     %if ~isempty(find(stim == 3))
%     str{i1}.invalid = sprintf('Run %d: invalid: ', i1);
%     %end
%     
%     for i2 = 1 : numel(stim)
%         if stim(i2) == 1    % Valid speech trial 
%             spOnset = tOK.voiceOnset{i1}(i2 - 1);
%             if isnan(spOnset)
%                 error('NaN value found in spOnset');
%             end
% %             str{i1}.speech = [str{i1}.speech, ...
% %                                 num2str((TA + (i2 - 1) * TR) / 1e3 + spOnset), ' '];
%             str{i1}.speech = [str{i1}.speech, ...
%                                 num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
%         elseif stim(i2) == 2    % Baseline trial
%             str{i1}.baseline = [str{i1}.baseline, ...
%                                 num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
%         elseif stim(i2) == 3    % Invalid speech trial
%             str{i1}.invalid = [str{i1}.invalid, ...
%                                 num2str((TA + (i2 - 1) * TR) / 1e3), ' '];
%         end
%     end
% end
% 
% % Write to file
% text_sched_fn = fullfile(DATA_dir, subjID, 'funcloc_sched.txt');
% f1 = fopen(text_sched_fn, 'w');
% for i1 = 1 : numel(str)
%     fprintf(f1, '%s\n', str{i1}.speech);
%     fprintf(f1, '%s\n', str{i1}.baseline);
%     %if isfield(str{i1}, 'invalid')
%         fprintf(f1, '%s\n', str{i1}.invalid);
%     %end
% end
% fclose(f1);
% fprintf('Stimulus schedlue written to file %s.\n', text_sched_fn);

return