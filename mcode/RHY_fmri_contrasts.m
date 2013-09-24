function RHY_fmri_contrasts(subjID)
% this script has to be executed inside the subject's directory.
if isequal(getHostName, 'ba3')
    dataDir = '/users/cais/RHY/DATA';
else
    dataDir = '/speechlab/5/scai/RHY/DATA';
end

check_dir(dataDir);

subjDataDir = fullfile(dataDir, strrep(subjID, 'MRI_', ''));
check_dir(subjDataDir);

%% Determine how many runs there are
desMatFN = fullfile(subjDataDir, 'fmri_model.mat');
check_file(desMatFN);

load(desMatFN);
assert(exist('sess', 'var') == 1);

nRuns = length(sess);
fprintf(1, 'INFO: Found %d runs in design matrix %s\n', ...
        nRuns, desMatFN);

%% example T contrasts cell array

% RegNameTContrasts{1} = {...
%     {'SvBL',    'Sn(1) R1', 1, 'Sn(1) R2', -1, ...
%                        'Sn(2) R1', 1, 'Sn(2) R2', -1, ...
%                        'Sn(3) R1', 1, 'Sn(3) R2', -1 ...
%                        }
%     };

%--- Speech vs Baseline ---%
RegNameTContrasts{1} = {'SvBL'};
for i1 = 1 : nRuns
   RegNameTContrasts{1}{end + 1} = sprintf('Sn(%d) R1', i1);
   RegNameTContrasts{1}{end + 1} = 0.5;
   
   RegNameTContrasts{1}{end + 1} = sprintf('Sn(%d) R2', i1);
   RegNameTContrasts{1}{end + 1} = 0.5;
   
   RegNameTContrasts{1}{end + 1} = sprintf('Sn(%d) R3', i1);
   RegNameTContrasts{1}{end + 1} = -1;
end

contrastFN = fullfile(subjDataDir, 'fmri_contrasts.mat');
save(contrastFN, 'RegNameTContrasts');
check_file(contrastFN);
fprintf('%s saved.\n', contrastFN);

return