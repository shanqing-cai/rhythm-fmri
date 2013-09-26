function [exptType] = guessExptType(expDir)
%% 
% Currently avaiable output types:
%     rand-twarp-fmt - randomized, time-warp & fmt (e.g., rhythm
%       behavioral experiment)
%     rand-RHY-fmri  - randomized, no perturbation, RHY paradigm (e.g., rhythm
%       fMRI experiment)
%     sust-fmt  - sustained perturbation, formant
%     rand-fmt  - randomized perturbation, formant (e.g., STUT_EH)
%%
if ~isdir(expDir)
    error_log(sprintf('Cannot find expDir: %s', expDir));
end

exptFN = fullfile(expDir, 'expt.mat');
if ~isfile(exptFN)
    error_log('Cannot find expt.mat file');
end

load(exptFN);   % gives expt

exptType = '';
stayDir = fullfile(expDir, 'stay');
run1Dir = fullfile(expDir, 'run1');

if isdir(stayDir);
    exptType = sprintf('%ssust-', exptType);
    
    %--- Find the first trial of the first rep in phase stay ---%
    rep1Dir = fullfile(stayDir, 'rep1');
    dfn = dir(fullfile(rep1Dir, 'trial-*-1.mat'));
    if length(dfn) == 0
        info_log('Experiment type guessing failed', '-warn');
        return
    end
    
    bFoundPertField = 0;
    for i1 = 1 : numel(dfn)
        load(fullfile(rep1Dir, dfn(i1).name));   % gives data
        
        if isfield(data.params, 'pertAmp') && ~isempty(find(data.params.pertAmp ~= 0))
            bFoundPertField = 1;
            break;
        end
        
        clear data;
    end
    
    if bFoundPertField == 1
        exptType = sprintf('%sfmt', exptType);
        return
    end
    
elseif isdir(run1Dir)
    exptType = sprintf('%srand-', exptType);
    
    fmtsPcf = dir(fullfile(run1Dir, 'fmt_*.pcf'));
    twarpPcf = dir(fullfile(run1Dir, 'twarp_*.pcf'));
    
    if ~isempty(fmtsPcf) && ~isempty(twarpPcf)
        exptType = sprintf('%stwarp-fmt', exptType);
        return
    elseif isempty(fmtsPcf) && isempty(twarpPcf)
        if isfield(expt.subject, 'trigByScanner') && expt.subject.trigByScanner == 1
            exptType = sprintf('%sRHY-fmri', exptType);
            return
        else
            info_log('Experiment type guessing failed', '-warn');
            return
        end
    else
        info_log('Experiment type guessing failed', '-warn');
        return
    end
else
    info_log('Experiment type guessing failed', '-warn');
    return
end

return