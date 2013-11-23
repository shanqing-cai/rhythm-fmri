function view_rhy_behav_trial(dataFN, bReproc, varargin)
%% 
% dataFN = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_9\run1\rep1\trial-2-1.mat';
% dataFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\run1\rep1\trial-6-2.mat';

%% Config
fontSize = 14;
lw = 1.5;
purple = [223, 0, 255] / 255;

%% Input options
bClean = ~isempty(fsic(varargin, '--clean'));
bShiftLeft = ~isempty(fsic(varargin, '--shift-left'));
bShowFmts = ~isempty(fsic(varargin, '--show-fmts'));

%% Process optional wild card
if ~isempty(strfind(dataFN, '*'))
    d = dir(dataFN);
    if length(d) == 0
        error('Wild card did not match any data files: %s', dataFN);
    else
        if length(d) > 1
            fprintf(1, 'WARNING: wild card matches more than one data files. Using only the first one: %s\n', d(1).name);
        end
        fDir = fileparts(dataFN);
        dataFN = fullfile(fDir, d(1).name);
    end
end

%%
load(dataFN); % Gives data

%% Search for any existing expt.mat files
[repDir, baseMatFN] = fileparts(dataFN);
[phaseDir, repBase] = fileparts(repDir);
[exptDir, phaseBase] = fileparts(phaseDir);

d_expt = dir(fullfile(exptDir, 'expt.mat'));
if length(d_expt) == 1
    load(fullfile(exptDir, 'expt.mat')); % Gives expt
    
    % -- Figure out the pertType -- %
    phase = phaseBase;
    repNum = str2double(strrep(repBase, 'rep', ''));
    splitted = strsplit(baseMatFN, '-');
    trialNum = str2double(splitted{2});
    trialType = str2double(splitted{3});
    
    pertType = expt.script.(phase).(['rep', num2str(repNum)]).pertType(trialNum);
    fprintf(1, 'INFO: trialType = %d; pertType = %d\n', trialType, pertType);
end

%% Search for any existing ost and pcf files
if trialType == 1
    d_ost = dir(fullfile(phaseDir, 'N.ost'));
elseif trialType == 2
    d_ost = dir(fullfile(phaseDir, 'R.ost'));
else
    error('Unexpected trialType: %d', trialType);
end

if bReproc
    if length(d_ost) == 1
        ost_fn = fullfile(phaseDir, d_ost(1).name);
    
        fprintf(1, 'INFO: using native ost file: %s\n', ost_fn);
        TransShiftMex(8, ost_fn, 1);
    else
        fprintf(1, 'WARNING: no native ost file is found\n');
    end
end

%% Change parameters
% data.params.bPitchShift = 0;
% data.params.bBypassFmt = 0;

if bReproc
    MexIO('init', data.params);
    MexIO('reset');
    
    fs = data.params.sr;

    sigIn = data.signalIn;

    sigIn = resample(sigIn, data.params.sr * data.params.downfact, fs);     
    sigInCell = makecell(sigIn, data.params.frameLen * data.params.downfact);

    for n = 1 : length(sigInCell)
        TransShiftMex(5, sigInCell{n});
    end

    data = MexIO('getData');
end

%%
figure('Position', [50, 100, 1400, 600]);
subplot('Position', [0.05, 0.5, 0.9, 0.4]);
set(gca, 'FontSize', fontSize);

show_spectrogram(data.signalIn, data.params.sr, 'noFig');
title(data.params.name, 'FontSize', fontSize * 1.2);

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0 : frameDur : frameDur * (size(data.fmts, 1) - 1);

if ~bClean
    plot(tAxis, data.ost_stat * 250, 'b-');
end
if ~bClean || bShowFmts
    plot(tAxis, data.fmts(:, 1 : 2), 'w--', 'LineWidth', lw); 
end

% if isfield(handles, 'FmtShiftStat0') && isfield(handles, 'FmtShiftStat1')
FmtShiftStat0 = 5;
FmtShiftStat1 = 9;
    
ys = get(gca, 'YLim');

if ~bClean
    if ~isempty(find(data.ost_stat == FmtShiftStat0, 1))
        plot(repmat(tAxis(find(data.ost_stat == FmtShiftStat0, 1)), 1, 2), ys, 'b--');
    end

    if ~isempty(find(data.ost_stat == FmtShiftStat1 + 1, 1))
        plot(repmat(tAxis(find(data.ost_stat == FmtShiftStat1 + 1, 1)), 1, 2), ys, 'b-');
    end
end
% end

ylabel('Frequency (Hz)');

%%
subplot('Position', [0.05, 0.1, 0.9, 0.4]);
set(gca, 'FontSize', fontSize);

if bShiftLeft
    data.signalOut = data.signalOut(data.params.pvocFrameLen + 1 : end);
end

show_spectrogram(data.signalOut, data.params.sr, 'noFig');

if ~bClean || bShowFmts
    plot(tAxis, data.fmts(:, 1 : 2), 'w--', 'LineWidth', lw);
    plot(tAxis, data.sfmts(:, 1 : 2), 'k--', 'LineWidth', lw);
end

xlabel('Time (s)');
ylabel('Frequency (Hz)');

%%
if ~isempty(fsic(varargin, '-p'))
    wavplay(data.signalIn, data.params.sr);
    wavplay(data.signalOut, data.params.sr);
end

return