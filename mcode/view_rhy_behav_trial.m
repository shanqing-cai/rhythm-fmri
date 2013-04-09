function view_rhy_behav_trial(dataFN, varargin)
%% 
% dataFN = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_9\run1\rep1\trial-2-1.mat';
% dataFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\run1\rep1\trial-6-2.mat';

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

%% Change parameters
% data.params.bPitchShift = 0;
% data.params.bBypassFmt = 0;

%%
figure('Position', [50, 100, 1400, 600]);
subplot('Position', [0.05, 0.5, 0.9, 0.4]);
show_spectrogram(data.signalIn, data.params.sr, 'noFig');
title(data.params.name);

frameDur = data.params.frameLen / data.params.sr;
tAxis = 0 : frameDur : frameDur * (size(data.fmts, 1) - 1);
plot(tAxis, data.fmts(:, 1 : 2), 'w-');



plot(tAxis, data.ost_stat * 250, 'b-');

% if isfield(handles, 'FmtShiftStat0') && isfield(handles, 'FmtShiftStat1')
FmtShiftStat0 = 5;
FmtShiftStat1 = 9;
    
ys = get(gca, 'YLim');
if ~isempty(find(data.ost_stat == FmtShiftStat0, 1))
    plot(repmat(tAxis(find(data.ost_stat == FmtShiftStat0, 1)), 1, 2), ys, 'b--');
end

if ~isempty(find(data.ost_stat == FmtShiftStat1 + 1, 1))
    plot(repmat(tAxis(find(data.ost_stat == FmtShiftStat1 + 1, 1)), 1, 2), ys, 'b-');
end
% end

ylabel('Frequency (Hz)');

%%
subplot('Position', [0.05, 0.1, 0.9, 0.4]);
show_spectrogram(data.signalOut, data.params.sr, 'noFig');

plot(tAxis, data.fmts(:, 1 : 2), 'w-');
plot(tAxis, data.sfmts(:, 1 : 2), 'c-');

xlabel('Time (s)');
ylabel('Frequency (Hz)');

%%
if ~isempty(fsic(varargin, '-p'))
    wavplay(data.signalIn, data.params.sr);
    wavplay(data.signalOut, data.params.sr);
end

return