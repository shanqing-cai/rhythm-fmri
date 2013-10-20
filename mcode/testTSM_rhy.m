function testTSM_rhy(varargin)
%% 
% dataFN = 'G:\DATA\RHYTHM-FMRI\TestExpt_behav_9\run1\rep1\trial-2-1.mat';
dataFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\run1\rep1\trial-6-1.mat';
% dataFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_M01\run2\rep2\trial-2-2 .mat';

ostFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\run1\R.ost';
pcfFN = 'G:\DATA\RHYTHM-FMRI\PILOT_ANS_F01\run1\twarp_R.pcf';

B_USE_OUT_FILT = 0;

%%
load(dataFN); % Gives data

%% Load filter coefficients
if B_USE_OUT_FILT
    addpath E:\speechres\signals\S14_filter;

    [h0, Fs] = load_filter('EQF_demoL.bin');
    h = resample(h0, data.params.sr, Fs);
end

%% Change parameters
data.params.bPitchShift = 1;
data.params.bBypassFmt = 0;

check_file(ostFN);
check_file(pcfFN);
TransShiftMex(8, ostFN);
TransShiftMex(9, pcfFN);

%%
fs = data.params.sr;

sigIn = data.signalIn;

sigIn = resample(sigIn, data.params.sr * data.params.downfact, fs);     
sigInCell = makecell(sigIn, data.params.frameLen * data.params.downfact);

TransShiftMex(6);   % Reset;\

% p.rmsClipThresh=0.01;
% p.bRMSClip=1;

%% Configure filter
if B_USE_OUT_FILT
    TransShiftMex(3, 'outfiltlen', length(h));
    TransShiftMex(3, 'outfircoeff', h);
    TransShiftMex(3, 'buseoutfilt', 1);
    fprintf(1, 'INFO: Configured and enabled output filter\n');
else
    TransShiftMex(3, 'buseoutfilt', 0);
end

%%
MexIO('init', data.params);

for n = 1 : length(sigInCell)
    TransShiftMex(5, sigInCell{n});
end

data1 = MexIO('getData');

% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
% [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));

%%
timeLims = [0.5, 2.0];

figure('Position', [50, 100, 1400, 600]);
subplot('Position', [0.05, 0.5, 0.9, 0.4]);
set(gca, 'FontSize', 20);
show_spectrogram(data1.signalIn, data.params.sr, 'noFig');
set(gca, 'XLim', timeLims);
xs = get(gca, 'XLim');

subplot('Position', [0.05, 0.1, 0.9, 0.4]);
set(gca, 'FontSize', 20);
fprintf(1, 'INFO: Compensating left shift is used.\n')
show_spectrogram(data1.signalOut(data.params.pvocHop * 3 + 1 : end), ...
                 data.params.sr, 'noFig');
set(gca, 'XLim', xs);

% figure;
% subplot(2, 1, 1);
% plot(data.signalIn);
% 
% subplot(2, 1, 2);
% plot(data.signalOut);

% soundsc(data1.signalIn, data.params.sr);
soundsc(data1.signalOut, data.params.sr);

return