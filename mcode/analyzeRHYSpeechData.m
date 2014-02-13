function varargout = analyzeRHYSpeechData(subjID, varargin)
%% 
rhyConds = {'N', 'R'};
rhyConds_long = {'Non-rhythmic', 'Rhythmic'};

pertTypes = {'noPert', 'F1Up', 'decel', ...
             'noPert_precNoPert', 'noPert_precF1Up', 'noPert_precDecel'};
pertTypes_main = pertTypes(1 : 3);

FRAME_DUR = 2e-3;   % Unit: s (WARNING: ad hoc)

FMT_INTERP_N = 50;

AVG_FALLOFF = 0.75;     % The trace counter ratio threshold for trace averaging

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

colors.noPert = [0, 0, 0];
colors.F1Up = [1, 0, 1];
% colors.decel = [1, 0.5, 0];
colors.decel = [0.75, 0.25, 0];

colors.noPert_precNoPert = [0, 0, 0];
colors.noPert_precF1Up = [1, 0, 1];
colors.noPert_precDecel = [1, 0.5, 0];

lineTypes.N = 'o--';
lineTypes.R = 'o-';

pltLW = 1;

corrSigSquareClr = [0, 0, 0];
corrSigSquareLW = 1;
corrSigSquareMarkerSize = 10;

hostName = lower(getHostName);
if isequal(hostName, 'cns-pc34')
    dacacheDir = 'D:/DATA/RHYTHM-FMRI/_dacache';
    rawDataDir = 'D:/DATA/RHYTHM-FMRI';
else
    dacacheDir = '../dacache';
    rawDataDir = 'G:/DATA/RHYTHM-FMRI';
end

pdataFN = fullfile(dacacheDir, [subjID, '.mat']);
check_file(pdataFN);

T_INT_P_THRESH = 0.05;
CHI_P_THRESH = 0.05;

prodRatingThresh = 0; % Most liberal: 0; More conservative: 1

timeIntervals = {'s', 't1', 1, 3; ...
                 't1', 'd', 3, 5; ...
                 'd', 'b1', 5, 7; ...
                 'b1', 'g', 7, 10; ...
                 'g', 'b2', 10, 13; ...
                 'b2', 't2', 13, 16; ...
                 't2', 'p1', 16, 18};

P_THRESH_UNC = 0.05;
P_THRESH_CORR = 0.05;
             
             
%% Visualization options
fontSize = 14;

%% Additional input options
bMan = ~isempty(fsic(varargin, '--man')); % Use manual time labels

if ~isempty(fsic(varargin, '--prodRatingThresh'))
    prodRatingThresh = varargin{fsic(varargin, '--prodRatingThresh') + 1};
end

bTotSentDur = ~isempty(fsic(varargin, '--tot-sent-dur'));

%%
load(pdataFN);  % gives pdata

%% Determine FRAME_DUR
if isfile(pdata.mainData.rawDataFNs{1})
    load(pdata.mainData.rawDataFNs{1});
    FRAME_DUR = data.params.frameLen / data.params.sr;
    
    clear data;
    
    info_log(sprintf('WARNING: Cannot find raw data files. Assuming FRAME_DUR = %f s.', FRAME_DUR));
else
    info_log(sprintf('WARNING: Cannot find raw data files. Assuming FRAME_DUR = %f s.', FRAME_DUR), '-warn');
end

%% Process noPert trials according to preceding pertType
exptFN = fullfile(pdata.subject.dataDir, pdata.subject.name, 'expt.mat');
if ~isfile(exptFN) % Change expt file name
    exptFN = fullfile(rawDataDir, pdata.subject.name, 'expt.mat');
    
    if ~isfile(exptFN)
        error('Cannot find expt.mat');
    end
end
% check_file(exptFN);

assert(exist('expt', 'var') == 0);
load(exptFN);
assert(exist('expt', 'var') == 1);

[pdata, idx_noPert_precNoPert, idx_noPert_precF1Up, idx_noPert_precDecel] = ...
    proc_pdata_preceding(pdata, expt);
clear('expt');

%% Check whether the pert screening is complete
for i0 = 2 : 3
    t_pertType = pertTypes{i0};
    
    idx_pert = find(pdata.mainData.pertType == i0 - 1);
    idx_noScreen = find(isnan(pdata.mainData.bPertOkay(idx_pert)));
    idx_noScreen = idx_pert(idx_noScreen);

    if length(idx_noScreen) > 0
        fprintf(2, 'ERROR: the perturbation screening step of preprocessing is incomplete.\n');

        for i1 = 1 : numel(idx_noScreen)
            fprintf(2, '\t%s: phase=%s; rep=%d; trialNum=%d\n', ...
                    t_pertType, ...
                    pdata.mainData.phases{idx_noScreen(i1)}, ...
                    pdata.mainData.blockNums(idx_noScreen(i1)), ...
                    pdata.mainData.trialNums(idx_noScreen(i1)));
        end
        
        return;
    end
end

%%
idxRuns = [];
rn = 1;
runName = sprintf('run%d', rn);
while ~isempty(fsic(pdata.mainData.phases, runName))
    idxRuns = [idxRuns, fsic(pdata.mainData.phases, runName)];
    
    rn = rn + 1;
    runName = sprintf('run%d', rn);
end
isRun = zeros(1, length(pdata.mainData.phases));
idxRun(idxRuns) = 1;

for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    if isequal(rc, 'N'), rci = 0;
    else                 rci = 1;
    end
    
    for i2 = 1 : 3
        pt = pertTypes{i2};
        assert(isempty(strfind(pt, '_prec')));
        
        if isequal(pt, 'noPert')
            pti = 0;
            po = ones(size(pdata.mainData.pertType));
        elseif isequal(pt, 'F1Up')
            pti = 1;
            po = pdata.mainData.bPertOkay == 1;
        elseif isequal(pt, 'decel')
            pti = 2;
            po = pdata.mainData.bPertOkay == 1;
        else
            error();
        end
       
        idx.(rc).(pt) = find(idxRun == 1 & ...
                            pdata.mainData.rating > prodRatingThresh & ...
                            pdata.mainData.bRhythm == rci & ...
                            pdata.mainData.pertType == pti & ...
                            po);
    end

    %--- For between-trial effects ---%
    for i2 = 4 : 6
        pt = pertTypes{i2};
        assert(length(strfind(pt, '_prec')) == 1);
        
        if isequal(pt, 'noPert_precNoPert')
            ppti = 0;
        elseif isequal(pt, 'noPert_precF1Up')
            ppti = 1;
        elseif isequal(pt, 'noPert_precDecel')
            ppti = 2;
        else
            error();
        end
        
        idx.(rc).(pt) = find(idxRun == 1 & ...
                                    pdata.mainData.rating > prodRatingThresh & ...
                                    pdata.mainData.bRhythm == rci & ...
                                    pdata.mainData.pertType == 0 & ...
                                    pdata.mainData.precPertType == ppti);
    end

                            
end

%% Get CV if IVI stats and other vowel acoustic measures
vidx = get_vowel_indices(strrep(pdata.subject.pertSent, '_', ' '));
vidx = vidx(2 : end) - 2;

cvIVI = struct;
meanIVI = struct;
mn_cvIVI = struct;
sd_cvIVI = struct;
mn_meanIVI = struct;
sd_meanIVI = struct;

if bTotSentDur
    totSentDur = struct; % Total sentence duration
    mn_totSentDur = struct;
    sd_totSentDur = struct;
end

for i1 = 1 : numel(rhyConds)    
    rc = rhyConds{i1};

    for i2 = 1 : numel(pertTypes_main)
        pt = pertTypes{i2};

        cvIVIs.(rc).(pt) = nan(size(idx.(rc).(pt)));
        meanIVI.(rc).(pt) = nan(size(idx.(rc).(pt)));
        totSentDur.(rc).(pt) = nan(size(idx.(rc).(pt)));

        for i3 = 1 : numel(idx.(rc).(pt))
            t_idx = idx.(rc).(pt)(i3);
            
            asrTBeg = pdata.mainData.asrTBeg(:, t_idx);
            ts_onset = asrTBeg(vidx);
            ts_offset = asrTBeg(vidx + 1);
            ts_mid = (ts_onset + ts_offset) / 2;
            ivis = diff(ts_mid);
            
            meanIVI.(rc).(pt)(i3) = mean(ivis);
            cvIVI.(rc).(pt)(i3) = std(ivis) / mean(ivis);
            
            if bTotSentDur
                %-- Total sentence duration --%
                trial_asrDir = strrep(pdata.mainData.rawDataFNs{t_idx}, '.mat', '_asr');
                check_dir(trial_asrDir);
                asrOutTxt = fullfile(trial_asrDir, 'julian_stdout.txt');
                check_file(asrOutTxt);
                wavFN = fullfile(trial_asrDir, 'speech.wav');
                check_file(wavFN);
            
                pa = parse_asr_out(asrOutTxt, wavFN);
    
                if pa.nphns == 24 && isequal(pa.phones{1}, 'sil') && isequal(pa.phones{2}, 'dh') ...
                   && isequal(pa.phones{end - 1}, 's') && isequal(pa.phones{end}, 'sil')
                    totSentDur.(rc).(pt)(i3) = pa.tend(end - 1) - pa.tbeg(2);
                end
            end
        end
        
        
        mn_meanIVI.(rc).(pt) = nanmean(meanIVI.(rc).(pt));
        sd_meanIVI.(rc).(pt) = nanstd(meanIVI.(rc).(pt));

        mn_cvIVI.(rc).(pt) = nanmean(cvIVI.(rc).(pt));
        sd_cvIVI.(rc).(pt) = nanstd(cvIVI.(rc).(pt));
        
        if bTotSentDur
            %-- Total sentence duration --%
            mn_totSentDur.(rc).(pt) = nanmean(totSentDur.(rc).(pt)(i3));
            sd_totSentDur.(rc).(pt) = nanstd(totSentDur.(rc).(pt)(i3));
        end
        
    end
    
end

%--- Visualization ---%
figure;
set(gca, 'FontSize', fontSize);
hold on;
for i1 = 1 : numel(rhyConds)   
    rc = rhyConds{i1};

    bar(i1, mn_cvIVI.(rc).noPert, ...
        'EdgeColor', 'k', 'FaceColor', colors.(rc));
    plot([i1, i1], mn_cvIVI.(rc).noPert + [-1, 1] * sd_cvIVI.(rc).noPert, 'k-');
end
xlabel('Rhythm condition');
ylabel('CV of inter-vowel intervals (mean\pm1 SD)');
set(gca, 'XTick', 1 : numel(rhyConds), 'XTickLabel', rhyConds_long);


%% Contingency table for dysfluencies
analyze_dysf_fraction(pdata, rhyConds, rhyConds_long, ...
                      pertTypes, CHI_P_THRESH, colors);

%% Comparing ASR results on input and output
figure('Position', [50, 150, 1500, 600]);
spCnt = 1;
for i1 = 1 : numel(rhyConds)    
    rc = rhyConds{i1};
    for i2 = 1 : numel(pertTypes)
        pt = pertTypes{i2};
        
        subplot(numel(rhyConds), numel(pertTypes), spCnt);
        spCnt = spCnt + 1;
        
        tbegs = pdata.mainData.asrTBeg(:, idx.(rc).(pt));
        tbegs_FB = pdata.mainData.asrTBeg_FB(:, idx.(rc).(pt));
        
        plot(1e3 * (tbegs_FB - tbegs));        
        
        set(gca, 'YLim', [-50, 150]);
        
        set(gca, 'XTick', 1 : size(tbegs_FB - tbegs, 1));
        set(gca, 'XTickLabel', pdata.mainData.asrPhns);
        title(sprintf('%s - %s', rc, pt));
        
        %--- The average amount of tShift in t1 (t in the word "steady") ---%
        if isequal(pt, 'decel')
            idx_t1 = fsic(pdata.mainData.asrPhns, 't1');
            assert(length(idx_t1) == 1);
            mn_decelPert_tShift_t1.(rc) = nanmean(tbegs_FB(idx_t1, :) - tbegs(idx_t1, :));
            sd_decelPert_tShift_t1.(rc) = nanstd(tbegs_FB(idx_t1, :) - tbegs(idx_t1, :));
        end
    end
    

end

%% Manual labels
% int_s_t1 = struct;
% int_s_d = struct;
% int_s_b1 = struct;
% int_s_g = struct;
% int_s_b2 = struct;
% for i1 = 1 : numel(rhyConds)
%     rc = rhyConds{i1};
%     
%     for i2 = 1 : numel(pertTypes)
%         pt = pertTypes{i2};
%         
%         int_s_t1.(rc).(pt) = pdata.mainData.t1OnsetTime(idx.(rc).(pt)) - ...
%                             pdata.mainData.sOnsetTime(idx.(rc).(pt));
%         int_s_d.(rc).(pt) = pdata.mainData.dOnsetTime(idx.(rc).(pt)) - ...
%                             pdata.mainData.sOnsetTime(idx.(rc).(pt));
%         int_s_b1.(rc).(pt) = pdata.mainData.b1OnsetTime(idx.(rc).(pt)) - ...
%                             pdata.mainData.sOnsetTime(idx.(rc).(pt));
%         int_s_g.(rc).(pt) = pdata.mainData.gOnsetTime(idx.(rc).(pt)) - ...
%                             pdata.mainData.sOnsetTime(idx.(rc).(pt));
%         int_s_b2.(rc).(pt) = pdata.mainData.b2OnsetTime(idx.(rc).(pt)) - ...
%                             pdata.mainData.sOnsetTime(idx.(rc).(pt));
%     end
% end

%% Manual time labels
for i1 = 1 : size(timeIntervals, 1)
    eval(sprintf('man_%s_%s = struct;', timeIntervals{i1, 1}, timeIntervals{i1, 2}));
end

for i0 = 1 : size(timeIntervals, 1)
    tp1 = timeIntervals{i0, 1};
    tp2 = timeIntervals{i0, 2};
    
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        for i2 = 1 : numel(pertTypes)
            pt = pertTypes{i2};

            eval(sprintf('man_%s_%s.(rc).(pt) = pdata.mainData.%sOnsetTime(idx.(rc).(pt)) - pdata.mainData.%sOnsetTime(idx.(rc).(pt));', ...
                         tp1, tp2, tp2, tp1));
            eval(sprintf('man_%s_%s.(rc).(pt) = man_%s_%s.(rc).(pt)(~isnan(man_%s_%s.(rc).(pt)));', ...
                         tp1, tp2, tp1, tp2, tp1, tp2));
        end
    end
end
%% ASR time labels
for i1 = 1 : size(timeIntervals, 1)
    eval(sprintf('asr_%s_%s = struct;', timeIntervals{i1, 1}, timeIntervals{i1, 2}));
end

for i0 = 1 : size(timeIntervals, 1)
    tp1 = timeIntervals{i0, 1};
    tp2 = timeIntervals{i0, 2};
    tidx1 = timeIntervals{i0, 3};
    tidx2 = timeIntervals{i0, 4};
    
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        for i2 = 1 : numel(pertTypes)
            pt = pertTypes{i2};

            eval(sprintf('asr_%s_%s.(rc).(pt) = pdata.mainData.asrTBeg(%d, idx.(rc).(pt)) - pdata.mainData.asrTBeg(%d, idx.(rc).(pt));', ...
                         tp1, tp2, tidx2, tidx1));
            eval(sprintf('asr_%s_%s.(rc).(pt) = asr_%s_%s.(rc).(pt)(~isnan(asr_%s_%s.(rc).(pt)));', ...
                         tp1, tp2, tp1, tp2, tp1, tp2));

            % --- Cumulative intervals --- %
    %         asr_s_d.(rc).(pt) = pdata.mainData.asrTBeg(5, idx.(rc).(pt)) - ...
    %                             pdata.mainData.asrTBeg(1, idx.(rc).(pt));
    %         asr_s_d.(rc).(pt) = asr_s_d.(rc).(pt)(~isnan(asr_s_d.(rc).(pt)));
    %         
    %         asr_s_b1.(rc).(pt) = pdata.mainData.asrTBeg(7, idx.(rc).(pt)) - ...
    %                             pdata.mainData.asrTBeg(1, idx.(rc).(pt));
    %         asr_s_b1.(rc).(pt) = asr_s_b1.(rc).(pt)(~isnan(asr_s_b1.(rc).(pt)));
        end
    end
    
end

%% 
time_ints = {'s_t1', 't1_d', 'd_b1', 'b1_g', 'g_b2', 'b2_t2', 't2_p1'};
content_phones = {'s', 't eh', 'd iy', 'b ae t', 'g ey v', 'b er th', 't uw'};
assert(length(time_ints) == length(content_phones));

tIntChgs_mn = struct;
tIntChgs_pse = struct;  % Pooled standard errors
tIntChgs_tp = struct;   % P-values from unpaired t-tests

for h0 = 1 : 2
    if h0 == 1
        vis_pertTypes = pertTypes([2, 1, 3]);
        name = 'Timing change';
        fs = fontSize;
    else
        vis_pertTypes = pertTypes([5, 1, 6]);
        name = 'Timing change: between-trial effects';
        fs = fontSize * 0.6;
    end
    
    figure('Position', [100, 100, 800, 600], 'Name', name);
    tBarW = 0.05;
    tBarSpace = 0.02;

    hsp = nan(1, numel(rhyConds));
    for i0 = 1 : numel(rhyConds)
        rc = rhyConds{i0};
        hsp(i0) = subplot(2, 1, i0);
        set(gca, 'FontSize', fs);
        hold on;
        
        if h0 == 1
            tIntChgs_mn.(rc) = struct;
            tIntChgs_pse.(rc) = struct;
            tIntChgs_tp.(rc) = struct;
        end

        if isequal(rc, 'N')
            title('Non-rhythmic');
        else
            title('Rhythmic');
        end

        cum_vals = cell(1, 3);
        cum_mean = zeros(length(time_ints) + 1, 3);

        for i1 = 1 : numel(vis_pertTypes)
            pt = vis_pertTypes{i1};
            
            if (h0 == 1) && ~isequal(pt, pertTypes_main{1})
                tIntChgs_mn.(rc).(pt) = nan(1, numel(time_ints));
                tIntChgs_pse.(rc).(pt) = nan(1, numel(time_ints));
                tIntChgs_tp.(rc).(pt) = nan(1, numel(time_ints));
            end

            cum_ts = [];
            for i2 = 1 : numel(time_ints)
                if ~bMan % - Use ASR time labels - %
                    eval(sprintf('tlens = 1e3 * asr_%s.(rc).(pt);', time_ints{i2}));
                else % - Use manual time labels - %
                    eval(sprintf('tlens = 1e3 * man_%s.(rc).(pt);', time_ints{i2}));
                end
                
                if (h0 == 1) && ~isequal(pt, pertTypes_main{1})
                    eval(sprintf('tIntChgs_mn.(rc).(pt)(i2) = 1e3 * (mean(asr_%s.(rc).(pt)) - mean(asr_%s.(rc).noPert));', ...
                                 time_ints{i2}, time_ints{i2}));
                    eval(sprintf('sd_p = 1e3 * std(asr_%s.(rc).(pt));', time_ints{i2}));
                    eval(sprintf('sd_np = 1e3 * std(asr_%s.(rc).noPert);', time_ints{i2}));
                    eval(sprintf('n_p = length(asr_%s.(rc).(pt));', time_ints{i2}));
                    eval(sprintf('n_np = length(asr_%s.(rc).noPert);', time_ints{i2}));
                    pooled_sd = sqrt(((n_p - 1) * sd_p * sd_p + (n_np - 1) * sd_np * sd_np) / ...
                                     (n_p + n_np -2));
                    tIntChgs_pse.(rc).(pt)(i2) = pooled_sd * sqrt(1 / n_p + 1 / n_np);
                    
                    eval(sprintf('[~, t_p] = ttest2(asr_%s.(rc).(pt), asr_%s.(rc).noPert);', ...
                                 time_ints{i2}, time_ints{i2}));
                    tIntChgs_tp.(rc).(pt)(i2) = t_p;
                end
                
                if isempty(cum_ts) 
                    cum_ts = tlens;
                else
                    cum_ts = cum_ts + tlens;
                end

                %-- Show the mean --%
                rectangle('Position', [cum_mean(i2, i1), (tBarW + tBarSpace) * (i1 - 1), mean(tlens), tBarW], ...
                          'EdgeColor', colors.(pt));
                      
                      
                
                
                cum_mean(i2 + 1, i1) = cum_mean(i2, i1) + mean(tlens);

                if isequal(pt, 'noPert')
                    text(cum_mean(i2, i1) + mean(tlens) * 0.2, ...
                         (tBarW + tBarSpace) * (i1 - 0.6), ...
                         strrep(content_phones{i2}, ' ', '+'), ...
                         'Color', 'k', 'FontSize', fontSize - 2);
                end

                if isempty(cum_vals{i1})
                    cum_vals{i1} = [cum_vals{i1}, tlens(:)];
                else
                    cum_vals{i1} = [cum_vals{i1}, cum_vals{i1}(:, end) + tlens(:)];
                end
                
                %-- Show the standard error --%
                t_mn = mean(cum_vals{i1}(:, end));
                t_se = std(cum_vals{i1}(:, end)) / sqrt(length(cum_vals{i1}(:, end)));
                plot(t_mn + [-1; 1] * t_se, ...
                     repmat((tBarW + tBarSpace) * (i1 - 1) + 0.5 * tBarW, 1, 2), ...
                     '-', 'Color', colors.(pt));
                plot(repmat(t_mn - 1 * t_se, 1, 2), ...
                     (tBarW + tBarSpace) * (i1 - 1) + [0.25, 0.75] * tBarW, ...
                     '-', 'Color', colors.(pt));
                plot(repmat(t_mn + 1 * t_se, 1, 2), ...
                     (tBarW + tBarSpace) * (i1 - 1) + [0.25, 0.75] * tBarW, ...
                     '-', 'Color', colors.(pt));
%                 plot(cum_mean(i2, i1) + mean(tlens) + [-1; 1] * std(cum_ts) / sqrt(length(cum_ts)), ...
%                      repmat((tBarW + tBarSpace) * (i1 - 1) + 0.5 * tBarW, 1, 2), ...
%                      '-', 'Color', colors.(pt));
                
            end

        end

        for i2 = 1 : numel(time_ints) + 1
            if i2 > 1
                [~, p] = ttest2(cum_vals{1}(:, i2 - 1), cum_vals{2}(:, i2 - 1));          
            else
                p = 1;
            end
            if p >= T_INT_P_THRESH
                lc = [0.6, 0.6, 0.6];
                lw = 0.5;
            else
                lc = [0.0, 0.0, 0.0];
                lw = 2.5;
            end

            plot([cum_mean(i2, 2), cum_mean(i2, 1)], (tBarW + tBarSpace) * 1 + [0, -tBarSpace], ...
                 'Color', lc, 'LineWidth', lw);        

            if i2 > 1
                [~, p] = ttest2(cum_vals{3}(:, i2 - 1), cum_vals{2}(:, i2 - 1));
            else
                p = 1;
            end
            if p >= T_INT_P_THRESH
                lc = [0.6, 0.6, 0.6];
                lw = 0.5;
            else
                lc = [0.0, 0.0, 0.0];
                lw = 2.5;
            end
            plot([cum_mean(i2, 2), cum_mean(i2, 3)], (tBarW + tBarSpace) * 2 + [-tBarSpace, 0], ...
                 'Color', lc, 'LineWidth', lw);
        end

        xlabel('Time (ms)');
        set(gca, ...
            'YTick', tBarW / 2 + [0, 1, 2] * (tBarW + tBarSpace), ...
            'YTickLabel', vis_pertTypes); 

        xlims{i0} = get(gca, 'XLim');   
    end

    for i1 = 1 : numel(hsp)
        set(gcf, 'CurrentAxes', hsp(i1));
        set(gca, 'XLim', [0, max(xlims{1}(2), xlims{2}(2))]);
    end
end

%% Visualize the duration changes
figure('Position', [50, 150, 800, 600], 'Name', 'Summary of time-interval changes');
set(gca, 'FontSize', fontSize);

hsp = nan(1, numel(pertTypes_main(2 : end)));
for i1 = 1 : numel(pertTypes_main(2 : end))
    hsp(i1) = subplot(2, 1, i1);
    set(gca, 'FontSize', fontSize);
    hold on;
    box on;
end

for i1 = 1 : length(rhyConds)
    rc = rhyConds{i1};

    for i2 = 1 : numel(pertTypes_main(2 : end))
        pts = pertTypes_main(2 : end);
        pt = pts{i2};
        
        set(gcf, 'CurrentAxes', hsp(i2));
        errorbar(1 : length(tIntChgs_mn.(rc).(pt)), ...
                 tIntChgs_mn.(rc).(pt), ...
                 tIntChgs_pse.(rc).(pt), ...
                 lineTypes.(rc), 'Color', colors.(rc), ...
                 'LineWidth', pltLW);
             
    end
    
end

%--- Determine the YLim ---%
yss = nan(numel(hsp), 2);
for i1 = 1 : numel(hsp)
    set(gcf, 'CurrentAxes', hsp(i1));
    yss(i1, :) = get(gca, 'YLim');
end
YLim(1) = min(yss(:, 1));
YLim(2) = max(yss(:, 2));

for i2 = 1 : numel(pertTypes_main(2 : end))
    pt = pertTypes_main{i2 + 1};
    set(gcf, 'CurrentAxes', hsp(i2));
    
    set(gca, 'XTick', 1 : 7 , ...
        'XTickLabel', {'s', 't+eh', 'd+iy', 'b+ae+t', 'g+ey+v', 'b+er+th', 't+uw'});
    xs = get(gca, 'XLim');
    plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    
    set(gca, 'YLim', YLim);
    xlabel('Segments and Syllables');
%     ylabel(sprintf('Time inteval change from %s (ms) (mean\\pm1 SEM)', baseType));
    ylabel(sprintf('Duration change (ms, mean\\pm1 SE)'), ...
           'FontSize', fontSize);
    
    if i2 == 1
        legend(rhyConds_long);
    end
   
    title(['Perturbation type: ', strrep(pt, '_', '\_')]);
end

%% Time-intevral changes: Mark significant changes
for i1 = 1 : length(rhyConds)
    rc = rhyConds{i1};
    
    for i2 = 1 : numel(pertTypes_main(2 : end))
        pts = pertTypes_main(2 : end);
        pt = pts{i2};
        
        set(gcf, 'CurrentAxes', hsp(i2));
        
        for i3 = 1 : length(tIntChgs_tp.(rc).(pt))
            if tIntChgs_tp.(rc).(pt)(i3) < P_THRESH_UNC
                plot(i3, tIntChgs_mn.(rc).(pt)(i3), 'o', ...
                     'MarkerFaceColor', colors.(rc), ...
                     'MarkerEdgeColor', colors.(rc));
            end
            
            
            if ~isempty(fsic(varargin, '--perm-tint')) %--- Show permutation test results ---%
                if corrps_tIntChg(i3, i1, i2) < P_THRESH_CORR
                    plot(i3, tIntChgs_mn.(rc).(pt)(i3), 's', ...
                         'MarkerEdgeColor',  corrSigSquareClr, 'MarkerFaceColor', 'none', ...
                         'LineWidth', corrSigSquareLW, ...
                         'MarkerSize', corrSigSquareMarkerSize);
                end
            else  %--- Bonferroni correction ---%
                if tIntChgs_tp.(rc).(pt)(i3) < P_THRESH_UNC / numel(tIntChgs_tp.(rc).(pt))
                    plot(i3, tIntChgs_mn.(rc).(pt)(i3), 's', ...
                         'MarkerEdgeColor',  corrSigSquareClr, 'MarkerFaceColor', 'none', ...
                         'LineWidth', corrSigSquareLW, ...
                         'MarkerSize', corrSigSquareMarkerSize);
                end
            end
        end
    end
end

%% Show legend for significance symbols
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
x = xs(1) + 0.725 * range(xs);
y = ys(1) + 0.750 * range(ys);
w = 0.266 * range(xs);
h = 0.210 * range(ys);

show_legend_sig(x, y, w, h, rhyConds, pltLW, colors, corrSigSquareClr, corrSigSquareLW, corrSigSquareMarkerSize, ...
                P_THRESH_UNC, P_THRESH_CORR, '--bonferroni');


%% Formant trajectory analysis
analyze_fmts_vwls = {'eh', 'iy', 'ae', 'ey', 'er', 'uw', 'ah'};

idxphns = nan(1, length(analyze_fmts_vwls));

F1s = struct;
tnormF1s = struct;  % Time normalize F1s: 
spect = struct;

for i0 = 1 : numel(analyze_fmts_vwls)
    vwl = analyze_fmts_vwls{i0};
    assert(length(strmatch(vwl, pdata.mainData.asrPhns, 'exact')) == 1);
    idxphns(i0) = strmatch(vwl, pdata.mainData.asrPhns, 'exact');
    assert(~isnan(idxphns(i0)));
      
    F1s.(vwl) = struct;
    tnormF1s = struct;
    spec.(vwl) = struct;

    figure('Name', sprintf('Vowel formants under: %s', vwl));
    for i1 = 1 : numel(rhyConds)    
        rc = rhyConds{i1};

        subplot(1, numel(rhyConds), i1);
        title(rc);
        
        for i2 = 1 : numel(pertTypes_main)
            pt = pertTypes{i2};

            F1s.(vwl).(rc).(pt) = {};
            tnormF1s.(vwl).(rc).(pt) = {};

            for i3 = 1 : numel(idx.(rc).(pt))
                t_idx = idx.(rc).(pt)(i3);

                if pdata.mainData.rating(t_idx) == 0
                    continue;
                end           
                if pdata.mainData.bASROkay(t_idx) == 0
                    continue
                end

                F1s.(vwl).(rc).(pt){end + 1} = pdata.mainData.vwlFmts{t_idx}.(vwl)(:, 1);
                
                %--- Determine the ASR-determined begin and end time of
                %the vowel ---%
                t0 = pdata.mainData.asrTBeg(idxphns(i0), t_idx);     % Vowel begin
                t1 = pdata.mainData.asrTBeg(idxphns(i0) + 1, t_idx); % Vowel end
                t2 = pdata.mainData.asrTBeg(idxphns(i0) + 2, t_idx);
                
                idxm = round((t1 - t0) / FRAME_DUR);
                tnf = interp1(1 : idxm, pdata.mainData.vwlFmts{t_idx}.(vwl)(1 : idxm, 1), ...
                              linspace(1, idxm, FMT_INTERP_N));
                L = length(pdata.mainData.vwlFmts{t_idx}.(vwl)(:, 1));
                tnf = [tnf, interp1(idxm + 1 : L, pdata.mainData.vwlFmts{t_idx}.(vwl)(idxm + 1 : end, 1), ...
                                    linspace(idxm + 1, L, FMT_INTERP_N))];
                                
                tnormF1s.(vwl).(rc).(pt){end + 1} = tnf;
            end
            
            %--- Not time-normalized ---%
            aF1s.(vwl).(rc).(pt) = avgTrace1(F1s.(vwl).(rc).(pt));
            
            nFallOff = find(aF1s.(vwl).(rc).(pt)(:, 3) > aF1s.(vwl).(rc).(pt)(1, 3) * AVG_FALLOFF, 1, 'last');
            mnF1s.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(1 : nFallOff, 1);
            seF1s.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(1 : nFallOff, 2);
                        
            %--- Time-normalized ---%
            aF1s_tnorm.(vwl).(rc).(pt) = avgTrace1(tnormF1s.(vwl).(rc).(pt));
            
            mnF1s_tnorm.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(1 : nFallOff, 1);
            seF1s_tnorm.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(1 : nFallOff, 2);
            
            hold on;
            
            t_axis = 1e3 * (0 : FRAME_DUR : FRAME_DUR * (nFallOff - 1));
            plot(t_axis, mnF1s.(vwl).(rc).(pt), 'Color', colors.(pt));
        end
        
        xlabel('Time (ms)');
        legend(pertTypes_main);


    end
end




%% Visualization
% meas = asr_s_t1;
% measName = 'asr_s_t1';
% meta_rc_pt_plot(meas, measName, colors);
% 
% meas = asr_s_d;
% measName = 'asr_s_d';
% meta_rc_pt_plot(meas, measName, colors);
% 
% meas = asr_s_b1;
% measName = 'asr_s_b1';
% meta_rc_pt_plot(meas, measName, colors);
% 
% meas = asr_s_g;
% measName = 'asr_s_g';
% meta_rc_pt_plot(meas, measName, colors);
% 
% meas = asr_s_b2;
% measName = 'asr_s_b2';
% meta_rc_pt_plot(meas, measName, colors);

% meas = asr_s_t1;
% measName = 'asr_s_t1';
% meta_rc_pt_plot(meas, measName, colors);

%% Output
varargout = {};

if nargout == 1
    res = struct;
    
    res.mn_cvIVI = mn_cvIVI;
    res.sd_cvIVI = sd_cvIVI;
    
    res.mn_decelPert_tShift_t1 = mn_decelPert_tShift_t1;
    res.sd_decelPert_tShift_t1 = sd_decelPert_tShift_t1;
    
    res.time_ints = struct('asr_s_t1', asr_s_t1, ...
                           'asr_t1_d', asr_t1_d, ...
                           'asr_d_b1', asr_d_b1, ...
                           'asr_b1_g', asr_b1_g, ...
                           'asr_g_b2', asr_g_b2, ...
                           'asr_b2_t2', asr_b2_t2, ...
                           'asr_t2_p1', asr_t2_p1);
    res.aF1s = aF1s;
    res.aF1s_tnorm = aF1s_tnorm;
    
    if bTotSentDur
        res.mn_totSentDur = mn_totSentDur;
        res.sd_totSentDur = sd_totSentDur;
    end
        
    varargout{1} = res;
end


return