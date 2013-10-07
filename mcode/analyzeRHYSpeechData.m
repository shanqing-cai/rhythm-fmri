function varargout = analyzeRHYSpeechData2(subjID, varargin)
%% 
rhyConds = {'N', 'R'};
pertTypes = {'noPert', 'F1Up', 'decel', ...
             'noPert_precNoPert', 'noPert_precF1Up', 'noPert_precDecel'};

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

colors.noPert = [0, 0, 0];
colors.F1Up = [1, 0, 1];
colors.decel = [1, 0.5, 0];

colors.noPert_precNoPert = [0, 0, 0];
colors.noPert_precF1Up = [1, 0, 1];
colors.noPert_precDecel = [1, 0.5, 0];


dacacheDir = '../dacache';

pdataFN = fullfile(dacacheDir, [subjID, '.mat']);
check_file(pdataFN);

T_INT_P_THRESH = 0.05;

prodRatingThresh = 0; % Most liberal: 0; More conservative: 1

timeIntervals = {'s', 't1', 1, 3; ...
                 't1', 'd', 3, 5; ...
                 'd', 'b1', 5, 7; ...
                 'b1', 'g', 7, 10; ...
                 'g', 'b2', 10, 13; ...
                 'b2', 't2', 13, 16; ...
                 't2', 'p1', 16, 18};

%% Additional input options
bMan = ~isempty(fsic(varargin, '--man')); % Use manual time labels

if ~isempty(fsic(varargin, '--prodRatingThresh'))
    prodRatingThresh = varargin{fsic(varargin, '--prodRatingThresh') + 1};
end

%%
load(pdataFN);  % gives pdata

%% Process noPert trials according to preceding pertType
exptFN = fullfile(pdata.subject.dataDir, pdata.subject.name, 'expt.mat');
check_file(exptFN);

assert(exist('expt') == 0);
load(exptFN);
assert(exist('expt') == 1);

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

fontSize = 14;

for h0 = 1 : 2
    if h0 == 1
        vis_pertTypes = pertTypes([2, 1, 3]);
        name = 'Timing change';
        fs = fontSize;
    else
        vis_pertTypes = pertTypes([5, 1, 6]);
        name = 'Timing change (between-trial effects';
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

        if isequal(rc, 'N')
            title('Non-rhythmic');
        else
            title('Rhythmic');
        end

        cum_vals = cell(1, 3);
        cum_mean = zeros(length(time_ints) + 1, 3);

        for i1 = 1 : numel(vis_pertTypes)
            pt = vis_pertTypes{i1};

            for i2 = 1 : numel(time_ints)
                if ~bMan % - Use ASR time labels - %
                    eval(sprintf('tlens = 1e3 * asr_%s.(rc).(pt);', time_ints{i2}));
                else % - Use manual time labels - %
                    eval(sprintf('tlens = 1e3 * man_%s.(rc).(pt);', time_ints{i2}));
                end

                rectangle('Position', [cum_mean(i2, i1), (tBarW + tBarSpace) * (i1 - 1), mean(tlens), tBarW], ...
                          'EdgeColor', colors.(pt));
                cum_mean(i2 + 1, i1) = cum_mean(i2, i1) + mean(tlens);

                if isequal(pt, 'noPert')
                    text(cum_mean(i2, i1) + mean(tlens) * 0.2, ...
                         (tBarW + tBarSpace) * (i1 - 0.6), content_phones{i2}, ...
                         'Color', 'k', 'FontSize', fontSize - 2);
                end

                if isempty(cum_vals{i1})
                    cum_vals{i1} = [cum_vals{i1}, tlens(:)];
                else
                    cum_vals{i1} = [cum_vals{i1}, cum_vals{i1}(:, end) + tlens(:)];
                end
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

%% Formant trajectory analysis
analyze_fmts_vwls = {'eh', 'iy', 'ae'};

F1s = struct;
spect = struct;

for i0 = 1 : numel(analyze_fmts_vwls)
    vwl = analyze_fmts_vwls{i0};
      
    F1s.(vwl) = struct;
    spec.(vwl) = struct;

    figure('Name', sprintf('Vowel formants under: %s', vwl));
    for i1 = 1 : numel(rhyConds)    
        rc = rhyConds{i1};

        subplot(1, numel(rhyConds), i1);
        title(rc);
        
        for i2 = 1 : numel(pertTypes)
            pt = pertTypes{i2};

            F1s.(vwl).(rc).(pt) = {};

            for i3 = 1 : numel(idx.(rc).(pt))
                t_idx = idx.(rc).(pt)(i3);

                if pdata.mainData.rating(t_idx) == 0
                    continue;
                end           
                if pdata.mainData.bASROkay(t_idx) == 0
                    continue
                end

                F1s.(vwl).(rc).(pt){end + 1} = pdata.mainData.vwlFmts{t_idx}.(vwl)(:, 1);
            end

            aF1s.(vwl).(rc).(pt) = avgTrace1(F1s.(vwl).(rc).(pt));
            mnF1s.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(:, 1);
            seF1s.(vwl).(rc).(pt) = aF1s.(vwl).(rc).(pt)(:, 2);

            hold on;
            plot(mnF1s.(vwl).(rc).(pt), 'Color', colors.(pt));
        end


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
    res.time_ints = struct('asr_s_t1', asr_s_t1, ...
                           'asr_t1_d', asr_t1_d, ...
                           'asr_d_b1', asr_d_b1, ...
                           'asr_b1_g', asr_b1_g, ...
                           'asr_g_b2', asr_g_b2, ...
                           'asr_b2_t2', asr_b2_t2, ...
                           'asr_t2_p1', asr_t2_p1);
    res.aF1s = aF1s;
    varargout{1} = res;
end


return