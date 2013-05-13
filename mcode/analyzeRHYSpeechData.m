function varargout = analyzeRHYSpeechData(subjID, varargin)
%% 
rhyConds = {'N', 'R'};
pertTypes = {'noPert', 'F1Up', 'decel'};

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

colors.noPert = [0, 0, 0];
colors.F1Up = [1, 0, 1];
colors.decel = [1, 0.5, 0];

dacacheDir = '../dacache';

pdataFN = fullfile(dacacheDir, [subjID, '.mat']);
check_file(pdataFN);

%%
load(pdataFN);  % gives pdata

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

idx.N.noPert = find(idxRun == 1 & ...
                    pdata.mainData.rating > 0 & ...
                    pdata.mainData.bRhythm == 0 & ...
                    pdata.mainData.pertType == 0);
idx.N.F1Up = find(idxRun == 1 & ...
                  pdata.mainData.rating > 0 & ...
                  pdata.mainData.bRhythm == 0 & ...
                  pdata.mainData.pertType == 1);
idx.N.decel = find(idxRun == 1 & ...
                   pdata.mainData.rating > 0 & ...
                   pdata.mainData.bRhythm == 0 & ...
                   pdata.mainData.pertType == 2);
                
idx.R.noPert = find(idxRun == 1 & ...
                    pdata.mainData.rating > 0 & ...
                    pdata.mainData.bRhythm == 1 & ...
                    pdata.mainData.pertType == 0);
idx.R.F1Up = find(idxRun == 1 & ...
                  pdata.mainData.rating > 0 & ...
                  pdata.mainData.bRhythm == 1 & ...
                  pdata.mainData.pertType == 1);
idx.R.decel = find(idxRun == 1 & ...
                   pdata.mainData.rating > 0 & ...
                   pdata.mainData.bRhythm == 1 & ...
                   pdata.mainData.pertType == 2);
%% Comparing ASR results on input and output
figure;
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
    end
end

%% Manual labels
int_s_t1 = struct;
int_s_d = struct;
int_s_b1 = struct;
int_s_g = struct;
int_s_b2 = struct;
for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    
    for i2 = 1 : numel(pertTypes)
        pt = pertTypes{i2};
        
        int_s_t1.(rc).(pt) = pdata.mainData.t1OnsetTime(idx.(rc).(pt)) - ...
                            pdata.mainData.sOnsetTime(idx.(rc).(pt));
        int_s_d.(rc).(pt) = pdata.mainData.dOnsetTime(idx.(rc).(pt)) - ...
                            pdata.mainData.sOnsetTime(idx.(rc).(pt));
        int_s_b1.(rc).(pt) = pdata.mainData.b1OnsetTime(idx.(rc).(pt)) - ...
                            pdata.mainData.sOnsetTime(idx.(rc).(pt));
        int_s_g.(rc).(pt) = pdata.mainData.gOnsetTime(idx.(rc).(pt)) - ...
                            pdata.mainData.sOnsetTime(idx.(rc).(pt));
        int_s_b2.(rc).(pt) = pdata.mainData.b2OnsetTime(idx.(rc).(pt)) - ...
                            pdata.mainData.sOnsetTime(idx.(rc).(pt));
    end
end

%% ASR time labels
asr_s_t1 = struct;
asr_t1_d = struct;
asr_d_b1 = struct;
asr_b1_g = struct;
asr_g_b2 = struct;
asr_b2_t2 = struct;
asr_t2_p1 = struct;

for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    
    for i2 = 1 : numel(pertTypes)
        pt = pertTypes{i2};
        
        asr_s_t1.(rc).(pt) = pdata.mainData.asrTBeg(3, idx.(rc).(pt)) - ...
                            pdata.mainData.asrTBeg(1, idx.(rc).(pt));
        asr_s_t1.(rc).(pt) = asr_s_t1.(rc).(pt)(~isnan(asr_s_t1.(rc).(pt)));
        
        asr_t1_d.(rc).(pt) = pdata.mainData.asrTBeg(5, idx.(rc).(pt)) - ...
                            pdata.mainData.asrTBeg(3, idx.(rc).(pt));
        asr_t1_d.(rc).(pt) = asr_t1_d.(rc).(pt)(~isnan(asr_t1_d.(rc).(pt)));
        
        asr_d_b1.(rc).(pt) = pdata.mainData.asrTBeg(7, idx.(rc).(pt)) - ...
                             pdata.mainData.asrTBeg(5, idx.(rc).(pt));
        asr_d_b1.(rc).(pt) = asr_d_b1.(rc).(pt)(~isnan(asr_d_b1.(rc).(pt)));
        
        asr_b1_g.(rc).(pt) = pdata.mainData.asrTBeg(10, idx.(rc).(pt)) - ...
                             pdata.mainData.asrTBeg(7, idx.(rc).(pt));
        asr_b1_g.(rc).(pt) = asr_b1_g.(rc).(pt)(~isnan(asr_b1_g.(rc).(pt)));
        
        asr_g_b2.(rc).(pt) = pdata.mainData.asrTBeg(13, idx.(rc).(pt)) - ...
                             pdata.mainData.asrTBeg(10, idx.(rc).(pt));
        asr_g_b2.(rc).(pt) = asr_g_b2.(rc).(pt)(~isnan(asr_g_b2.(rc).(pt)));
        
        asr_b2_t2.(rc).(pt) = pdata.mainData.asrTBeg(16, idx.(rc).(pt)) - ...
                              pdata.mainData.asrTBeg(13, idx.(rc).(pt));
        asr_b2_t2.(rc).(pt) = asr_b2_t2.(rc).(pt)(~isnan(asr_b2_t2.(rc).(pt)));
        
        asr_t2_p1.(rc).(pt) = pdata.mainData.asrTBeg(18, idx.(rc).(pt)) - ...
                              pdata.mainData.asrTBeg(16, idx.(rc).(pt));
        asr_t2_p1.(rc).(pt) = asr_t2_p1.(rc).(pt)(~isnan(asr_t2_p1.(rc).(pt)));
    end
end

%% 
asr_ints = {'s_t1', 't1_d', 'd_b1', 'b1_g', 'g_b2', 'b2_t2', 't2_p1'};
content_phones = {'s', 't eh', 'd iy', 'b ae t', 'g ey v', 'b er th', 't uw'};
assert(length(asr_ints) == length(content_phones));

figure('Position', [100, 100, 800, 600]);
tBarW = 0.05;
tBarSpace = 0.02;
vis_pertTypes = pertTypes([2, 1, 3]);

hsp = nan(1, numel(rhyConds));
for i0 = 1 : numel(rhyConds)
    rc = rhyConds{i0};
    hsp(i0) = subplot(2, 1, i0);
    hold on;
    
    if isequal(rc, 'N')
        title('Non-rhythmic');
    else
        title('Rhythmic');
    end
    
    cum_vals = cell(1, 3);
    cum_mean = zeros(length(asr_ints) + 1, 3);
    
    for i1 = 1 : numel(vis_pertTypes)
        pt = vis_pertTypes{i1};

        for i2 = 1 : numel(asr_ints)
            eval(sprintf('tlens = 1e3 * asr_%s.(rc).(pt);', asr_ints{i2}));

            rectangle('Position', [cum_mean(i2, i1), (tBarW + tBarSpace) * (i1 - 1), mean(tlens), tBarW], ...
                      'EdgeColor', colors.(pt));
            cum_mean(i2 + 1, i1) = cum_mean(i2, i1) + mean(tlens);
            
            if isequal(pt, 'noPert')
                text(cum_mean(i2, i1) + mean(tlens) * 0.2, ...
                     (tBarW + tBarSpace) * (i1 - 0.6), content_phones{i2}, ...
                     'Color', 'k');
            end
            
            if isempty(cum_vals{i1})
                cum_vals{i1} = [cum_vals{i1}, tlens(:)];
            else
                cum_vals{i1} = [cum_vals{i1}, cum_vals{i1}(:, end) + tlens(:)];
            end
        end
        
    end
    
    for i2 = 1 : numel(asr_ints) + 1
        if i2 > 1
            [~, p] = ttest2(cum_vals{1}(:, i2 - 1), cum_vals{2}(:, i2 - 1));          
        else
            p = 1;
        end
        if p >= 0.05
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
        if p >= 0.05
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


%% Formant trajectory analysis
analyze_fmts_vwls = {'eh', 'iy', 'ae'};

for i0 = 1 : numel(analyze_fmts_vwls)
    vwl = analyze_fmts_vwls{i0};
    
    F1s.(vwl) = struct;

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
meas = int_s_t1;
measName = 'int_s_t1';
meta_rc_pt_plot(meas, measName, colors);

meas = int_s_d;
measName = 'int_s_d';
meta_rc_pt_plot(meas, measName, colors);

meas = int_s_b1;
measName = 'int_s_b1';
meta_rc_pt_plot(meas, measName, colors);

meas = int_s_g;
measName = 'int_s_g';
meta_rc_pt_plot(meas, measName, colors);

meas = int_s_b2;
measName = 'int_s_b2';
meta_rc_pt_plot(meas, measName, colors);

meas = asr_s_t1;
measName = 'asr_s_t1';
meta_rc_pt_plot(meas, measName, colors);

%% Output
varargout = {};

if nargout == 1
    res = struct;
    res.asr_ints = struct('asr_s_t1', asr_s_t1, ...
                          'asr_t1_d', asr_t1_d, ...
                          'asr_d_b1', asr_d_b1, ...
                          'asr_b1_g', asr_b1_g, ...
                          'asr_g_b2', asr_g_b2, ...
                          'asr_b2_t2', asr_b2_t2, ...
                          'asr_t2_p1', asr_t2_p1);
    varargout{1} = res;
end


return