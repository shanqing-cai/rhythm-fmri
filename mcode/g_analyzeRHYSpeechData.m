function g_analyzeRHYSpeechData(subjListFN, varargin)
%%
% Optional input arguments:
%    -r | --reload:     reload data from individual subjects
%    --between-trial:   perform between-trial (adaptation) analysis  
%    --tnorm:           use time-normalized formants
%
%% Config
rhyConds = {'N', 'R'};
pertTypes_within = {'noPert', 'F1Up', 'decel'};
pertTypes_between = {'noPert_precNoPert', 'noPert_precF1Up', 'noPert_precDecel'};

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

colors.noPert = [0, 0, 0];
colors.F1Up = [1, 0, 1];
colors.decel = [1, 0.5, 0];

colors.noPert_precNoPert = [0, 0, 0];
colors.noPert_precF1Up = [1, 0, 1];
colors.noPert_precDecel = [1, 0.5, 0];

gray = [0.5, 0.5, 0.5];

P_THRESH_UNC = 0.05;

FMT_ANA_VWLS = {'eh', 'iy', 'ae', 'ey'};
MAX_FMT_LEN = 256;

T_STEP = 0.002;     % Unit: s

%% Read subject list
check_file(subjListFN);
subjIDs = read_struct_from_text(subjListFN);
subjIDs = splitstring(subjIDs.SUBJ_IDS, ',');
ns = numel(subjIDs);

%% Analysis type: within-trial or between-trial
if ~isempty(fsic(varargin, '--between-trial'))
    baseType = 'noPert_precNoPert';
    pertTypes = pertTypes_between;
else
    baseType = 'noPert';
    pertTypes = pertTypes_within;
end

bFmtTNorm = ~isempty(fsic(varargin, '--tnorm'));

%%
npt = numel(pertTypes);

tInts = struct;
% ==== Edit me ==== %
tIntItems = {'s_t1', 't1_d', 'd_b1', 'b1_g', 'g_b2', 'b2_t2', 't2_p1'};
% ==== ~Edit me ==== %

for i0 = 1 : numel(tIntItems)
    ti = tIntItems{i0};
    
    tInts.(ti) = struct;
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};
        tInts.(ti).(rc) = nan(ns, npt);
    end
end

%% Load data
grpCacheFN = sprintf('%s_%s_grpCache.mat', mfilename, ...
                           strrep(subjListFN, '.txt', ''));

bReload = ~isempty(fsic(varargin, '-r')) || ~isempty(fsic(varargin, '--reload'));
if bReload
    sres = cell(1, ns);
    for i1 = 1 : numel(subjIDs)
        sID = subjIDs{i1};
        sres{i1} = analyzeRHYSpeechData(sID);
        
        close all hidden;
        drawnow;
    end
    
    save(grpCacheFN, 'sres');
    check_file(grpCacheFN);
    fprintf(1, 'INFO: Saved sres to file: %s\n', grpCacheFN);
else
    check_file(grpCacheFN);
    load(grpCacheFN); % gives sres
    fprintf(1, 'INFO: Loaded sres to file: %s\n', grpCacheFN);
end

%% Time intervals
assert(isequal(pertTypes{1}, baseType));

for i0 = 1 : numel(tIntItems)
    ti = tIntItems{i0};
    
    for i1 = 1 : numel(subjIDs)
        sID = subjIDs{i1};

        for j1 = 1 : numel(rhyConds)
            rc = rhyConds{j1};

            for j2 = 1 : numel(pertTypes)
                pt = pertTypes{j2};

                tInts.(ti).(rc)(i1, j2) = nanmean(sres{i1}.time_ints.(['asr_', ti]).(rc).(pt));
            end
        end
    end
end

%% Formant frequenciesavg
avgVwlF1 = struct;
steVwlF1 = struct;

if bFmtTNorm
    aF1_fld = 'aF1s_tnorm';
else
    aF1_fld = 'aF1s';
end

for h1 = 1 : length(FMT_ANA_VWLS)
    t_vwl = FMT_ANA_VWLS{h1};
    
    avgVwlF1.(t_vwl) = struct;
    chgAvgVwlF1s.(t_vwl) = struct; % Average across Ss, then calculate the pert-noPert difference
        
    vwlF1s.(t_vwl) = struct;
    chgVwlF1s.(t_vwl) = struct;
    
    avgChgVwlF1s.(t_vwl) = struct;  % Calculate the pert-noPert difference, then average across Ss (probably better)
    steChgVwlF1s.(t_vwl) = struct;
    ptChgVwlF1s.(t_vwl) = struct; % P-values from t-test of the F1 changes
    
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        avgVwlF1s.(t_vwl).(rc) = struct;
        steVwlF1s.(t_vwl).(rc) = struct;        

        vwlF1s.(t_vwl).(rc) = struct;
        chgAvgVwlF1s.(t_vwl).(rc) = struct;
        
        avgChgVwlF1s.(t_vwl).(rc) = struct;
        steChgVwlF1s.(t_vwl).(rc) = struct;
        ptChgVwlF1s.(t_vwl).(rc) = struct;
        
        for i2 = 1 : numel(pertTypes)
            pt = pertTypes{i2};

            avgVwlF1s.(t_vwl).(rc).(pt) = nan(MAX_FMT_LEN, length(sres));
            steVwlF1s.(t_vwl).(rc).(pt) = nan(MAX_FMT_LEN, length(sres));

            vwlF1s.(t_vwl).(rc).(pt) = nan(MAX_FMT_LEN, length(sres));
            if i2 ~= 1
                chgVwlF1s.(t_vwl).(rc).(pt) = nan(MAX_FMT_LEN, length(sres));
            end

            %--- Collect formant data from subjects ---%
            for j1 = 1 : length(sres);
                len = size(sres{j1}.(aF1_fld).(t_vwl).(rc).(pt), 1);
                vwlF1s.(t_vwl).(rc).(pt)(1 : len, j1) = sres{j1}.(aF1_fld).(t_vwl).(rc).(pt)(:, 1);
                
                if i2 ~= 1
                    len = min([length(sres{j1}.(aF1_fld).(t_vwl).(rc).(pt)(:, 1)), length(sres{j1}.(aF1_fld).(t_vwl).(rc).noPert(:, 1))]);
                    chgVwlF1s.(t_vwl).(rc).(pt)(1 : len, j1) = sres{j1}.(aF1_fld).(t_vwl).(rc).(pt)(1 : len, 1) - sres{j1}.(aF1_fld).(t_vwl).(rc).noPert(1 : len, 1);                   
                end
            end

            avgVwlF1s.(t_vwl).(rc).(pt) = mean(vwlF1s.(t_vwl).(rc).(pt), 2);
            steVwlF1s.(t_vwl).(rc).(pt) = std(vwlF1s.(t_vwl).(rc).(pt), [], 2) / sqrt(size(vwlF1s.(t_vwl).(rc).(pt), 2));

            glen = find(~isnan(avgVwlF1s.(t_vwl).(rc).(pt)), 1, 'last');    % Group-level length
            avgVwlF1s.(t_vwl).(rc).(pt) = avgVwlF1s.(t_vwl).(rc).(pt)(1 : glen);
            steVwlF1s.(t_vwl).(rc).(pt) = steVwlF1s.(t_vwl).(rc).(pt)(1 : glen);
            
            if i2 ~= 1
                %-- Calculate the pert-noPert difference, then average across Ss--%
                avgChgVwlF1s.(t_vwl).(rc).(pt) = mean(chgVwlF1s.(t_vwl).(rc).(pt), 2);
                steChgVwlF1s.(t_vwl).(rc).(pt) = std(chgVwlF1s.(t_vwl).(rc).(pt), [], 2) / ...
                                                 sqrt(size(chgVwlF1s.(t_vwl).(rc).(pt), 2));
                                   
                glen = find(~isnan(avgChgVwlF1s.(t_vwl).(rc).(pt)), 1, 'last');
                avgChgVwlF1s.(t_vwl).(rc).(pt) = avgChgVwlF1s.(t_vwl).(rc).(pt)(1 : glen);
                steChgVwlF1s.(t_vwl).(rc).(pt) = steChgVwlF1s.(t_vwl).(rc).(pt)(1 : glen);
                
                %-- Calculation of p-values (t-test) --%
                ptChgVwlF1s.(t_vwl).(rc).(pt) = nan(glen, 1);
                for i3 = 1 : glen
                    [~, t_p] = ttest(chgVwlF1s.(t_vwl).(rc).(pt)(i3, :));
                    ptChgVwlF1s.(t_vwl).(rc).(pt)(i3) = t_p;
                end
            end
        end
        
        assert(isequal(pertTypes{1}, baseType));
        
        for i2 = 2 : numel(pertTypes)
            pt = pertTypes{i2};
            clen = min([size(avgVwlF1s.(t_vwl).(rc).(pertTypes{1}), 1), ...
                        size(avgVwlF1s.(t_vwl).(rc).(pt), 1)]);
            chgAvgVwlF1s.(t_vwl).(rc).(pt) = avgVwlF1s.(t_vwl).(rc).(pt)(1 : clen) - ...
                                             avgVwlF1s.(t_vwl).(rc).(pertTypes{1})(1 : clen);
        end
        
        
    end

end

%% 
tIntChgs = struct;

for i0 = 1 : numel(tIntItems)
    ti = tIntItems{i0};
    tIntChgs.(ti) = struct;
    
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        tIntChgs.(ti).(rc) = tInts.(ti).(rc)(:, 2 : 3) - repmat(tInts.(ti).(rc)(:, 1), 1, 2);
    end
end

%% Visualization: time intervals changes
figure('Name', sprintf('Changes in time intervals from %s: condition', baseType), ...
       'Position', [100, 100, 1200, 600]);

for i0 = 1 : numel(tIntItems)
    ti = tIntItems{i0};
    
%     figure('Name', sprintf('Changes in %s', ti), ...
%            'Position', [100, 100, 400, 300]);
    subplot(2, 4, i0);
    hold on;
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        errorbar(1 : npt - 1, 1e3 * mean(tIntChgs.(ti).(rc)), 1e3 * ste(tIntChgs.(ti).(rc)), ...
                 'o-', 'Color', colors.(rc));
    end
    legend(rhyConds, 'Location', 'Northwest');
    set(gca, 'XLim', [0, npt]);
    set(gca, 'XTick', [1, 2], 'XTickLabel', pertTypes(2 : end));
    xlim = get(gca, 'XLim');
    plot(xlim, [0, 0], '-', 'Color', gray);
    ylabel(sprintf('Change in time interval %s from %s (ms)', ...
                   baseType, strrep(ti, '_', '\_')));
	title(strrep(sprintf('Changes in %s', ti), '_', '\_'));
end

%% Visualization: time intervals changes 2
figure('Position', [50, 150, 800, 400], 'Name', 'Summary of time-interval changes');
hsp = nan(1, numel(pertTypes(2 : end)));
for i1 = 1 : numel(pertTypes(2 : end))
    hsp(i1) = subplot(1, 2, i1);
    hold on;
end

p_vals_tIntChg = struct;

for i1 = 1 : length(rhyConds)
    rc = rhyConds{i1};
    
    mn_TIntChgs.(rc) = nan(0, 2);
    se_TIntChgs.(rc) = nan(0, 2);
    p_vals_tIntChg.(rc) = nan(0, 2); % Columns: perturbation types; Col 1 - F1Up, Col 2 - Decel
    for j1 = 1 : numel(tIntItems)
        t_item = tIntItems{j1};
        
        mn_TIntChgs.(rc) = [mn_TIntChgs.(rc); mean(tIntChgs.(t_item).(rc))];
        se_TIntChgs.(rc) = [se_TIntChgs.(rc); ste(tIntChgs.(t_item).(rc))];
        
        t_p_vals = nan(1, size(tIntChgs.(t_item).(rc), 2));
        
        for j2 = 1 : size(tIntChgs.(t_item).(rc), 2)
            [~, t_p_vals(j2)] = ttest(tIntChgs.(t_item).(rc)(:, j2));
        end
        
        p_vals_tIntChg.(rc) = [p_vals_tIntChg.(rc); t_p_vals];
    end

    for i2 = 1 : numel(pertTypes(2 : end))
        set(gcf, 'CurrentAxes', hsp(i2));
        errorbar(1 : size(mn_TIntChgs.(rc), 1), ...
                 1e3 * mn_TIntChgs.(rc)(:, i2), ...
                 1e3 * se_TIntChgs.(rc)(:, i2), ...
                 'o-', 'Color', colors.(rc));
             
    end
    
end

for i2 = 1 : numel(pertTypes(2 : end))
    pt = pertTypes{i2 + 1};
    set(gcf, 'CurrentAxes', hsp(i2));
    
    set(gca, 'XTick', 1 : 7 , ...
        'XTickLabel', {'s', 't-eh', 'd-iy', 'b-ae-t', 'g-ey-v', 'b-er-th', 't-uw'});
    xs = get(gca, 'XLim');
    plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    
    set(gca, 'YLim', [-20, 40]);
    xlabel('Phones');
    ylabel(sprintf('Time inteval change from %s (ms) (mean\\pm1 SEM)', baseType));
    legend(rhyConds);
   
    title(strrep(pt, '_', '\_'));
end

%% Time-intevral changes: Mark significant changes
for i1 = 1 : length(rhyConds)
    rc = rhyConds{i1};
    
    for i2 = 1 : numel(pertTypes(2 : end))
        set(gcf, 'CurrentAxes', hsp(i2));
        
        for i3 = 1 : size(p_vals_tIntChg.(rc), 1)
            if p_vals_tIntChg.(rc)(i3, i2) < P_THRESH_UNC
                plot(i3, 1e3 * mn_TIntChgs.(rc)(i3, i2), 'o', ...
                     'MarkerFaceColor', colors.(rc));
            end
            
            if p_vals_tIntChg.(rc)(i3, i2) < P_THRESH_UNC / numel(p_vals_tIntChg.(rc))
                    plot(i3, 1e3 * mn_TIntChgs.(rc)(i3, i2), 's', ...
                         'MarkerEdgeColor',  colors.(rc), 'MarkerFaceColor', 'none');
            end
            
        end
    end
end



%% Visualization: formant changes
spWspc = 0.06;
spW = 0.4;

spHspc = 0.1;
spH1 = 0.275;
spH2 = 0.1;

for h1 = 1 : numel(rhyConds)
    rc = rhyConds{h1};
    
    figure('Name', sprintf('Changes in F1 from %s under rhythm condition: %s', baseType, rc), ...
           'Position', [100, 100, 900, 600]);
    for i1 = 1 : numel(FMT_ANA_VWLS)
        t_vwl = FMT_ANA_VWLS{i1};
        
        if i1 == 1
            rowN = 2; colN = 1;
        elseif i1 == 2
            rowN = 2; colN = 2;
        elseif i1 == 3
            rowN = 1; colN = 1;
        elseif i1 == 4
            rowN = 1; colN = 2;
        else
            rowN = NaN; colN = NaN; % TODO
        end
        
%         subplot(2, 2, i1);
        hsps = [];
        hsps(end + 1) = subplot('Position', [spWspc + (colN - 1) * (spW + spWspc), spHspc + (rowN - 1) * (spH1 + spH2 + spHspc) + spH2, spW, spH1]);
        hold on;

        assert(isequal(pertTypes{1}, baseType));
%         for i2 = 2 : length(pertTypes)
        for i2 = 2 : 2
            pt = pertTypes{i2};
            tAxis = 0 : T_STEP : T_STEP * (length(avgChgVwlF1s.(t_vwl).(rc).(pt)) - 1);
            plot(tAxis * 1e3, avgChgVwlF1s.(t_vwl).(rc).(pt), ...
                 '-', 'Color', colors.(pt));
             
            
        end
        xs = get(gca, 'XLim');
        plot(xs, [0, 0], '-', 'Color', gray);

%         legend(pertTypes(2 : end), 'Location', 'Northwest');
        legend(pertTypes(2 : 2), 'Location', 'Northwest');
        title(strrep(sprintf('%s: %s', rc, t_vwl), '_', '\_'));
        
%         for i2 = 3 : length(pertTypes)
        for i2 = 2 : 2
            pt = pertTypes{i2};
            tAxis = 0 : T_STEP : T_STEP * (length(avgChgVwlF1s.(t_vwl).(rc).(pt)) - 1);
            plot(tAxis * 1e3, avgChgVwlF1s.(t_vwl).(rc).(pt) + steChgVwlF1s.(t_vwl).(rc).(pt), ...
                 '--', 'Color', colors.(pt));
            plot(tAxis * 1e3, avgChgVwlF1s.(t_vwl).(rc).(pt) - steChgVwlF1s.(t_vwl).(rc).(pt), ...
                 '--', 'Color', colors.(pt));
        end
        
        set(gca, 'XTickLabel', []);
        xs = get(gca, 'XLim');
        ylabel('F1 change (Hz)');
        draw_xy_axes;
        
        
        %--- p-value / sig-value plot ---%
        hsps(end + 1) = subplot('Position', [spWspc + (colN - 1) * (spW + spWspc), spHspc + (rowN - 1) * (spH1 + spH2 + spHspc), spW, spH2]);
        hold on;
%         for i2 = 3 : length(pertTypes)
        for i2 = 2 : 2
            pt = pertTypes{i2};
            tAxis = 0 : T_STEP : T_STEP * (length(ptChgVwlF1s.(t_vwl).(rc).(pt)) - 1);
            
            plot(tAxis * 1e3, -log10(ptChgVwlF1s.(t_vwl).(rc).(pt)));
            
            if bFmtTNorm                
                for i3 = 1 : numel(hsps)
                    set(gcf, 'CurrentAxes', hsps(i3));
                    ys = get(gca, 'YLim');
                    plot(repmat(100, 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
                    set(gca, 'YLim', ys);
                end
                
                xlabel('Time (normalized)');
            else
                xlabel('Time (ms)');
            end
            ylabel('Sig. val');
            
            plot(xs, repmat(-log10(0.05), 1, 2), '-', 'Color', [0.5, 0.5, 0.5]);
            set(gca, 'XLim', xs);
        end
        
    end
end

return