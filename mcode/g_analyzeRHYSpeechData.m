function g_analyzeRHYSpeechData(subjListFN, varargin)
%% Config
rhyConds = {'N', 'R'};
pertTypes = {'noPert', 'F1Up', 'decel'};

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

gray = [0.5, 0.5, 0.5];

P_THRESH_UNC = 0.05;

%% Read subject list
check_file(subjListFN);
subjIDs = read_struct_from_text(subjListFN);
subjIDs = splitstring(subjIDs.SUBJ_IDS, ',');
ns = numel(subjIDs);

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

bReload = ~isempty(fsic(varargin, '-r'));
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

%%
assert(isequal(pertTypes{1}, 'noPert'));

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

%% Visualization
for i0 = 1 : numel(tIntItems)
    ti = tIntItems{i0};
    
    figure('Name', sprintf('Changes in %s', ti), ...
           'Position', [100, 100, 400, 300]);
    hold on;
    for i1 = 1 : numel(rhyConds)
        rc = rhyConds{i1};

        errorbar(1 : npt - 1, 1e3 * mean(tIntChgs.(ti).(rc)), 1e3 * ste(tIntChgs.(ti).(rc)), ...
                 'o-', 'Color', colors.(rc));
    end
    legend(rhyConds);
    set(gca, 'XLim', [0, npt]);
    set(gca, 'XTick', [1, 2], 'XTickLabel', pertTypes(2 : end));
    xlim = get(gca, 'XLim');
    plot(xlim, [0, 0], '-', 'Color', gray);
    ylabel(sprintf('Change in time interval %s from noPert (ms)', ...
           strrep(ti, '_', '\_')))
end

%% 
figure('Position', [50, 150, 800, 400]);
hsp = nan(1, numel(pertTypes(2 : end)));
for i1 = 1 : numel(pertTypes(2 : end))
    hsp(i1) = subplot(1, 2, i1);
    hold on;
end

for i1 = 1 : length(rhyConds)
    rc = rhyConds{i1};
    
    mn_TIntChgs.(rc) = nan(0, 2);
    se_TIntChgs.(rc) = nan(0, 2);
    p_vals = nan(0, 2);    
    for j1 = 1 : numel(tIntItems)
        t_item = tIntItems{j1};
        
        mn_TIntChgs.(rc) = [mn_TIntChgs.(rc); mean(tIntChgs.(t_item).(rc))];
        se_TIntChgs.(rc) = [se_TIntChgs.(rc); ste(tIntChgs.(t_item).(rc))];
        
        t_p_vals = nan(1, size(tIntChgs.(t_item).(rc), 2));
        
        for j2 = 1 : size(tIntChgs.(t_item).(rc), 2)
            [~, t_p_vals(j2)] = ttest(tIntChgs.(t_item).(rc)(:, j2));
        end
        
        p_vals = [p_vals; t_p_vals];
    end
                 
    for i2 = 1 : numel(pertTypes(2 : end))
        set(gcf, 'CurrentAxes', hsp(i2));
        errorbar(1 : size(mn_TIntChgs.(rc), 1), ...
                 1e3 * mn_TIntChgs.(rc)(:, i2), ...
                 1e3 * se_TIntChgs.(rc)(:, i2), ...
                 'o-', 'Color', colors.(rc));
             
    end
    
    for i2 = 1 : numel(pertTypes(2 : end))
        for i3 = 1 : size(p_vals, 1)
            if p_vals(i3, i2) < P_THRESH_UNC
                plot(i3, 1e3 * mn_TIntChgs.(rc)(i3, i2), 'o', ...
                     'MarkerFaceColor', colors.(rc));
            end
            
            if p_vals(i3, i2) < P_THRESH_UNC / numel(p_vals)
                    plot(i3, 1e3 * mn_TIntChgs.(rc)(i3, i2), 's', ...
                         'MarkerEdgeColor',  colors.(rc), 'MarkerFaceColor', 'none');
            end
            
        end
    end
end

for i2 = 1 : numel(pertTypes(2 : end))
    pt = pertTypes{i2 + 1};
    set(gcf, 'CurrentAxes', hsp(i2));
    
    set(gca, 'XTick', 1 : 7 , ...
        'XTickLabel', {'s', 't-eh', 'd-iy', 'b-ae-t', 'g-ey-v', 'b-er-th', 't-uw'});
    xs = get(gca, 'XLim');
    plot(xs, [0, 0], '-', 'Color', [0.5, 0.5, 0.5]);
    
    set(gca, 'YLim', [-20, 30]);
    xlabel('Phones');
    ylabel('Time inteval change from noPert (ms) (mean\pm1 SEM)');
    legend(rhyConds);
   
    title(pt);
end

return