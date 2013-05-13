function g_analyzeRHYSpeechData(subjListFN, varargin)
%% Config
rhyConds = {'N', 'R'};
pertTypes = {'noPert', 'F1Up', 'decel'};

colors.N = [0, 0.5, 0];
colors.R = [0, 0, 1];

gray = [0.5, 0.5, 0.5];

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

                tInts.(ti).(rc)(i1, j2) = nanmean(sres{i1}.asr_ints.(['asr_', ti]).(rc).(pt));
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


return