function g_demographicAnalysis(demoXls, mode)
% Input arguments:
%   demoXls: File name of the xls file that contains the demographic
%            information
%   mode:    Mode of data analysis: {'behav', 'MRI', 'MRI_AWSbehav'}
%% Constants
titleStr_SID = 'Subject ID';
titleStr_DoB = 'date of birth';
titleStr_DoBS = 'date of behavioral session';
titleStr_DoMS = 'date of MRI session';

DAYS_IN_YEAR = 365.2442;

grps = {'ANS', 'AWS'};
colors.ANS = 'b';
colors.AWS = 'r';

%% 
check_file(demoXls);

%%
[N, T] = xlsread(demoXls);

%% Get the column indicies
titles = T(1, :);
for i1 = 1 : numel(titles)
    titles{i1} = lower(titles{i1});
end

c_SID = get_column_index(titles, titleStr_SID);
c_DoB = get_column_index(titles, titleStr_DoB);
c_DoBS = get_column_index(titles, titleStr_DoBS);
c_DoMS = get_column_index(titles, titleStr_DoMS);

%% 
T = T(2 : end, :);
a_SID = T(:, c_SID);
a_DoB = T(:, c_DoB);    % Date of birth
a_DoBS = T(:, c_DoBS);  % Date of behavioral session
a_DoMS = T(:, c_DoMS);  % Date of MRI session

%-- Categorize the subjects by group --%
a_grps = cell(1, length(a_SID));
for i1 = 1 : numel(a_SID)
    if ~isempty(strfind(a_SID{i1}, 'ANS_'))
        a_grps{i1} = 'ANS';
    else
        a_grps{i1} = 'AWS';
    end
end

a_DoB_nums = nan(1, length(a_DoB));
a_DoBS_nums = nan(1, length(a_DoB));
a_DoMS_nums = nan(1, length(a_DoB));
for i1 = 1 : numel(a_DoB)
    a_DoB_nums(i1) = datenum(a_DoB{i1});
    a_DoBS_nums(i1) = datenum(a_DoBS{i1});
    
    if ~isempty(a_DoMS{i1})
        a_DoMS_nums(i1) = datenum(a_DoMS{i1});
    end
        
end

a_ageBS = (a_DoBS_nums - a_DoB_nums) / DAYS_IN_YEAR;
a_ageMS = (a_DoMS_nums - a_DoB_nums) / DAYS_IN_YEAR;

%% Group-wise data
for i1 = 1 : numel(grps)
    grp = grps{i1};
    idxGrp = strmatch(grp, a_grps, 'exact');
    
    SID.(grp) = a_SID(idxGrp);
    ageBS.(grp) = a_ageBS(idxGrp);
    ageMS.(grp) = a_ageMS(idxGrp);
end

%% Visualization
grpKeep = struct;
if isequal(mode, 'behav')    
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        grpKeep.(grp) = find(~isnan(ageBS.(grp)));
        
        anaAge.(grp) = ageBS.(grp)(grpKeep.(grp));
    end
elseif isequal(lower(mode), 'mri')
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        grpKeep.(grp) = find(~isnan(ageMS.(grp)));
        
        anaAge.(grp) = ageMS.(grp)(grpKeep.(grp));
    end
elseif isequal(lower(mode), lower('MRI_AWSbehav'))
    for i1 = 1 : numel(grps)
        grp = grps{i1};
        
        if isequal(grp, 'AWS')
            grpKeep.(grp) = find(~isnan(ageBS.(grp)));
            anaAge.(grp) = ageBS.(grp)(grpKeep.(grp));
        else
            grpKeep.(grp) = find(~isnan(ageMS.(grp)));
            anaAge.(grp) = ageMS.(grp)(grpKeep.(grp));
        end
        
        
    end
else
    error_log(sprintf('Unrecognized mode: %s',  mode));
end

figure;
hold on;
for i1 = 1 : numel(grps)
    grp = grps{i1};

    plot(repmat(i1, 1, length(anaAge.(grp))), anaAge.(grp), ...
         'o', 'Color', colors.(grp));
    plot(i1 + 0.1, mean(anaAge.(grp)), 'o', 'Color', colors.(grp));
    plot(repmat(i1 + 0.1, 1, 2), mean(anaAge.(grp)) + [-1, 1] * std(anaAge.(grp)), ...
         '-', 'Color', colors.(grp));  
end
set(gca, 'XLim', [0, length(grps) + 1]);
set(gca, 'XTick', [1 : length(grps)]);
set(gca, 'XTickLabel', grps);
xlabel('Group');
ylabel('Age (y.o.)');

%--- Perform statistical comparisons ---%
[p_rs, h_rs] = ranksum(anaAge.(grps{1}), anaAge.(grps{2}));
[h_t, p_t] = ttest2(anaAge.(grps{1}), anaAge.(grps{2}));
[h_ks, p_ks] = kstest2(anaAge.(grps{1}), anaAge.(grps{2}))

xs = get(gca, 'XLim'); 
ys = get(gca, 'YLim');
for i1 = 1 : numel(grps)
    grp = grps{i1};
    text(xs(1) + 0.05 * range(xs), ys(2) - 0.06 * i1 * range(ys), ...
         sprintf('Group %s: N = %d', grp, length(anaAge.(grp))), ...
         'Color', colors.(grp));
end
text(xs(1) + 0.05 * range(xs), ys(2) - 0.18 * range(ys), ...
     sprintf('Rank-sum test: p = %f', p_rs));
text(xs(1) + 0.05 * range(xs), ys(2) - 0.24 * range(ys), ...
     sprintf('Two-sample t-test: p = %f', p_t));
text(xs(1) + 0.05 * range(xs), ys(2) - 0.30 * range(ys), ...
     sprintf('K-S test: p = %f', p_ks));

return

%% Subroutines
function cn = get_column_index(titles, titleStr)
cn = fsic(titles, lower(titleStr));

if length(cn) ~= 1
    error_log(sprintf('Cannot find exactly one entry for %s in the title column', titleStr));
end
return