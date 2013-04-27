function scan_repDir(repDir)
check_dir(repDir);
dfns_1 = dir(fullfile(repDir, 'trial-*-1.mat'));
dfns_2 = dir(fullfile(repDir, 'trial-*-2.mat'));

dfns = [dfns_1; dfns_2];

tns = nan(1, length(dfns));
for i1 = 1 : numel(dfns)
    t_items = splitstring(dfns(i1).name, '-');
    tns(i1) = str2double(t_items{2});
end

[tns, idxsrt] = sort(tns, 'ascend');
dfns = dfns(idxsrt);
return