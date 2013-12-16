function corrps = perm_test_tIntChgs(tIntChgs, nPerm, varargin)
%%
bVerbose = 1;

%%
tIntNames = fields(tIntChgs);
nti = numel(tIntNames);     % Number of time intervals

rhythmConds = fields(tIntChgs.(tIntNames{1}));

nrc = length(rhythmConds); % Number of rhythm conditions
npc = size(tIntChgs.(tIntNames{1}).(rhythmConds{1}), 2); % Number of perturbation conditions
ns = size(tIntChgs.(tIntNames{1}).(rhythmConds{1}), 1); % Number of subjects

%% Iterate through all combinations of rhythm conditions and perturbation conditions
npt = nrc * npc;    % Number of permutation tests

corrps = nan(nti, nrc, npc);

for i1 = 1 : nrc
    rc = rhythmConds{i1};
    
    for i2 = 1 : npc
        if bVerbose
            info_log(sprintf('Performing permutation test on rhythm condition %s / perturbation condition %d of %d...', ...
                             rc, i2, npc));
        end
        
        %--- Create the data matrix ---%
        dat = nan(ns, nti);
        for j1 = 1 : nti
            ti = tIntNames{j1};
            
            dat(:, j1) = tIntChgs.(ti).(rc)(:, i2);
        end
        
        %--- Get the original sig-values ---%
        osig = get_sig_(dat);
        
        %--- Permutation main loop ---%
        psigs = nan(nPerm, nti);
        for j1 = 1 : nPerm
            rsgn = 2 * (rand(ns, 1) >= 0.5) - 1;
            
            pdat = dat .* repmat(rsgn, 1, nti);
            psigs(j1, :) = get_sig_(pdat);
        end
        
        max_psigs = max(abs(psigs), [], 2);
        
        t_corrps = nan(nti, 1);
        for j1 = 1 : nti
            t_corrps(j1) = numel(find(abs(osig(j1)) <= max_psigs)) / nPerm;
        end
        
        corrps(:, i1, i2) = t_corrps;
    end
end

return

%% Subroutines
function sig = get_sig_(dat)
nti = size(dat, 2);

sig = nan(1, nti);
for j1 = 1 : nti
    [~, t_p] = ttest(dat(:, j1));
    t_sgn = sign(mean(dat(:, j1)));
    sig(j1) = t_sgn * -log10(t_p);
end
return