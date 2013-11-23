function fmtPermRes = perm_test_fmtChgs(chgVwlF1s, ...
                                        nPermFmt, P_THRESH_UNC, anaPerts, ...
                                        varargin)
%% Configurations
bVerbose = 1;

%% Additional options
tail = 'two';  % {'two', 'left', 'right'}
if ~isempty(fsic(varargin, '--tail'))
    tail = varargin{fsic(varargin, '--tail') + 1};
end

%%
vwls = fields(chgVwlF1s);
nvs = length(vwls); % Number of vowels

rhyConds = fields(chgVwlF1s.(vwls{1}));
nrc = length(rhyConds); % Number of rhythm conditions

%%
fmtPermRes = struct;

nTotSteps = numel(anaPerts) * nvs * nrc;
stepCnt = 0;
for i0 = 1 : numel(anaPerts)
    pt = anaPerts{i0};
    
    fmtPermRes.(pt) = struct;
    
    for i1 = 1 : nvs % Iterate through all vowels
        v = vwls{i1};
    
        fmtPermRes.(pt).(v) = struct;
    
        for i2 = 1 : nrc
            stepCnt = stepCnt + 1;
            rc = rhyConds{i2};
            
            info_log(sprintf('(Step %d/%d) Performing random permutation test on formants from: pert=%s, vowel=%s, rhythm=%s ...', ...
                             stepCnt, nTotSteps, pt, v, rc));
        
            fmtPermRes.(pt).(v).(rc) = struct;
        
            glen = find(~isnan(sum(chgVwlF1s.(v).(rc).(pt), 2)), 1, 'last');
            t_chgVwlF1s = chgVwlF1s.(v).(rc).(pt)(1 : glen, :);
            
            %--- Unperturbed results ---%
            [idx_on_pos_0, idx_off_pos_0, idx_on_neg_0, idx_off_neg_0] = ...
                get_sig_stretches(t_chgVwlF1s, P_THRESH_UNC);
            len_pos = idx_off_pos_0 - idx_on_pos_0 + ones(size(idx_off_pos_0));
            len_neg = idx_off_neg_0 - idx_on_neg_0 + ones(size(idx_off_neg_0));
            
            %--- Permute ---%
            perm_max_lens_pos = nan(1, nPermFmt);
            perm_max_lens_neg = nan(1, nPermFmt);
            for k1 = 1 : nPermFmt
                rnd_sign = sign(rand(1, size(t_chgVwlF1s, 2)) - 0.5);
                perm_chgVwlF1s = t_chgVwlF1s .* repmat(rnd_sign, size(t_chgVwlF1s, 1), 1);
                
                [idx_on_pos, idx_off_pos, idx_on_neg, idx_off_neg] = ...
                    get_sig_stretches(perm_chgVwlF1s, P_THRESH_UNC);
                
                if ~isempty(idx_on_pos)
                    perm_max_lens_pos(k1) = max(idx_off_pos - idx_on_pos + ones(size(idx_off_pos)));
                else
                    perm_max_lens_pos(k1) = 0;
                end
                
                if ~isempty(idx_on_neg)
                    perm_max_lens_neg(k1) = max(idx_off_neg - idx_on_neg + ones(size(idx_off_neg)));
                else
                    perm_max_lens_neg(k1) = 0;
                end
                
            end
        
            if isequal(tail, 'two')
                perm_max_lens = max([perm_max_lens_pos; perm_max_lens_neg], [], 1);
            elseif isequal(tail, 'left')
                perm_max_lens = perm_max_lens_neg;
            elseif isequal(tail, 'right')
                perm_max_lens = perm_max_lens_pos;
            else
                error_log('Unrecognized tail: %s', tail);
            end
            
            ions = [idx_on_pos_0, idx_on_neg_0];
            ioffs = [idx_off_pos_0, idx_off_neg_0];
            lens = [len_pos, len_neg];
            sgns = [ones(size(len_pos)), -ones(size(len_neg))];
            if isequal(tail, 'left')
                ions = ions(sgns < 0);
                ioffs = ioffs(sgns < 0);
                lens = lens(sgns < 0);
                sgns = sgns(sgns < 0);
            elseif isequal(tail, 'right')
                ions = ions(sgns < 0);
                ioffs = ioffs(sgns < 0);
                lens = lens(sgns > 0);
                sgns = sgns(sgns > 0);
            end
            corrps = nan(size(sgns));
            
            for k1 = 1 : numel(lens)
                corrps(k1) = numel(find(lens(k1) < perm_max_lens)) / nPermFmt;
            end
            
            %--- Furnish output results ---
            fmtPermRes.(pt).(v).(rc).ions = ions;
            fmtPermRes.(pt).(v).(rc).ioffs = ioffs;
            fmtPermRes.(pt).(v).(rc).lens = lens;
            fmtPermRes.(pt).(v).(rc).sgns = sgns;
            fmtPermRes.(pt).(v).(rc).corrps = corrps;
            
        end
    end
end

return