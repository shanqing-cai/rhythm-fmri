function [idx_on_pos, idx_off_pos, idx_on_neg, idx_off_neg] = ...
         get_sig_stretches(chgVwlF1s, P_THRESH_UNC)
glen = size(chgVwlF1s, 1);
sigs = nan(glen, 1);
for i1 = 1 : glen
    [~, p] = ttest(chgVwlF1s(i1, :));
    sigs(i1) = sign(mean(chgVwlF1s(i1, :))) * -log10(p);
end

[idx_on_pos, idx_off_pos] = get_cont_stretches(sigs > -log10(P_THRESH_UNC));
[idx_on_neg, idx_off_neg] = get_cont_stretches(sigs < log10(P_THRESH_UNC));
return