function meta_rc_pt_plot(meas, measName, colors)
figure('Name', measName);

rhyConds = fields(meas);
pertTypes = fields(meas.(rhyConds{1}));

for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    errorbar([length(pertTypes) * (i1 - 1) + 1 : length(pertTypes) * i1], ...
             [mean(meas.(rc).noPert), mean(meas.(rc).F1Up), mean(meas.(rc).decel)], ...
             [ste(meas.(rc).noPert), ste(meas.(rc).F1Up), ste(meas.(rc).decel)], ...
             'o-', 'Color', colors.(rc));
	hold on;
    
    [~, p_decel_v_noPert] = ttest2(meas.(rc).decel, meas.(rc).noPert);
    fprintf(1, '%s: rc = %s: p_decel_v_noPert = %f\n', ...
            measName, rc, p_decel_v_noPert);
end
return