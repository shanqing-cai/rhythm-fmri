function analyze_timing_response_ratio(rhyConds, rhyConds_long, ...
                                       mn_decelPert_tShifts_t1, tIntChgs, ...
                                       colors, pltLW, fontSize)
figure;
set(gca, 'FontSize', fontSize);
hold on; box on;

%---  Ratio of compensation ---%
ratio_decel_tCompens = struct;

for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    
    ratio_decel_tCompens.(rc) = ...
        (tIntChgs.s_t1.(rc)(:, 2) + tIntChgs.t1_d.(rc)(:, 2) + tIntChgs.d_b1.(rc)(:, 2)) ./ ...
        mn_decelPert_tShifts_t1.(rc)(:);
% 	ratio_decel_tCompens.(rc) = ...
%         (tIntChgs.t1_d.(rc)(:, 2)) ./ ...
%         mn_decelPert_tShifts_t1.(rc)(:);
    
    plot(1e3 * mn_decelPert_tShifts_t1.(rc), ...
         1e3 * (tIntChgs.s_t1.(rc)(:, 2) + tIntChgs.t1_d.(rc)(:, 2) + tIntChgs.d_b1.(rc)(:, 2)), 'o', 'Color', colors.(rc));
end
xlabel('Mean duration of the syllable s+t+eh (ms)');
ylabel('Mean duration response under decel (ms)');

%--- Visualize the ratio of compensation ---%
figure('Position', [150, 150, 400, 300]);
hold on; box on;
set(gca, 'FontSize', fontSize);
for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    
    fprintf(1, 'Ratio of duration response under rhythm condition %s: mean = %f%%, SEM = %f%%\n', ...
            rc, 1e2 * mean(ratio_decel_tCompens.(rc)), 1e2 * ste(ratio_decel_tCompens.(rc)));
    bar(i1, 1e2 * mean(ratio_decel_tCompens.(rc)), ...
            'EdgeColor', 'k', 'FaceColor', colors.(rc), 'LineWidth', pltLW);
    plot([i1, i1], 1e2 * (mean(ratio_decel_tCompens.(rc)) + [-1, 1] * ste(ratio_decel_tCompens.(rc))), ...
         'Color', 'k', 'LineWidth', pltLW);
end
set(gca, 'XLim', [0.5, length(rhyConds) + 0.5], ...
         'XTick', 1 : length(rhyConds), 'XTickLabel', rhyConds_long);
ylabel('Ratio of timing response (%, mean\pm1 SEM)', 'FontSize', fontSize - 1);

%--- Comparison of compensation ratio ---%
[h, p, ~, tstats] = ttest(ratio_decel_tCompens.N, ratio_decel_tCompens.R);
xs = get(gca, 'XLim'); ys = get(gca, 'YLim');
text(xs(1) + 0.05 * range(xs), ys(2) - 0.075 * range(ys), ...
     sprintf('Paired t-test: t = %f, p = %f', tstats.tstat, p), ...
     'FontSize', fontSize - 1);

return