function show_legend_sig(x, y, w, h, ...
                         rhyConds, pltLW, colors, ...
                         corrSigSquareClr, corrSigSquareLW, corrSigSquareMarkerSize, ...
                         P_THRESH_UNC, P_THRESH_CORR, varargin)
%% 
rh2 = 0.82;
rh1 = 0.44;
rh0 = 0.18;


%%
hold on;
rectangle('Position', [x, y, w, h], 'EdgeColor', 'k');

for i1 = 1 : numel(rhyConds)
    rc = rhyConds{i1};
    plot(x + 0.1 * w + (i1 - 1) * 0.1 * w, y + rh2 * h, 'o', ...
        'MarkerEdgeColor', colors.(rc), 'MarkerFaceColor', colors.(rc));
    
    plot(x + 0.1 * w + (i1 - 1) * 0.1 * w, y + rh1 * h, 'o', ...
        'MarkerEdgeColor', colors.(rc), 'MarkerFaceColor', colors.(rc));
    plot(x + 0.1 * w + (i1 - 1) * 0.1 * w, y + rh1 * h, 's', ...
         'MarkerEdgeColor',  corrSigSquareClr, 'MarkerFaceColor', 'none', ...
         'LineWidth', corrSigSquareLW, ...
         'MarkerSize', corrSigSquareMarkerSize);
end

text(x + 0.1 * w + numel(rhyConds) * 0.1 * w, y + rh2 * h, ...
     sprintf('Uncorrected p<%.2f', P_THRESH_UNC));
text(x + 0.1 * w + numel(rhyConds) * 0.1 * w, y + rh1 * h, ...
     sprintf('Corrected p<%.2f', P_THRESH_CORR)); 
 
corrMeth = 'Permutation test';
if ~isempty(fsic(varargin, '--bonferroni'))
    corrMeth = 'Bonferroni corr.';
end
text(x + 0.1 * w + numel(rhyConds) * 0.1 * w, y + rh0 * h, ...
     sprintf('%s', corrMeth));

return