function hFig = analyze_dysf_fraction(pdata, rhyConds, pertTypes, CHI_P_THRESH, colors)
%% Constants
markerVPad = 0.125;

%%
dysMat = nan(3, 2, 2); 

% Dimension 1: Three perturbation types
% Dimension 2: Column 1: fluent, Column 2: not fluent
% Dimension 3: Two rhythmic conditions {N, R} and total: {N, R, T}

bDysf = nan(1, length(pdata.mainData.fluencyCode));
for i1 = 1 : numel(pdata.mainData.fluencyCode)
    bDysf(i1) = ~isempty(pdata.mainData.fluencyCode{i1});
end

for i0 = 1 : 2
    for i1 = 1 : 3  % Three perturbation types
        for i2 = 1 : 2 % [Fluent, Dysfluent]
            dysMat(i1, i2, i0) = ...
                length(find(pdata.mainData.bRhythm == i0 - 1 & ...
                            pdata.mainData.pertType == i1 - 1 & ...
                            bDysf == i2 - 1));
        end
    end
end

dysMat(:, :, 3) = sum(dysMat(:, :, [1, 2]), 3);

hFig = figure('Name', 'Fraction of dysfluent productions');
for i1 = 1 : 3
    subplot(2, 2, i1);
    dm = dysMat(:, :, i1);
    ratioDysf = dm(:, 2) ./ sum(dm, 2);
    
    if i1 < 3
        rcond = rhyConds{i1};
        clr = colors.(rcond);
    else
        rcond = 'Total';
        clr = [0.5, 0.5, 1.0];
    end
    
    bar(ratioDysf, 'FaceColor', clr);
    hold on;
    for j1 = 1 : length(ratioDysf)
        text(j1 - 0.35, ratioDysf(j1) + 0.05, ...
             sprintf('%d/%d', dm(j1, 2), dm(j1, 1) + dm(j1, 2)), ...
             'Color', [0, 0, 0]);
    end
    
    set(gca, 'YLim', [0, 1]);
    set(gca, 'XTickLabel', pertTypes(1 : 3));
    title(rcond);
    ylabel('Fraction of dysfluent trials');
    
    ys = get(gca, 'YLim');
    
    
    for i2 = 2 : 3
        p_chi2 = chi2test(dm([1, i2], :));
        
        yBar = ys(2) - (markerVPad * 2 + (i2 - 2) * markerVPad) * range(ys);
        yBar1 = yBar - 0.025 * range(ys);
        yBarStar = yBar + 0.05 * range(ys);
        
        if p_chi2 < CHI_P_THRESH
            plot(mean([1, i2]), yBarStar, 'k*');
            fontWeight = 'bold';
        else
            fontWeight = 'normal';
        end

        plot([1, i2], repmat(yBar, 1, 2), 'k-');
        plot([1, 1], [yBar, yBar1], 'k-');
        plot([i2, i2], [yBar, yBar1], 'k-');
        text(mean([1, i2]) + 0.05 * (i2 - 1), yBarStar, sprintf('p=%f', p_chi2), ...
             'Color', 'k', 'FontWeight', fontWeight);
    end
    
end
return