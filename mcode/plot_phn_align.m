function plot_phn_align(pa)
visTPad = 0.1; % s

% xs = get(gca, 'XLim');
ys = get(gca, 'YLim');
for k1 = 1 : pa.nphns
    plot(repmat(pa.tbeg(k1), 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(pa.tbeg(k1), ys(2) - 0.06 * range(ys), pa.phones{k1});
end
set(gca, 'YLim', ys);

if pa.nphns > 2
    xlims(1) = max([pa.tbeg(2) - visTPad, pa.tbeg(1)]);
    xlims(2) = min([pa.tend(end - 1) + visTPad, pa.tend(end)]);
else
    xlims = get(gca, 'XLim');
end
set(gca, 'XLim', xlims);
return