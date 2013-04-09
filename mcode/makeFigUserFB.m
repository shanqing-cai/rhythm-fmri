function figDat = makeFigUserFB()
%% CONFIG
bot = 0.2;
hgt = 0.7;


%%
figDat = struct;

figDat.fid = figure('Position', [200, 200, 320, 240], ...
                    'ToolBar', 'none', 'MenuBar', 'none', ...
                    'Name', 'Subject feedback', 'NumberTitle', 'off');

% cv_ivis feedback
figDat.axes_cv_ivis = subplot('Position', [0.1, bot, 0.275, hgt]);
set(figDat.axes_cv_ivis, 'XTick', [], 'YTick', []);
box on;
title('Rhythmicity');

% Speaking rate feedback
figDat.axes_rate = subplot('Position', [0.375, bot, 0.275, hgt]);
set(figDat.axes_rate, 'XTick', [], 'YTick', []);
box on;
title('Speed');

% Intensity feedback
figDat.axes_int = subplot('Position', [0.65, bot, 0.275, hgt]);
set(figDat.axes_int, 'XTick', [], 'YTick', []);
box on;
title('Intensity');

% Axes for warning message display
figDat.axes_warn = subplot('Position', [0.2, 0.025, 0.6, 0.1], 'Color', [0.25, 0.25, 0.25]);
box on;
set(gca, 'XTick', [], 'YTick', []);
return