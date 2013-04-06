function subjectMessageDialog(title, msg, varargin)
%% Options
figW = 600;
figH = 400;

msgFontSize = 13;

bgClr = [0.6, 0.6, 0.6];

%%
ms = get(0, 'MonitorPosition');

if size(ms, 1) == 2
    scr2W = ms(2, 3) - ms(2, 1) + 1;
    scr2H = ms(2, 4) - ms(2, 2) + 1;
    scr2X0 = ms(2, 1);
    scr2Y0 = ms(2, 2);
    
    figX = scr2X0 + (scr2W - figW) * 0.5;
    figY = ms(1, 4) - ms(2, 4) + (scr2H - figH) / 2;
end

hfig = figure('Position', [figX, figY, figW, figH], ...
              'Menu', 'none', 'ToolBar', 'none', ...
              'Name', title, 'NumberTitle', 'off', ...
              'Color', bgClr);
htxt = uicontrol('Parent', hfig, 'Style', 'text', ...
                 'HorizontalAlignment', 'left', ...
                 'Units', 'normalized', 'Position', [0.1, 0.125, 0.8, 0.825], ...
                 'FontSize', msgFontSize, ...
                 'BackgroundColor', bgClr);
set(htxt, 'String', msg);

hButtOK = uicontrol('Parent', hfig, 'Style', 'Pushbutton', ...
                    'Units', 'normalized', 'Position', [0.7, 0.025, 0.2, 0.075], ...
                    'FontSize', msgFontSize, 'String', 'OK', ...
                    'BackgroundColor', bgClr, ...
                    'Callback', {@buttOK_cbk, hfig});
return

function buttOK_cbk(src, evnt, hfig)
close(hfig);
return