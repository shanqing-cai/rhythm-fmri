function show_fb(h, varargin)
% Input: h: handles
% 
%% Visualization options: 
barW = 0.6;
nBlink_rhythm = 4;
nBlink_rate = 3;
nBlink_int = 3;
nBlink_warn = 3;
blinkPeriod = 0.25; % Unit: s
blinkPeriod_warn = 0.4; % Unit: s

clrs = {[0, 0.5, 0], [0, 0, 1]};

warnFontWeight = 'Bold';
warnFontSize = 13;

%%
% --- Related to cv_ivis --- %
cv_ivis_thresh = 0.25;
min_cv_ivis = 0.05;

reci_cv_ivis_thresh = 1 / cv_ivis_thresh;
max_reci_cv_ivis = 1 / min_cv_ivis;

xlim_rhythm = [0, 3];

% --- Related to cv_ivis --- %
sylDur_LB = 0.15;
sylDur_UB = 0.5;

rate_LB = 1 / sylDur_UB;
rate_UB = 1 / sylDur_LB;

% --- Related to intensity -- %
int_LB = 60;
int_UB = 90;

xlim_rate = [0, 1];
xlim_int = [0, 1];

%% -- Show "try again" under ASR-detected speech error -- %%
% TODO

%% --- cv_ivis --- %%
udat = guidata(h.figUFBDat.fid);

set(0, 'CurrentFigure', h.figUFBDat.fid);
set(gcf, 'CurrentAxes', h.figUFBDat.axes_cv_ivis);

set(gca, 'XLim', xlim_rhythm);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'NR', 'R'});

set(gca, 'YTick', []);

hold on;

if ~isempty(fsic(h.showRhythmicityFB_phases, h.phase))
    if isfield(udat, 'crossBars')
        delete(udat.crossBars(1));
        delete(udat.crossBars(2));
        udat = rmfield(udat, 'crossBars');
    end
    
    if h.showRhythmicityFB_onlyRhythm
        xlim_rhythm = [1.5, 2.5];
        set(gca, 'XLim', xlim_rhythm);
        set(gca, 'XTick', [2], 'XTickLabel', {'R'});
    end
    
    if ~isfield(udat, 'rhythmDivBar')
        set(gca, 'YLim', [0, 1]);
        
        udat.rhythmDivBar = ...
            plot(xlim_rhythm, repmat(reci_cv_ivis_thresh / max_reci_cv_ivis, 1, 2), ...
            'k--', 'LineWidth', 2);
        udat.hRhythmTitle = ...
            title('Rhythmicity', 'FontWeight', 'Bold');
    end

    if ~isnan(h.t_cv_ivis)
        if (h.showRhythmicityFB_onlyRhythm == 1 && h.trialType == 2) || ...
                h.showRhythmicityFB_onlyRhythm == 0  
            for i1 = 1 : nBlink_rhythm
                    if i1 < nBlink_rhythm
                        t_clr = clrs{h.trialType};
                    else
                        t_clr = [0.5, 0.5, 0.5];
                    end
                    
                    hbar = plot(h.trialType + [-barW / 2, barW / 2], ...
                         repmat((1 / h.t_cv_ivis) / max_reci_cv_ivis, 1, 2), ...
                         '-', 'Color', t_clr, 'LineWidth', 2);
                    drawnow;
                    ys = get(gca, 'YLim');
                    if  (1 / h.t_cv_ivis) / max_reci_cv_ivis > ys(2)
                        ys(2) = 1.1 * (1 / h.t_cv_ivis) / max_reci_cv_ivis;
                        set(gca, 'YLim', ys);
                    end

                    if i1 < nBlink_rhythm
                        pause(blinkPeriod);
                        delete(hbar);
                        pause(blinkPeriod);
                        drawnow;  
                    end
            end
            
            if ~isfield(udat, 'hBarsRhythm');
                udat.hBarsRhythm = {[], []};
            end
            
            udat.hBarsRhythm{h.trialType}(end + 1) = hbar;        
        end
    else
        htxt = text(xlim_int(1) + 0.1 * range(xlim_int), 0.5, 'Cannot determine rhythmicity');
        pause(1);
        delete(htxt);
    end
else
    cla;
    hold on;
%     if ~isfield(udat, 'crossBars')
    udat.crossBars(1) = plot(xlim_rhythm, [0, 1], 'k-');
    udat.crossBars(2) = plot(xlim_rhythm, [1, 0], 'k-');
%     end
    if isfield(udat, 'rhythmDivBar')
        udat = rmfield(udat, {'rhythmDivBar', 'hBarsRhythm'});
    end
end


%% --- Speaking rate --- %%
set(0, 'CurrentFigure', h.figUFBDat.fid);
set(gcf, 'CurrentAxes', h.figUFBDat.axes_rate);

cla; hold on;

set(gca, 'XLim', xlim_rate);
set(gca, 'XTick', []);
set(gca, 'YLim', [0, 1]);
set(gca, 'YTick', []);

t_rate_max = 1 / h.minSylDur;
t_rate_min = 1 / h.maxSylDur;

if ~isnan(h.t_mean_ivi)
    t_rate = 1 / h.t_mean_ivi;
end
    
if ~isempty(fsic(h.showRateFB_phases, h.phase))
    udat.rateRangeBars = nan(1, 2);
    
    udat.rateRangeBars(1) = ...
        plot(xlim_rate, ...
             repmat((t_rate_min - rate_LB) / (rate_UB - rate_LB), 1, 2), ...
             'k--', 'LineWidth', 2);
    udat.rateRangeBars(2) = ...
        plot(xlim_rate, ...
             repmat((t_rate_max - rate_LB) / (rate_UB - rate_LB), 1, 2), ...
             'k--', 'LineWidth', 2);
    udat.hRateTitle = ...
        title('Speed', 'FontWeight', 'Bold');        

    text(xlim_rate(2) * 0.3, 0.05, 'SLOW');
    text(xlim_rate(2) * 0.3, 0.95, 'FAST');
    
    if ~isnan(h.t_mean_ivi)
        
        for i1 = 1 : nBlink_rate
            hbar = plot(xlim_rate(2) * 0.5 + [-barW / 2, barW / 2], ...
                        repmat((t_rate - rate_LB) / (rate_UB - rate_LB), 1, 2), ...
                        '-', 'Color', 'k', 'LineWidth', 2);
            drawnow;

            if i1 < nBlink_rate
                pause(blinkPeriod);
                delete(hbar);
                pause(blinkPeriod);
                drawnow;  
            end
        end
        
    else
        text(xlim_rate(1) + 0.1 * range(xlim_rate), 0.5, 'Cannot determine speed');
    end

    set(gca, 'XLim', xlim_rate);
else
    plot([0, 1], [0, 1], 'k-');
    plot([0, 1], [1, 0], 'k-');
end

if ~isempty(fsic(h.showRateWarn_phases, h.phase)) 
    set(gcf, 'CurrentAxes', h.figUFBDat.axes_warn);
    cla;
    if t_rate > t_rate_max
        warnMsg = 'Slower, please';
    elseif t_rate < t_rate_min
        warnMsg = 'Faster, please';
    else
        warnMsg = '';
    end

    if ~isempty(warnMsg)
        for i1 = 1 : nBlink_warn
            hwarn = text(0.2, 0.5, warnMsg, 'Color', 'm', ...
                         'FontWeight', warnFontWeight, 'FontSize', warnFontSize);

            if i1 < nBlink_rhythm
                pause(blinkPeriod_warn);
                delete(hwarn);
                pause(blinkPeriod_warn);
                drawnow;  
            end
        end
    end
end

%% --- Intensity --- %%
set(0, 'CurrentFigure', h.figUFBDat.fid);
set(gcf, 'CurrentAxes', h.figUFBDat.axes_int);

cla; hold on;

set(gca, 'XLim', xlim_int);
set(gca, 'XTick', []);
set(gca, 'YLim', [0, 1]);
set(gca, 'YTick', []);

if ~isempty(fsic(h.showIntFB_phases, h.phase))
    udat.intRangeBars = nan(1, 2);
    udat.intRangeBars(1) = ...
        plot(xlim_int, ...
             repmat((h.minVwlLevel - int_LB) / (int_UB - int_LB), 1, 2), ...
             'k--', 'LineWidth', 2);
    udat.intRangeBars(2) = ...
        plot(xlim_int, ...
             repmat((h.maxVwlLevel - int_LB) / (int_UB - int_LB), 1, 2), ...
             'k--', 'LineWidth', 2);
    udat.hRateTitle = ...
        title('Volume', 'FontWeight', 'Bold');
    % end
    
    text(xlim_int(2) * 0.3, 0.05, 'SOFT');
    text(xlim_int(2) * 0.3, 0.95, 'LOUD');

    if ~isnan(h.t_mean_vwl_lv)
        for i1 = 1 : nBlink_int
            hbar = plot(xlim_int(2) * 0.5 + [-barW / 2, barW / 2], ...
                        repmat((h.t_mean_vwl_lv - int_LB) / (int_UB - int_LB), 1, 2), ...
                        '-', 'Color', 'k', 'LineWidth', 2);
            drawnow;

            if i1 < nBlink_int
                pause(blinkPeriod);
                delete(hbar);
                pause(blinkPeriod);
                drawnow;  
            end
        end
        
        
        
    else
        text(xlim_int(1) + 0.1 * range(xlim_int), 0.5, 'Cannot determine volume');
    end

    set(gca, 'XLim', xlim_int);
else
    plot([0, 1], [0, 1], 'k-');
    plot([0, 1], [1, 0], 'k-');
end

if ~isempty(fsic(h.showIntWarn_phases, h.phase)) 
    set(gcf, 'CurrentAxes', h.figUFBDat.axes_warn);
    cla;
    if h.t_mean_vwl_lv > h.maxVwlLevel
        warnMsg = 'Softer, please';
    elseif h.t_mean_vwl_lv < h.minVwlLevel
        warnMsg = 'Louder, please';
    else
        warnMsg = '';
    end

    if ~isempty(warnMsg)
        for i1 = 1 : nBlink_warn
            hwarn = text(0.2, 0.5, warnMsg, 'Color', 'm', ...
                         'FontWeight', warnFontWeight, 'FontSize', warnFontSize);

            if i1 < nBlink_rhythm
                pause(blinkPeriod_warn);
                delete(hwarn);
                pause(blinkPeriod_warn);
                drawnow;
            end
        end
    end
end

%%
guidata(h.figUFBDat.fid, udat);

return