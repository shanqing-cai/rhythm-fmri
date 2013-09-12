function playSig_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls, mode)
i1 = get(uihdls.hlist, 'Value');

% load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

hostName=deblank(getHostName);
if isequal(hostName,'smcg-w510') || isequal(hostName,'smcgw510') || isequal(hostName,'smcg_w510')
    driveLet = 'G:';
else
    error('Unsupported host: %s', hostName);
end

rfn = getRawFN_(state.rawDataDir,state.trialList.fn{i1});
if ~isequal(rfn(1 : 2), driveLet)
    rfn(1 : 2) = driveLet;
end

load(rfn); % gives data

if isequal(mode, 'in')
    soundsc(data.signalIn, data.params.sr);
elseif isequal(mode, 'out')
    soundsc(data.signalOut, data.params.sr);
elseif isequal(mode, 'in/out')
    set(uihdls.hfig_aux, 'Visible', 'on');
    set(0, 'CurrentFigure', uihdls.hfig_aux);
    
    clf;
    subplot('Position', [0.1, 0.5, 0.8, 0.4]);
    show_spectrogram(data.signalIn, data.params.sr, 'noFig');
    set(gca, 'YTick', [0 : 500 : 4000]);
    grid on;
    
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    
    subplot('Position', [0.1, 0.1, 0.8, 0.4]);
    show_spectrogram(data.signalOut, data.params.sr, 'noFig');
    set(gca, 'YTick', [0 : 500 : 4000]);
    grid on;
    
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    drawnow;
    
%     soundsc(data.signalIn, data.params.sr);
%     soundsc(data.signalOut, data.params.sr); 
end

return
