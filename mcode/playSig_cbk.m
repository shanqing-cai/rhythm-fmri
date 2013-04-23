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
else
    soundsc(data.signalOut, data.params.sr);
end

return
