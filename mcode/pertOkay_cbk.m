function pertOkay_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

dataFld = 'mainData';

idx_trial = state.trialList.allOrderN(i1);

val = get(uihdls.pm_pertOkay, 'Value');
items = get(uihdls.pm_pertOkay, 'String');

if isequal(items{val}, 'N/A')
    pdata.(dataFld).bPertOkay(idx_trial) = NaN;
else
    pdata.(dataFld).bPertOkay(idx_trial) = isequal(items{val}, 'Good');
end

save(dacacheFN, 'pdata');
fprintf('Saved to %s\n', dacacheFN);

fn = state.trialList.fn{i1};
fprintf('INFO: trial: %s: bPertOkay -> %d\n', fn, pdata.(dataFld).bASROkay(idx_trial));
return