function comments_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
i1 = get(uihdls.hlist, 'Value');
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state

% if isfield(state.trialList, 'isOther') && state.trialList.isOther(i1) == 1
%     dataFld = 'otherData';
% elseif state.trialList.isRand(i1) == 1
%     dataFld = 'randData';
% elseif state.trialList.isSust(i1) == 1
%     dataFld = 'sustData';
% end
dataFld = 'mainData';

idx_trial = state.trialList.allOrderN(i1);

str = get(uihdls.edit_comments, 'String');
pdata.(dataFld).comments{idx_trial} = str;
save(dacacheFN, 'pdata');

fn = state.trialList.fn{i1};
fprintf('INFO: trial: %s: comments -> [%s]\n', fn, pdata.(dataFld).comments{idx_trial});
return