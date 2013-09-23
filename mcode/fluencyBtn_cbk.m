function fluencyBtn_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
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

for i1 = 1 : numel(uihdls.utterWords)
    uWord = uihdls.utterWords{i1};    
    btnName = sprintf('bt_%s', uWord);
    
    if hObject == uihdls.(btnName)
        if isequal(get(uihdls.(btnName), 'ForegroundColor'), [0, 1, 0]);
            set(uihdls.(btnName), 'ForegroundColor', [1, 0, 0]);
        else
            set(uihdls.(btnName), 'ForegroundColor', [0, 1, 0]);
        end
    end
end

% --- Get fluency code --- %
t_fluencyCode = [];
for i1 = 1 : numel(uihdls.utterWords)
    uWord = uihdls.utterWords{i1};    
    btnName = sprintf('bt_%s', uWord);
    
    if isequal(get(uihdls.(btnName), 'ForegroundColor'), [1, 0, 0])
        t_fluencyCode(end + 1) = i1;
    end
end

str = get(uihdls.edit_comments, 'String');
pdata.(dataFld).fluencyCode{idx_trial} = t_fluencyCode;
save(dacacheFN, 'pdata');

fn = state.trialList.fn{i1};
    fprintf(1, 'INFO: trial: %s: fluencyCode = \n\t', fn);
    disp(t_fluencyCode);
    
    fprintf(1, '\n');
return