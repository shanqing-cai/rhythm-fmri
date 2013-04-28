function gen_vwl_fmts_all_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN);    % gives pdata
load(stateFN);      % gives state
    
bFoundMissing = 0;
flds = {'mainData'};

for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if ~isfield(pdata.(fld), 'asrTBeg')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).rating)                
        if pdata.(fld).rating(h2) == 0
            continue;
        end

%                 if isnan(pdata.(fld).sOnsetTime(h2)) || ...
%                         isnan(pdata.(fld).p2OnsetTime(h2))
        if isnan(pdata.(fld).rating(h2)) || ...
           isnan(pdata.(fld).bASROkay(h2))
            bFoundMissing = 1;
            break;
        end
        
        if pdata.(fld).bASROkay(h2) == 1 && ~isempty(find(isnan(pdata.(fld).asrTBeg(:, h2))))
            bFoundMissing = 1;
            break;
        end
    end            
end

if bFoundMissing
    msgbox('gen_vwl_fmts_all_cbk all cannot proceed at this time because rating and/or ASR are not finished yet', ...
           'gen_vwl_fmts_all_cbk proceed', 'error');
    return
end

%%
dataFld = 'mainData';

clear TransShiftMex;
i1s = get(uihdls.hlist, 'String'); % Individual trial
for n = 1 : length(i1s)    
    set(uihdls.hlist, 'Value', n);
    drawnow;
    
    idx_trial = state.trialList.allOrderN(n);
    
    if pdata.(dataFld).rating(idx_trial) == 0 || pdata.(dataFld).bDiscard(idx_trial) == 1
        fprintf(1, 'INFO: Trial #%d in the list has a rating of 0 and/or a bDiscard of 1. gen_vwl_fmts will not be extracted for this trial.\n', n);
        continue
    end
    
    if pdata.(dataFld).bASROkay(idx_trial) == 0
        fprintf(1, 'INFO: Trial #%d in the list has a bASROkay ==  0. gen_vwl_fmts will not be extracted for this trial.\n', n);
        continue;
    end

    gen_vwl_fmts_trial_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls);
end

return