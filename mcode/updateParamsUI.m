function updateParamsUI(uihdls, varargin)
%% 
DEFAULT_RATING = 2;
DEFAULT_ASR_OKAY= 1;
DEFAULT_OST_OKAY = 1;
DEFAULT_PERT_OKAY = 1;
DEFAULT_FLUENCY_CODE = [];

%%
red = [1, 0, 0];
green = [0, 0.5, 0];

if nargin == 2 % data
    data = varargin{1};
    
    rmsThresh = data.params.rmsThresh;
    nLPC = round(data.params.nLPC);
    fn1 = data.params.fn1;
    fn2 = data.params.fn2;
    aFact = data.params.aFact;
    bFact = data.params.bFact;
    gFact = data.params.gFact;
    bCepsLift = data.params.bCepsLift;
    cepsWinWidth = data.params.cepsWinWidth;
    
    srt_nLPCs = [];
    str_srt_nLPCs = {};
    idx_nLPC = NaN;
    
    rating = DEFAULT_RATING;
    bOSTOkay = DEFAULT_OST_OKAY;
    bASROkay = DEFAULT_ASR_OKAY;
    bPertOkay = DEFAULT_PERT_OKAY;
    
    comments = '';
    
    fluencyCode = DEFAULT_FLUENCY_CODE;
elseif nargin == 4 % pdata
    pdata = varargin{1};
    dataFld = varargin{2};
    idx = varargin{3};
    
    rmsThresh = pdata.(dataFld).rmsThresh(idx);
    nLPC = round(pdata.(dataFld).nLPC(idx));
    
    if isfield(pdata.(dataFld), 'srt_nLPCs') ...
            && ~isempty(pdata.(dataFld).srt_nLPCs{idx})
        srt_nLPCs = pdata.(dataFld).srt_nLPCs{idx};
        str_srt_nLPCs = cell(size(srt_nLPCs));
        
        for k1 = 1 : length(srt_nLPCs)
            str_srt_nLPCs{k1} = sprintf('%d', srt_nLPCs(k1));
        end
        
        idx_nLPC = find(srt_nLPCs == nLPC);        
    else
        srt_nLPCs = [];
        str_srt_nLPCs = {};
        idx_nLPC = NaN;
    end
    
    fn1 = pdata.(dataFld).fn1(idx);
    fn2 = pdata.(dataFld).fn2(idx);
    aFact = pdata.(dataFld).aFact(idx);
    bFact = pdata.(dataFld).bFact(idx);
    gFact = pdata.(dataFld).gFact(idx);
    bCepsLift = pdata.(dataFld).bCepsLift(idx);
    cepsWinWidth = pdata.(dataFld).cepsWinWidth(idx);
    
    rating = pdata.(dataFld).rating(idx);
    bOSTOkay = pdata.(dataFld).bOSTOkay(idx);
    bASROkay = pdata.(dataFld).bASROkay(idx);
    
    if ~isfield(pdata.(dataFld), 'bPertOkay')
        pdata.(dataFld).bPertOkay = ones(size(pdata.(dataFld).bASROkay));
    end
    bPertOkay = pdata.(dataFld).bPertOkay(idx);
    
    comments = pdata.(dataFld).comments{idx};
    
    fluencyCode = pdata.(dataFld).fluencyCode{idx};
else
    return;
end

% set(uihdls.edit_rmsThresh, 'String', sprintf('%.5f', rmsThresh));
set(uihdls.edit_rmsThresh, 'String', sprintf('%.8f', rmsThresh));
set(uihdls.edit_nLPC, 'String', sprintf('%d', nLPC));

if ~isempty(srt_nLPCs)
    set(uihdls.lst_srt_nLPCs, 'String', str_srt_nLPCs, 'Enable', 'on', ...
        'Value', idx_nLPC);
else
    set(uihdls.lst_srt_nLPCs, 'String', str_srt_nLPCs, 'Enable', 'off');
end

set(uihdls.edit_fn1, 'String', sprintf('%.1f', fn1));
set(uihdls.edit_fn2, 'String', sprintf('%.1f', fn2));
set(uihdls.edit_aFact, 'String', sprintf('%.1f', aFact));
set(uihdls.edit_bFact, 'String', sprintf('%.1f', bFact));
set(uihdls.edit_gFact, 'String', sprintf('%.1f', gFact));
set(uihdls.edit_bCepsLift, 'String', sprintf('%d', bCepsLift));
set(uihdls.edit_cepsWinWidth, 'String', sprintf('%d', cepsWinWidth));

items = get(uihdls.pm_rating, 'String');
for i1 = 1 : numel(items)
    if isequal(str2num(items{i1}), rating)
        break;
    end
end
set(uihdls.pm_rating, 'Value', i1);

items = get(uihdls.pm_ostOkay, 'String');
if bOSTOkay == 1
    set(uihdls.pm_ostOkay, 'Value', fsic(items, 'Good'));
else
    set(uihdls.pm_ostOkay, 'Value', fsic(items, 'Bad'));
end

items = get(uihdls.pm_asrOkay, 'String');
if bASROkay == 1
    set(uihdls.pm_asrOkay, 'Value', fsic(items, 'Good'));
else
    set(uihdls.pm_asrOkay, 'Value', fsic(items, 'Bad'));
end

set(uihdls.edit_comments, 'String', comments);

items = get(uihdls.pm_pertOkay, 'String');
if isnan(bPertOkay)
    set(uihdls.pm_pertOkay, 'Value', fsic(items, 'N/A'));
elseif bPertOkay == 1 
    set(uihdls.pm_pertOkay, 'Value', fsic(items, 'Good'));
else
    set(uihdls.pm_pertOkay, 'Value', fsic(items, 'Bad'));
end

for i1 = 1 : numel(uihdls.utterWords)
    uWord = uihdls.utterWords{i1};    
    btnName = sprintf('bt_%s', uWord);
    
    if ~isempty(find(fluencyCode == i1, 1))
        set(uihdls.(btnName), 'ForegroundColor', [1, 0, 0]);
    else
        set(uihdls.(btnName), 'ForegroundColor', [0, 1, 0]);
    end
    
end

return