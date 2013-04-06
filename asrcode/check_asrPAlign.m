function errCode = check_asrPAlign(apa)
% Output
%       errCode
%           0 - Good (no error detected)
%           1 - Some phones are too long
%           2 - Some consecutive phones are too short
%           10 - No sil detected at the two ends
%           20 - Number of phones <= 2
%           100 - Unequal lengths of phones, tbeg and tend

%% CONFIG
MAX_PHONE_LEN = 0.36; % Unit: s
MIN_CONSEC_PHONE_LEN = 0.04; % Unit: s
MIN_CONSEC_PHONE_LEN_CNT_MAX = 4;

%% 
if length(apa.phones) <= 2
    errCode = 20;
    return
end

if length(apa.phones) ~= length(apa.tbeg) || length(apa.phones) ~= length(apa.tend)
    errCode = 100;
    return
end
    
if ~(isequal(apa.phones{1}, 'sil') && isequal(apa.phones{end}, 'sil'))
    errCode = 10;
    return
end

pdurs = apa.tend(2 : end - 1) - apa.tbeg(2 : end - 1);

if ~isempty(find(pdurs > MAX_PHONE_LEN))
    errCode = 1;
    return
end

bshort = pdurs < MIN_CONSEC_PHONE_LEN;
[shortIdx0, shortIdx1] = get_cont_stretches(bshort);
if ~isempty(shortIdx0)
    shortStretches = shortIdx1 - shortIdx0 + 1;
    if ~isempty(find(shortStretches > MIN_CONSEC_PHONE_LEN_CNT_MAX))
        errCode = 2;
        return
    end
end

errCode = 0;


return