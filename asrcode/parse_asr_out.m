function pa = parse_asr_out(asrOutFN, wavFN, varargin)
aot = textread(asrOutFN, '%s', 'delimiter', '\n'); % ASR output text

bFoundWav = 0;
for i1 = 1 : length(aot)
    tline = aot{i1};
    if length(tline) < 18
        continue
    end
    
    if isequal(tline(1 : 18), 'input speechfile: ')
        if isequal(tline(end - length(wavFN) + 1 : end), wavFN)
            bFoundWav = 1;
            break;
        end
    end
end

if ~bFoundWav
    error('Cannot find the entry for wave file %s in the ASR output text', wavFN);
end


bFoundPABeg = 0; % Search for the "phoneme alignment beging" line
for i2 = 1 + 1 : length(aot)
    tline = aot{i2};
    
    if length(strfind(tline, ' phoneme alignment begin ')) == 1
        bFoundPABeg = 1;
        break;     
    end
end

if ~bFoundPABeg
    error('Cannot find phoneme alignment begin line for %s in ASR output text', wavFN);
end

bFoundPAEnd = 0;
for i3 = i2 + 1 : length(aot)
    tline = aot{i3};
    
    if length(strfind(tline, ' phoneme alignment end ')) == 1
        bFoundPAEnd = 1;
        break;     
    end
end

if ~bFoundPABeg
    error('Cannot find phoneme alignment end line for %s in ASR output text', wavFN);
end

palines = aot(i2 + 3 : i3 - 2);
nphns = length(palines);

pa = struct;
pa.nphns = nphns;
pa.phones = cell(1, nphns);
pa.tbeg = nan(1, nphns);
pa.tend = nan(1, nphns);

tStep = 1e2;
if ~isempty(fsic(varargin, 'frameLen'))
    frameLen = varargin{fsic(varargin, 'frameLen') + 1};
    tStep = 0.01 / frameLen * tStep;
end

for k1 = 1 : nphns
    t_line = palines{k1};
    t_items = splitstring(t_line);
    
    t_phn = t_items{5};
    if length(strfind(t_phn, '[')) == 1
        t_phn = t_phn(1 : strfind(t_phn, '[') - 1);
    end
    t_phn = strrep(strrep(t_phn, '{', ''), '}', '');
    
    if ~isempty(strfind(t_phn, '-'))
        t_phn = t_phn(strfind(t_phn, '-') + 1 : end);
    end
    
    if ~isempty(strfind(t_phn, '+'))
        t_phn = t_phn(1 : strfind(t_phn, '+') - 1);
    end
    
    pa.phones{k1} = t_phn;
    pa.tbeg(k1) = str2double(t_items{2}) / tStep;
    pa.tend(k1) = str2double(t_items{3}) / tStep;
end

return