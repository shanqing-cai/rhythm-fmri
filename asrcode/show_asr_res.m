function show_asr_res(inDir, wavFN)
asrOutFN = fullfile(inDir, 'asr_out');
if ~isfile(asrOutFN)
    error('Cannot find asr_out file: %s', asrOutFN);
end

aot = textread(asrOutFN, '%s', 'delimiter', '\n'); % ASR output text

if ~isequal(wavFN(end - 3 : end), '.wav')
   wavFN = [wavFN, '.wav']; 
end
    
wavFullFN = fullfile(inDir, wavFN);
if ~isfile(wavFullFN)
    error('Cannot find wav file: %s', wavFullFN);
end

%% Parse asr output
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
for i2 = i1 + 1 : length(aot)
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
    pa.tbeg(k1) = str2double(t_items{2}) / 1e2;
    pa.tend(k1) = str2double(t_items{3}) / 1e2;
end

%% Visualization
figure('Position', [50, 100, 1350, 450]);
subplot('Position', [0.05, 0.1, 0.925, 0.8]);
[w, wfs] = wavread(wavFullFN);
show_spectrogram(w, wfs, 'noFig');
hold on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

ys = get(gca, 'YLim');
for k1 = 1 : pa.nphns
    plot(repmat(pa.tbeg(k1), 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
    text(pa.tbeg(k1), ys(2) - 0.06 * range(ys), pa.phones{k1});
end
set(gca, 'YLim', ys);

title(strrep(wavFN, '_', '\_'));
drawnow;

soundsc(w, wfs);
return