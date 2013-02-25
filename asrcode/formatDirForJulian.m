function formatDirForJulian(inDir, outDir)
%% CONFIG: File name pattern matching 
FNWC = 'main/rep%d/trial-9-6.mat';
REPN = [1 : 20];

DICTFN = 'E:\speechres\rhythm-fmri\asrcode\cmudict.0.7a';

SAMPFREQ = 16000;
NBITS = 16;

%% 
if ~isdir(inDir)
    error('Cannot find input directory: %s', inDir);
end

if ~isdir(outDir)
    mkdir(outDir)
else
    fprintf(1, 'Directory already exists: %s\n', inDir);
end

%% Scan for mat files
matFList = {};
repn = [];
for i1 = 1 : numel(REPN)
    fp2 = sprintf(FNWC, REPN(i1));
    
    matfn = fullfile(inDir, fp2);
    if isfile(matfn)
        matFList{end + 1} = matfn;
        repn(end + 1) = i1;
    else
        fprintf(1, 'WARNING: Cannot find mat file: %s\n', matfn);
    end
end

%% Iterate through the files and do the work
wavFList = {};
sent = {};
for i1 = 1 : numel(matFList)
    load(matFList{i1});
    
    w = resample(data.signalIn, SAMPFREQ, data.params.sr);
    wavfn = fullfile(outDir, sprintf('rep_%d.wav', repn(i1)));
    wavwrite(w, SAMPFREQ, NBITS, wavfn);
    
    wavFList{end + 1} = wavfn;
    
    fprintf(1, '%s --> %s\n', matFList{i1}, wavfn);
    
    sent{end + 1} = data.params.name;
end

%% Construct .voca file 
dtxt = textread(DICTFN, '%s', 'delimiter', '\n');

vocaWords = {};
vocaPronun = {};

for i1 = 1 : numel(sent)
    t_words = upper(splitstring(sent{i1}));
    for i2 = 1 : numel(t_words)
        t_word = t_words{i2};
        t_word = strrep(strrep(strrep(t_word, ',', ''), '.', ''), ';', '');
        
        if isempty(fsic(vocaWords, t_word))
            t_phns = get_dict_pronun(t_word, dtxt);
            if isempty(t_phns)
                error('Failed to find the pronunciation of word %s in dictionary file %s', t_word, DICTFN);
            end
            
            vocaWords{end + 1} = t_word;
            vocaPronun{end + 1} = t_phns;
        end
    end
end

vocaFN = fullfile(outDir, 'gram.voca');
vocaF = fopen(vocaFN, 'wt');
fprintf(vocaF, '%% NS_B\n');
fprintf(vocaF, '<s>    sil\n\n');
fprintf(vocaF, '%% NS_E\n');
fprintf(vocaF, '</s>    sil\n\n');

for i1 = 1 : numel(vocaWords)
    fprintf(vocaF, '%% %s\n', upper(vocaWords{i1}));
    fprintf(vocaF, '%s    ', upper(vocaWords{i1}));
    for i2 = 1 : length(vocaPronun{i1})
        fprintf(vocaF, '%s', translate_phn_cmu2julius(lower(vocaPronun{i1}{i2})));
        if i2 < length(vocaPronun{i1})
            fprintf(vocaF, ' ');
        else
            fprintf(vocaF, '\n\n');
        end
    end
end
fclose(vocaF);

fprintf(1, 'Written to .voca file: %s\n', vocaFN);

%% Create grammar file
grammarFN = fullfile(outDir, 'gram.grammar');
grammarF = fopen(grammarFN, 'wt');

for i1 = 1 : numel(sent)
    t_sent = upper(deblank(strrep(sent{i1}, '.', '')));
    fprintf(grammarF, 'S : NS_B %s NS_E\n', ...
            strrep(strrep(strrep(t_sent, ',', ''), '.', ''), ';', ''));
end

fclose(grammarF);


return