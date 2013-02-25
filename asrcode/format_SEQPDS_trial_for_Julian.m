function format_SEQPDS_trial_for_Julian(inFile, inWord, inPhones, outDir, varargin)
%% CONFIG: File name pattern matching 
% FNWC = 'main/rep%d/trial-9-6.mat';
% REPN = [1 : 20];

DICTFN = 'E:\speechres\rhythm-fmri\asrcode\cmudict.0.7a';

SAMPFREQ = 16000;
NBITS = 16;

%% 
if ~isfile(inFile)
    error('Cannot find input file: %s', inFile);
end

if ~isdir(outDir)
    mkdir(outDir)
else
    fprintf(1, 'Directory already exists: %s\n', outDir);
end

%% Iterate through the files and do the work
wavFList = {};

[x, wfs] = wavread(inFile);

if ~isempty(fsic(varargin, 'wiener'))
    x = WienerScalart96(x, wfs);
end

w = resample(x, SAMPFREQ, wfs);

[tpath, tfn] = fileparts(inFile);
if isempty(strfind(tfn, '.wav'))
    tfn = [tfn, '.wav'];
end
wavfn = fullfile(outDir, tfn);
wavwrite(w, SAMPFREQ, NBITS, wavfn);

wavFList{end + 1} = wavfn;

fprintf(1, '%s --> %s\n', inFile, wavfn);        


%% Construct .voca file 
% dtxt = textread(DICTFN, '%s', 'delimiter', '\n');

vocaWords = {inWord};
vocaPronun = {inPhones};
sent = {inWord};

% Make an erroneous word: reverse the final two phones
% vocaWords{end + 1} = [inWord, '_W1'];
% vocaPronun{end + 1} = [inPhones(1 : end - 2), inPhones(end), inPhones(end - 1)];
% sent{end + 1} = [inWord, '_W1'];

% for i1 = 1 : numel(sent)
%     t_words = upper(splitstring(sent{i1}));
%     for i2 = 1 : numel(t_words)
%         t_word = t_words{i2};
%         
%         if isempty(fsic(vocaWords, t_word))
%             t_phns = get_dict_pronun(t_word, dtxt);
%             if isempty(t_phns)
%                 error('Failed to find the pronunciation of word %s in dictionary file %s', t_word, DICTFN);
%             end
%             
%             vocaWords{end + 1} = t_word;
%             vocaPronun{end + 1} = t_phns;
%         end
%     end
% end

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
    fprintf(grammarF, 'S : NS_B %s NS_E\n', t_sent);
end

fclose(grammarF);


return