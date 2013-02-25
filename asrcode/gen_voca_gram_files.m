function gen_voca_gram_files()
%% CONFIG
DICTFN = 'E:\speechres\rhythm-fmri\asrcode\cmudict.0.7a';
OUTDIR = 'E:\speechres\rhythm-fmri\asr';

%%
stcs = sents();
% stcs = sents(1 : 3); % DEBUG

dtxt = textread(DICTFN, '%s', 'delimiter', '\n');

vocaWords = {};
vocaPronun = {};

sentPhns = cell(1, numel(stcs));
vowelIndices = cell(1, numel(stcs));
for i1 = 1 : numel(stcs)
    fprintf('Analyzing sentence %d of %d: %s...\n', i1, numel(stcs), stcs{i1});
%     if ~isempty(strfind(stcs{i1}, ''''))
%         fprintf(1, 'INFO: skipping sentence: %s\n', stcs{i1});
%         continue;
%        
%     end
    stcs{i1} = strrep(stcs{i1}, '''', '');
    
    t_words = upper(splitstring(stcs{i1}));
    t_sentPhns = {};
    
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
        else
            idxw = fsic(vocaWords, t_word);
            t_phns = vocaPronun{idxw};
        end
        
        t_sentPhns = [t_sentPhns, t_phns];
    end
    
    t_isVowels = zeros(1, length(t_sentPhns));
    for i3 = 1 : numel(t_sentPhns)
        t_isVowels(i3) = is_vowel_cmu(t_sentPhns{i3});
    end
    
    sentPhns{i1} = t_sentPhns;
    vowelIndices{i1} = find(t_isVowels);
    
end

% --- Save the vowel indices information to file --- %
vowelIndicesFN = fullfile(OUTDIR, 'vowel_indices.mat');
save(vowelIndicesFN, 'stcs', 'sentPhns', 'vowelIndices');

vocaFN = fullfile(OUTDIR, 'gram.voca');
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

%% Create grammar file
grammarFN = fullfile(OUTDIR, 'gram.grammar');
grammarF = fopen(grammarFN, 'wt');

for i1 = 1 : numel(stcs)
%     if ~isempty(strfind(stcs{i1}, ''''))
%         fprintf(1, 'INFO: skipping sentence: %s', stcs{i1});
%         continue;        
%     end
    stcs{i1} = strrep(stcs{i1}, '''', '');
    
    t_sent = upper(deblank(strrep(stcs{i1}, '.', '')));
    fprintf(grammarF, 'S : NS_B %s NS_E\n', ...
            strrep(strrep(strrep(t_sent, ',', ''), '.', ''), ';', ''));
end

fclose(grammarF);

%% Create individual sentence .voca, .grammar, as well as .dict, .term and .dfa files
%% Create individual sentence .voca, .grammar, as well as .dict, .term and .dfa files
MKDFA_BIN = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/bin/mkdfa.pl';
DOS2UNIX_BIN = 'C:\Programs\dos2unix\bin\dos2unix.exe';
UNIX2DOS_BIN = 'C:\Programs\dos2unix\bin\unix2dos.exe';

for i1 = 1 : numel(stcs)
    sent = stcs{i1};
    sent = upper(strrep(strrep(sent, '''', ''), ',', ''));
    fprintf('Generating files for individual sentence %d of %d: %s\n', i1, numel(stcs), sent)

    sentBase = sprintf('sent_%s', strrep(sent, ' ', '_'));

    % --- .grammar files --- 
    sentGrammarFN = fullfile(OUTDIR, [sentBase, '.grammar']);
    sentGrammarF = fopen(sentGrammarFN, 'wt');
    fprintf(sentGrammarF, 'S : NS_B %s NS_E\n', sent);
    fclose(sentGrammarF);
    
    % --- .voca file ---
    sentWords = splitstring(sent, ' ');    
    
    sentVocaFN = fullfile(OUTDIR, [sentBase, '.voca']);
    
    newVocaTxt = {};
    vocaTxt = textread(vocaFN, '%s', 'delimiter', '\n');
    for i1 = 1 : 3 : length(vocaTxt)
        t_line = vocaTxt{i1};
        t_items = splitstring(t_line, ' ');
        t_word = t_items{2};

        if ~isempty(fsic(sentWords, t_word)) ...
                || isequal(t_word, 'NS_B') || isequal(t_word, 'NS_E')
            newVocaTxt{end + 1} = vocaTxt{i1};
            newVocaTxt{end + 1} = vocaTxt{i1 + 1};
            newVocaTxt{end + 1} = vocaTxt{i1 + 2};
        end 
    end
    
    sentVocaF = fopen(sentVocaFN, 'wt');
    for i1 = 1 : numel(newVocaTxt)
        fprintf(sentVocaF, '%s\n', newVocaTxt{i1});
    end
    fclose(sentVocaF);

    % --- Use mkdfa.pl to generate the .dfa, .term and .dict files ---
    dfaFN = fullfile(OUTDIR, [sentBase, '.dfa']);
    termFN = fullfile(OUTDIR, [sentBase, '.term']);
    dictFN =  fullfile(OUTDIR, [sentBase, '.dict']);
    
    [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, sentGrammarFN))');
    [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, sentVocaFN))');
    
    [~] = evalc('system(sprintf(''perl %s %s'', MKDFA_BIN, fullfile(OUTDIR, sentBase)))');
    
end

return

