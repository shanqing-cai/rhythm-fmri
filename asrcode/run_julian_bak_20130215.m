function pa = run_julian(dataFN, varargin)
%% CONFIG
BASE_JCONF = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/julian.jconf';
HMMDEF = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/acoustic_model_files_build726/hmmdefs';
TIEDLIST = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/acoustic_model_files_build726/tiedlist';
MKDFA_BIN = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/bin/mkdfa.pl';
JULIAN_BIN = 'E:/speechres/rhythm-fmri/julius-3.5.2-quickstart-windows/bin/julian.exe';
% MKDFA_BIN = 'E:/speechres/rhythm-fmri/asrcode/mkdfa.pl';

GRAM_VOCA = 'E:/speechres/rhythm-fmri/asr/gram.voca';
GRAM_GRAMMAR = 'E:/speechres/rhythm-fmri/asr/gram.grammar';

DOS2UNIX_BIN = 'C:\Programs\dos2unix\bin\dos2unix.exe';
UNIX2DOS_BIN = 'C:\Programs\dos2unix\bin\unix2dos.exe';

wfs = 16e3;

%%
GRAM_BASE = strrep(GRAM_VOCA, '.voca', '');

if isstruct(dataFN) % Input is data, not data file name
    data = dataFN;
    
    if isempty(fsic(varargin, 'outDir'))
        error('%s: outDir must be supplied when the input is data (not data file name)');
    end
    outDir = varargin{fsic(varargin, 'outDir') + 1};
else
    if ~isfile(dataFN)
        error('Cannot find input data file: %s', dataFN);
    end
    load(dataFN);

    outDir = strrep(dataFN, '.mat', '_asr');
end

if ~isdir(outDir)
    mkdir(outDir);
end

% --- Create the wav file --- %
y = data.signalIn;
y = resample(y, wfs, data.params.sr);
wavFN = fullfile(outDir, 'speech.wav');
wavwrite(y, wfs, wavFN);

wavFList = fullfile(outDir, 'wavlist');
wavFList_f = fopen(wavFList, 'wt');
fprintf(wavFList_f, '%s\n', wavFN);
fclose(wavFList_f);

%% Generate the trial-specific grammar file
sent = upper(data.params.name);
sent = strrep(sent, '''', '');
sent = strrep(sent, ',', '');
sent = strrep(sent, '.', '');
sent = strrep(sent, ';', '');

tic;
trialGrammarFN = fullfile(outDir, 'gram.grammar');
trialGrammarF = fopen(trialGrammarFN, 'wt');
fprintf(trialGrammarF, 'S : NS_B %s NS_E\n', sent);
fclose(trialGrammarF);

if ~isfile(trialGrammarFN)
    error('Failed to generate file: %s', trialGrammarF);
end

%% Copy over gram.voca; Remove the words that do not exist in the sentence
trialVocaFN = fullfile(outDir, 'gram.voca');
% copyfile(GRAM_VOCA, trialVocaFN);

sentWords = splitstring(sent, ' ');

newVocaTxt = {};
vocaTxt = textread(GRAM_VOCA, '%s', 'delimiter', '\n');
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

trialVocaF = fopen(trialVocaFN, 'wt');
for i1 = 1 : numel(newVocaTxt)
    fprintf(trialVocaF, '%s\n', newVocaTxt{i1});
end
fclose(trialVocaF);

if ~isfile(trialVocaFN)
    error('Failed to generate file: %s', trialVocaFN);
end

%% Run dos2unix on the gram.voca and gram.grammar files (required by the mkdfa.pl script)
dfaFN = [GRAM_BASE, '.dfa'];
termFN = [GRAM_BASE, '.term'];
dictFN = [GRAM_BASE, '.dict'];

if ~isempty(fsic(varargin, '-a'))
    bRedoBase = ~isfile(dfaFN) || ~isfile(termFN) || ~isfile(dictFN);
    
    if bRedoBase
        [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, GRAM_GRAMMAR))');
        [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, GRAM_VOCA))');
    end
else
    [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, trialGrammarFN))');
    [~] = evalc('dos(sprintf(''%s %s'', DOS2UNIX_BIN, trialVocaFN))');
end

% dos(sprintf('%s %s', UNIX2DOS_BIN, trialGrammarFN));
% dos(sprintf('%s %s', UNIX2DOS_BIN, trialVocaFN));
gramFiles_time = toc;


%% Call mkdfa.pl to generate the dfa file
% perl(MKDFA_BIN, fullfile(outDir, 'gram'));
tic;
if ~isempty(fsic(varargin, '-a'))
    if bRedoBase
        [~] = evalc('system(sprintf(''perl %s %s'', MKDFA_BIN, GRAM_BASE))');
    end    
else    
    [~] = evalc('system(sprintf(''perl %s %s'', MKDFA_BIN, fullfile(outDir, ''gram'')))');
    
    dfaFN = fullfile(outDir, 'gram.dfa');
    termFN = fullfile(outDir, 'gram.term');
    dictFN = fullfile(outDir, 'gram.dict');
end
% -- Check output files -- % 


if ~isfile(dfaFN) || ~isfile(termFN) || ~isfile(dictFN)
    error('Missing output file(s) from %s', MKDFA_BIN);
end
mkdfa_time = toc;

%% Generate the jconf file
tic; 
trialJConf = fullfile(outDir, 'jconf');
trialJConf_f = fopen(trialJConf, 'wt');

bjct = textread(BASE_JCONF, '%s', 'delimiter', '\n');
for i1 = 1 : numel(bjct)
    t_line = bjct{i1};
    
    if length(t_line) > 3 && isequal(t_line(1 : 3), '-h ')
        bjct{i1} = sprintf('-h %s', HMMDEF);
    elseif length(t_line) > 7 && isequal(t_line(1 : 7), '-hlist ')
        bjct{i1} = sprintf('-hlist %s', TIEDLIST);
    elseif length(t_line) > 5 && isequal(t_line(1 : 5), '-dfa ')
        bjct{i1} = sprintf('-dfa %s', dfaFN);
    elseif length(t_line) > 3 && isequal(t_line(1 : 3), '-v ')
        bjct{i1} = sprintf('-v %s', dictFN);
    end
    
    fprintf(trialJConf_f, '%s\n', bjct{i1});
end
fclose(trialJConf_f);

if ~isfile(trialJConf)
    error('Creation of config file %s failed', trialJConf);
end

jconf_time = toc;

%% First pass of running julian
jcmd1 = sprintf('%s -input rawfile -filelist %s -C %s -palign', ...
                JULIAN_BIN, wavFList, trialJConf);
% jcmd1 = sprintf('%s -input rawfile -filelist %s -C %s -palign -1pass', ...
%                 JULIAN_BIN, wavFList, trialJConf);
tic;
[so] = evalc('system(jcmd1)');
julian_time = toc;
sout = splitstring(so);

% Other possible options: -1pass, -progout, -b 0
% -n N
% -nomultigramout

% -- Scan for missing triphones --
if ~isempty(fsic(sout, 'Missing')) && ~isempty(fsic(sout, 'phones:'))    
    iline0 = fsic(sout, 'Missing');
    iline1 = fsic(sout, '//////////////////////');
    mtt = sout(iline0 + 2 : iline1 - 1);
    
    tlt = textread(TIEDLIST, '%s', 'delimiter', '\n');
    tlt1 = {};
    
%     mtt = mtt(1 : 111);
    
    for i0 = 1 : numel(mtt)
        mt = mtt{i0};
        
        [p1, p2, p3] = get_triphone_parts(mt);
        
        cands = {};
        cscores = [];
        
        for i1 = 1 : numel(tlt)
            tline = tlt{i1};
            
            if isempty(tline)
                continue;
            end
            if isempty(deblank(tline))
                continue;
            end
            if ~isempty(strfind(tline, ' or '))
                continue;
            end
            
            [t_p1, t_p2, t_p3] = get_triphone_parts(tline);
            
            if isequal(t_p2, p2)
                cands{end + 1} = tline;
                if isequal(t_p1, p1)
                    cscores(end + 1) = 2;
                elseif isequal(t_p3, p3)
                    cscores(end + 1) = 1;
                else
                    cscores(end + 1) = 0;
                end
            end
        end
        
        if isempty(cands)
            error('Cannot find a suitable substitute for missing triphone: %s', mt);
        end
        
        [~, idx_max] = max(cscores);
        st = cands{idx_max};
        
        fprintf(1, '%s --> %s\n', mt, st);
        tlt1{end + 1} = sprintf('%s %s', mt, st);
    end
    
    tlf = fopen(TIEDLIST, 'at');
    for i1 = 1 : numel(tlt1)
        fprintf(tlf, '%s\n', tlt1{i1});    
    end
    fclose(tlf);
    
    % --- Second pass --- %
    jcmd2 = sprintf('%s -input rawfile -filelist %s -C %s -palign', ...
                    JULIAN_BIN, wavFList, trialJConf);
    [so] = evalc('system(jcmd2)');
    
    asrout = so;
else
    asrout = so;
end

asrOutFN = fullfile(outDir, 'asrout');
asrOutF = fopen(asrOutFN, 'wt');
fprintf(asrOutF, '%s', asrout);
fclose(asrOutF);

pa = parse_asr_out(asrOutFN,  wavFN);

%% Optional visualization
if ~isempty(fsic(varargin, '-s'))
    figure('Position', [50, 100, 1350, 450]);
    subplot('Position', [0.05, 0.1, 0.925, 0.8]);
    % [w, wfs] = wavread(wavFullFN);
    show_spectrogram(y, wfs, 'noFig');
    hold on;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    ys = get(gca, 'YLim');
    for k1 = 1 : pa.nphns
        plot(repmat(pa.tbeg(k1), 1, 2), ys, '-', 'Color', [0.5, 0.5, 0.5]);
        text(pa.tbeg(k1), ys(2) - 0.06 * range(ys), pa.phones{k1});
    end
    set(gca, 'YLim', ys);

    title([strrep(strrep(wavFN, '\', '\\'), '_', '\_'), sprintf(': "%s"', data.params.name)])
    drawnow;
end

if ~isempty(fsic(varargin, '-p'))
    soundsc(y, wfs);
end
    
return