function varargout = run_julian(dataFN, varargin)
%% CONFIG
BASE_JCONF = 'julius-3.5.2-quickstart-windows/julian.jconf';
HMMDEF = 'julius-3.5.2-quickstart-windows/acoustic_model_files_build726/hmmdefs';
TIEDLIST = 'julius-3.5.2-quickstart-windows/acoustic_model_files_build726/tiedlist';
MKDFA_BIN = 'julius-3.5.2-quickstart-windows/bin/mkdfa.pl';
JULIAN_BIN = 'julius-3.5.2-quickstart-windows/bin/julian.exe';
ASR_DATA_DIR = 'asr';
% MKDFA_BIN = 'E:/speechres/rhythm-fmri/asrcode/mkdfa.pl';

GRAM_VOCA = 'asr/gram.voca';
GRAM_GRAMMAR = 'asr/gram.grammar';

DOS2UNIX_BIN = 'C:\Programs\dos2unix\bin\dos2unix.exe';
UNIX2DOS_BIN = 'C:\Programs\dos2unix\bin\unix2dos.exe';

wfs = 16e3;

%% Set all absolute paths 
cwd = pwd;
while ~(length(cwd) > 11 && isequal(cwd(end - 10 : end), 'rhythm-fmri'));
    cwd = fileparts(cwd);
end

if ~(length(cwd) > 11 && isequal(cwd(end - 10 : end), 'rhythm-fmri'))
    error('Unrecognized current working directory: %s', pwd);
end

BASE_JCONF = fullfile(cwd, BASE_JCONF);
HMMDEF = fullfile(cwd, HMMDEF);
TIEDLIST = fullfile(cwd, TIEDLIST);
MKDFA_BIN = fullfile(cwd, MKDFA_BIN);
JULIAN_BIN = fullfile(cwd, JULIAN_BIN);
ASR_DATA_DIR = fullfile(cwd, ASR_DATA_DIR);

GRAM_VOCA = fullfile(cwd, GRAM_VOCA);
GRAM_GRAMMAR = fullfile(cwd, GRAM_GRAMMAR);

check_file(JULIAN_BIN);

%% Optional input arguments
bPrep = ~isempty(fsic(varargin, 'prep'));

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
if isempty(fsic(varargin, 'output'))
    y = data.signalIn;
else
    y = data.signalOut;
end
y = resample(y, wfs, data.params.sr);
wavFN = fullfile(outDir, 'speech.wav');
wavwrite(y, wfs, wavFN);

wavFList = fullfile(outDir, 'wavlist');
wavFList_f = fopen(wavFList, 'wt');
fprintf(wavFList_f, '%s\n', wavFN);
fclose(wavFList_f);



%% Generate the trial-specific grammar file
if isempty(fsic(varargin, 'sent'))
    sent = data.params.name;
else
    sent = varargin{fsic(varargin, 'sent') + 1};
end
sent = upper(sent);
sent = strrep(sent, '''', '');
sent = strrep(sent, ',', '');
sent = strrep(sent, '.', '');
sent = strrep(sent, ';', '');
sent = strrep(sent, '-', ' ');

srcBaseFN = sprintf('sent_%s', strrep(sent, ' ', '_'));

trialGrammarFN = fullfile(outDir, 'gram.grammar');

srcGrammarFN = fullfile(ASR_DATA_DIR, [srcBaseFN, '.grammar']);
if ~isfile(srcGrammarFN)
    error('Canont find source grammar file: %s', srcGrammarFN);
end

copyfile(srcGrammarFN, trialGrammarFN);

%% Copy over gram.voca; Remove the words that do not exist in the sentence
trialVocaFN = fullfile(outDir, 'gram.voca');

srcVocaFN = fullfile(ASR_DATA_DIR, [srcBaseFN, '.voca']);
if ~isfile(srcVocaFN)
    error('Canont find source voca file: %s', srcVocaFN);
end

copyfile(srcVocaFN, trialVocaFN);

%% Run dos2unix on the gram.voca and gram.grammar files (required by the mkdfa.pl script)
% dfaFN = [GRAM_BASE, '.dfa'];
% termFN = [GRAM_BASE, '.term'];
% dictFN = [GRAM_BASE, '.dict'];
dfaFN =  fullfile(outDir, 'gram.dfa');
termFN =  fullfile(outDir, 'gram.term');
dictFN = fullfile(outDir, 'gram.dict');

srcDfaFN  = fullfile(ASR_DATA_DIR, [srcBaseFN, '.dfa']);
if ~isfile(srcDfaFN)
    error('Canont find source dfa file: %s', srcDfaFN);
end

copyfile(srcDfaFN, dfaFN);

srcTermFN  = fullfile(ASR_DATA_DIR, [srcBaseFN, '.term']);
if ~isfile(srcTermFN)
    error('Canont find source term file: %s', srcTermFN);
end

copyfile(srcTermFN, termFN);

srcDictFN  = fullfile(ASR_DATA_DIR, [srcBaseFN, '.dict']);
if ~isfile(srcDictFN)
    error('Canont find source dict file: %s', srcDictFN);
end

copyfile(srcDictFN, dictFN);

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
            
if ~isempty(fsic(varargin, 'frameLen'))
    frameLen = varargin{fsic(varargin, 'frameLen') + 1}; % Unit: s
    fshift = round(frameLen * wfs);
    jcmd1 = sprintf('%s -fshift %d', jcmd1, fshift);
end
            
if bPrep
    stdOutFN = fullfile(outDir, 'julian_stdout.txt');
    jcmd1 = [jcmd1, ' > ', stdOutFN];
    varargout{1} = jcmd1;
    
    
    %% Clean up
%     if isfile(wavFN)
%         delete(wavFN);
%     end
    
    return;
end


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
varargout{1} = pa;

if nargout >= 2
    varargout{2} = outDir;
end


%% Clean up
if isfile(wavFN)
    delete(wavFN);
end

%% Optional visualization
if ~isempty(fsic(varargin, '-s'))
    hfig = figure('Position', [50, 100, 1350, 450]);
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
    
    if nargout == 3
        varargout{3} = hfig;
    end
end

if ~isempty(fsic(varargin, '-p'))
    soundsc(y, wfs);
end
    
return