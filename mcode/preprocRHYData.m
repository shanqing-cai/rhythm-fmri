function preprocRHYData(subjID, varargin)
%% Configs
hostName = getHostName;

MAIN_UTTER = 'The steady bat gave birth to pups';

dacacheDir='E:/speechres/rhythm-fmri/dacache';
rawDataDir='G:/DATA/RHYTHM-FMRI/';

mvaWinWidth=21;     % 21 * 1.333 = 28 (ms)
fineParseWin=100e-3;	% sec

ylim=[0,5000];
rmsOneSide=0.01;  % Unit: sec
rmsLBRatio=0.75;
shiraOneSide=0.025; % Unit: sec

%%
expDir=fullfile(rawDataDir,subjID);

if ~isdir(expDir)
    error('Cannot find directory %s. Terminated.', expDir);
end
if ~isfile(fullfile(expDir,'expt.mat'))
    error('Cannot find expt.mat in directory %s. Terminated.', expDir);    
end

%%
load(fullfile(expDir,'expt.mat'));  % gives expt
fprintf(1, 'Subject ID: \t%s\n', expt.subject.name);
fprintf(1, 'Subject gender: \t%s\n', expt.subject.sex);

pdata = struct;

pdata.subject = expt.subject;
pdata.mvaWinWidth = mvaWinWidth;

dacacheFN = fullfile(dacacheDir, [pdata.subject.name, '.mat']);
stateFN = fullfile(dacacheDir, [pdata.subject.name, '_state.mat']);

%% Calculate the total number of speech trials
procPhases = {'pract1', 'pract2', 'pre', 'run1', 'run2', 'run3'};
mainData = init_data('', MAIN_UTTER, procPhases, expt, rawDataDir, subjID);

pdata.mainData = mainData;

%% Build a list of all trials
if isfile(stateFN)
    fprintf('Found state.mat at %s.\n',stateFN);
    if isempty(fsic(varargin, 'phase')) && isempty(fsic(varargin, 'pdata_fixAutoRMS'))
        a = input('Resume? (0/1): ');
    else
        a = 1;
    end
    if a==1
        load(stateFN);  % gives state
        load(dacacheFN);    % gives pdata
        trialList=state.trialList;
        trialListPert=state.trialListPert;
        bNew=0;
    else
        a = input('Are you sure you want to start the screening process over? (0/1): ');
        if a == 1
            bNew=1;
        else
            return
        end
    end
else
    bNew=1;
end

if bNew
    if isfile(dacacheFN)
        delete(dacacheFN);
        fprintf('%s deleted.\n',dacacheFN);
    end
    
    fprintf('\nstate.mat not found.\n');
    fprintf('Getting information about all trials...\n');
        
    trialList.fn={};
    trialList.phase={};
    trialList.block=[];
    trialList.trialN=[];   
    trialList.word={};
    trialList.allOrderN=[];
    
    fprintf('Randomizing the order of all trials...\n');
    
    idxN = find(mainData.bRhythm ~= 1);
    idxR = find(mainData.bRhythm == 1);
    
    idxRandPerm_N = randperm(numel(idxN));
    idxRandPerm_R = randperm(numel(idxR));
    
    % --- Decide at random whether the N or R trials will be processed
    % first --- 
    if rand(1) < 0.5
        idxRandPerm = [idxN(idxRandPerm_N), idxR(idxRandPerm_R)];
    else
        idxRandPerm = [idxR(idxRandPerm_R), idxN(idxRandPerm_N)];
    end
    
    trialList.fn = mainData.rawDataFNs(idxRandPerm);
    trialList.phase = mainData.phases(idxRandPerm);
    trialList.block = mainData.blockNums(idxRandPerm);
    trialList.trialN = mainData.trialNums(idxRandPerm);
    trialList.word = mainData.words(idxRandPerm);
    trialList.bRhythm = mainData.bRhythm(idxRandPerm);
    trialList.pertType = mainData.pertType(idxRandPerm);
    
    idxOrder_o = 1 : numel(mainData.trialNums);
    trialList.allOrderN = idxOrder_o(idxRandPerm);
    
    % The perturbed trials in the rand phase
    idx_pert = find(trialList.pertType ~= 0);
    trialListPert.fn = trialList.fn(idx_pert);
    trialListPert.phase = trialList.phase(idx_pert);
    trialListPert.block = trialList.block(idx_pert);
    trialListPert.trialN = trialList.trialN(idx_pert);
    trialListPert.word = trialList.word(idx_pert);
    trialListPert.bRhythm = trialList.bRhythm(idx_pert);
    trialListPert.pertType = trialList.pertType(idx_pert);
    
    state.trialList = trialList;
    state.trialListPert = trialListPert;
    
    state.rawDataDir = rawDataDir;
    state.dacacheDir = dacacheDir;
    state.expDir = expDir;
    
    state.persist_rmsThresh = NaN;
    state.bFirstTime = 1;
    
    state.stats = zeros(1, numel(state.trialList.fn));
    state.statsPert = zeros(1, numel(state.trialListPert.fn));   
    
    save(stateFN, 'state');
    save(dacacheFN, 'pdata');
    
    fprintf('Saved state info to %s\n', stateFN);
    fprintf('Saved pdata to %s\n', dacacheFN);
else
    load(stateFN);
    load(dacacheFN);
    
    state.bFirstTime = 1;
end

%%
uihdls = struct;
uihdls.dacacheFN = dacacheFN;

uihdls.hfig = figure('Position', [20, 150, 1560, 600]);
hlist_title = uicontrol('Style', 'text', ...
                  'Unit', 'Normalized', ...
                  'Position', [0.02, 0.93, 0.18, 0.05], ...
                  'String', 'Trial list: (*: Vowel bounds done)', ...
                  'HorizontalAlignment', 'left');
uihdls.hlist_title = hlist_title;

hlist = uicontrol('Style', 'listbox', ...
                  'Unit', 'Normalized', ...
                  'Position', [0.02, 0.15, 0.15, 0.8], ...
                  'BackgroundColor', [1, 1, 1]);
uihdls.hlist = hlist;

uihdls.hmenu_nLPC = uimenu('Parent', uihdls.hfig, 'Label', 'nLPC');
set(uihdls.hfig, 'MenuBar', 'none');
uihdls.hmenu_nLPC_show_overall_best = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Show overall best');
uihdls.hmenu_nLPC_set_overall_best = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Set all trials to overall best');
uihdls.hmenu_nLPC_set_list_1st = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Set all trials to 1st in list', 'Separator', 'on');
uihdls.hmenu_nLPC_restore_user = uimenu(uihdls.hmenu_nLPC, ...
                                        'Label', 'Restore user selections', 'Separator', 'on');                                   
                                    
uihdls.hmenu_comments = uimenu('Parent', uihdls.hfig, 'Label', 'Comments');
uihdls.hmenu_comments_recover = uimenu(uihdls.hmenu_comments, ...
                                       'Label', 'Recover from file...');
                                   
uihdls.hmenu_rmsThresh = uimenu('Parent', uihdls.hfig, 'Label', 'rmsThresh');
uihdls.hmenu_rmsThresh_scan = uimenu('Parent', uihdls.hmenu_rmsThresh, ...
                                     'Label', 'Scan for trials with gaps');

hreveal = uicontrol('Style', 'pushbutton', ...
                    'Unit', 'Normalized', ...
                    'Position', [0.02, 0.04, 0.18, 0.04], ...
                    'String', 'Reveal trial details');
uihdls.hreveal = hreveal;

hShowComments = uicontrol('Style', 'pushbutton', ...
                          'Unit', 'Normalized', ...
                          'Position', [0.02, 0.09, 0.18, 0.04], ...
                          'String', 'Show comments');
uihdls.hShowComments = hShowComments;

% htitle = uicontrol('Style', 'text', ...
%                    'Unit', 'Normalized', ...
%                    'Position', [0.24, 0.93, 0.3, 0.04], ...
%                    'String', 'Title', ...
%                    'FontSize', 12);
% uihdls.htitle = htitle;

haxes1 = axes('Unit', 'Normalized', 'Position', [0.21, 0.26, 0.63, 0.68]);
uihdls.haxes1 = haxes1;

haxes2 = axes('Unit', 'Normalized', 'Position', [0.21, 0.13, 0.63, 0.15]);
uihdls.haxes2 = haxes2;

% Zoom buttons
hzo = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.24, 0.025, 0.10, 0.045], ...
                'String', 'Zoom out');
uihdls.hzo = hzo;

hzi = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.34, 0.025, 0.10, 0.045], ...
                'String', 'Zoom in');
uihdls.hzi = hzi;

hpleft = uicontrol('Style', 'pushbutton', ...
                   'Unit', 'Normalized', ...
                   'Position', [0.46, 0.025, 0.10, 0.045], ...
                   'String', 'Pan left');
uihdls.hpleft = hpleft;

hpright = uicontrol('Style', 'pushbutton', ...
                   'Unit', 'Normalized', ...
                   'Position', [0.56, 0.025, 0.10, 0.045], ...
                   'String', 'Pan right');
uihdls.hpright = hpright;

hzd = uicontrol('Style', 'pushbutton', ...
                'Unit', 'Normalized', ...
                'Position', [0.68, 0.025, 0.10, 0.045], ...
                'String', 'Default zoom');
uihdls.hzd = hzd;

bt_playSigIn = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [0.87, 0.88, 0.07, 0.04], ...
                         'String', 'Play sigIn');
uihdls.bt_playSigIn = bt_playSigIn;
bt_playSigOut = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [0.87, 0.83, 0.07, 0.04], ...
                         'String', 'Play sigOut');
uihdls.bt_playSigOut = bt_playSigOut;

rb_alwaysPlaySigIn = uicontrol('Style', 'radiobutton', ...
                               'Unit', 'Normalized', ...
                               'Position', [0.945, 0.88, 0.055, 0.04], ...
                               'String', 'Always play');
uihdls.rb_alwaysPlaySigIn = rb_alwaysPlaySigIn;

lblLeft = 0.845;
lblWidth = 0.04;
editLeft = 0.89;
editWidth = 0.04;
lbl_rmsThresh = uicontrol('Style', 'text', ...
                         'Unit', 'Normalized', ...
                         'Position', [lblLeft, 0.78, lblWidth, 0.04], ...
                         'String', 'rmsThresh: ');
uihdls.lbl_rmsThresh = lbl_rmsThresh;
edit_rmsThresh = uicontrol('Style', 'edit', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft, 0.78, editWidth, 0.04], ...
                         'String', 'rmsThresh', 'HorizontalAlignment', 'left');
uihdls.edit_rmsThresh = edit_rmsThresh;

bt_auto_rmsThresh = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [editLeft + editWidth + 0.004, 0.78, editWidth * 1.65, 0.04], ...
                             'String', 'Auto rmsThresh', 'FontSize', 8);
uihdls.bt_auto_rmsThresh = bt_auto_rmsThresh;

bt_auto_rmsThresh_all = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [editLeft + editWidth + 0.004, 0.73, editWidth * 1.65, 0.04], ...
                             'String', 'Auto rmsThresh all', 'FontSize', 8);
uihdls.bt_auto_rmsThresh_all = bt_auto_rmsThresh_all;

                     
lbl_nLPC = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.72, lblWidth, 0.04], ...
                     'String', 'nLPC: ');
uihdls.lbl_nLPC = lbl_nLPC;
edit_nLPC = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.72, editWidth, 0.04], ...
                      'String', 'nLPC', 'HorizontalAlignment', 'left');
uihdls.edit_nLPC = edit_nLPC;


bt_auto_nLPC = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.65, editWidth * 1.5, 0.04], ...
                         'String', 'Auto nLPC', ...
                         'FontSize', 8);
uihdls.bt_auto_nLPC = bt_auto_nLPC;

lbl_srt_nLPCs = uicontrol('Style', 'text', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.59, editWidth * 1.5, 0.04], ...
                         'String', 'Sorted nLPCs: ', 'HorizontalAlignment', 'left', ...
                         'FontSize', 8);
uihdls.lbl_srt_nLPCs = lbl_srt_nLPCs; 

lst_srt_nLPCs = uicontrol('Style', 'listbox', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.44, editWidth * 1.5, 0.16], ...
                         'String', {}, 'Enable', 'off', ...
                         'FontSize', 8);
uihdls.lst_srt_nLPCs = lst_srt_nLPCs;

bt_auto_nLPC_all = uicontrol('Style', 'pushbutton', ...
                         'Unit', 'Normalized', ...
                         'Position', [editLeft + editWidth + 0.004, 0.38, editWidth * 1.5, 0.04], ...
                         'String', 'Auto nLPC all', ...
                         'FontSize', 8);
uihdls.bt_auto_nLPC_all = bt_auto_nLPC_all;

% bt_best_nLPC = uicontrol('Style', 'pushbutton', ...
%                          'Unit', 'Normalized', ...
%                          'Position', [editLeft + editWidth + 0.004, 0.32, editWidth * 1.5, 0.04], ...
%                          'String', 'Overall best nLPC', ...
%                          'FontSize', 8);
% uihdls.bt_best_nLPC = bt_best_nLPC;

lbl_fn1 = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.67, lblWidth, 0.04], ...
                     'String', 'fn1: ');
uihdls.lbl_fn1 = lbl_fn1;
edit_fn1 = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.67, editWidth, 0.04], ...
                      'String', 'fn1', 'HorizontalAlignment', 'left');
uihdls.edit_fn1 = edit_fn1;
                  
lbl_fn2 = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.62, lblWidth, 0.04], ...
                     'String', 'fn2: ');
uihdls.lbl_fn2 = lbl_fn2;
edit_fn2 = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.62, editWidth, 0.04], ...
                      'String', 'fn2', 'HorizontalAlignment', 'left');
uihdls.edit_fn2 = edit_fn2;
                  
lbl_aFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.56, lblWidth, 0.04], ...
                     'String', 'aFact: ');
uihdls.lbl_aFact = lbl_aFact;
edit_aFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.56, editWidth, 0.04], ...
                      'String', 'aFact', 'HorizontalAlignment', 'left');
uihdls.edit_aFact = edit_aFact;

lbl_bFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.51, lblWidth, 0.04], ...
                     'String', 'bFact: ');
uihdls.lbl_bFact = lbl_bFact;
edit_bFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.51, editWidth, 0.04], ...
                      'String', 'bFact', 'HorizontalAlignment', 'left');
uihdls.edit_bFact = edit_bFact;
                  
lbl_gFact = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.45, lblWidth, 0.04], ...
                     'String', 'gFact: ');
uihdls.lbl_gFact = lbl_gFact;
edit_gFact = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.45, editWidth, 0.04], ...
                      'String', 'gFact', 'HorizontalAlignment', 'left');
uihdls.edit_gFact = edit_gFact;

lbl_bCepsLift = uicontrol('Style', 'text', ...
                          'Unit', 'Normalized', ...
                          'Position', [lblLeft, 0.38, lblWidth, 0.04], ...
                          'String', 'bCepsLift: ');
uihdls.lbl_bCepsLift = lbl_bCepsLift;
edit_bCepsLift = uicontrol('Style', 'edit', ...
                           'Unit', 'Normalized', ...
                           'Position', [editLeft, 0.38, editWidth, 0.04], ...
                           'String', 'bCepsLift', 'HorizontalAlignment', 'left');
uihdls.edit_bCepsLift = edit_bCepsLift;

lbl_cepsWinWidth = uicontrol('Style', 'text', ...
                          'Unit', 'Normalized', ...
                          'Position', [lblLeft, 0.32, lblWidth, 0.04], ...
                          'String', 'cepsWinWidth: ');
uihdls.lbl_cepsWinWidth = lbl_cepsWinWidth;
edit_cepsWinWidth = uicontrol('Style', 'edit', ...
                           'Unit', 'Normalized', ...
                           'Position', [editLeft, 0.32, editWidth, 0.04], ...
                           'String', 'cepsWinWidth', 'HorizontalAlignment', 'left');
uihdls.edit_cepsWinWidth = edit_cepsWinWidth;
                  
bt_reproc = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft, 0.27, 0.10, 0.04], ...
                      'String', 'Reprocess');
uihdls.bt_reproc = bt_reproc;

bt_relabel = uicontrol('Style', 'pushbutton', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft, 0.96, 0.10, 0.035], ...
                      'String', 'ReLabel');
uihdls.bt_relabel = bt_relabel;

bt_relabel_focus = uicontrol('Style', 'pushbutton', ...
                             'Unit', 'Normalized', ...
                             'Position', [lblLeft, 0.925, 0.10, 0.035], ...
                             'String', 'Focus and ReLabel');
uihdls.bt_relabel_focus = bt_relabel_focus;
                  
lbl_rating = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.22, lblWidth, 0.04], ...
                     'String', 'Prod. rating: ');
uihdls.lbl_rating = lbl_rating;
pm_rating = uicontrol('Style', 'popupmenu', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.22, lblWidth, 0.04], ...
                      'String', {'0', '1', '2'}, ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');
uihdls.pm_rating = pm_rating;

lbl_ostOkay = uicontrol('Style', 'text', ...
                        'Unit', 'Normalized', ...
                        'Position', [lblLeft, 0.18, lblWidth, 0.04], ...
                        'String', 'OST: ');
uihdls.lbl_ostOkay = lbl_ostOkay;
pm_ostOkay = uicontrol('Style', 'popupmenu', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.18, lblWidth, 0.04], ...
                      'String', {'Good', 'Bad'}, ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');
uihdls.pm_ostOkay = pm_ostOkay;

lbl_asrOkay = uicontrol('Style', 'text', ...
                        'Unit', 'Normalized', ...
                        'Position', [lblLeft, 0.14, lblWidth, 0.04], ...
                        'String', 'ASR: ');
uihdls.lbl_asrOkay = lbl_asrOkay;
pm_asrOkay = uicontrol('Style', 'popupmenu', ...
                      'Unit', 'Normalized', ...
                      'Position', [editLeft, 0.14, lblWidth, 0.04], ...
                      'String', {'Good', 'Bad'}, ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');
uihdls.pm_asrOkay = pm_asrOkay;

lbl_comments = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft, 0.11, lblWidth * 3, 0.03], ...
                     'String', 'Comments (e.g., unc {t1, d}):');
uihdls.lbl_comments = lbl_comments;
edit_comments = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft + 0.01, 0.08, 0.14, 0.03], ...
                      'String', 'Comments', ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');                  
uihdls.edit_comments = edit_comments;


lbl_fluency = uicontrol('Style', 'text', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft - lblWidth, 0.02, lblWidth, 0.09], ...
                     'String', 'Fluency comments:');
uihdls.lbl_fluency = lbl_fluency;
edit_fluency = uicontrol('Style', 'edit', ...
                      'Unit', 'Normalized', ...
                      'Position', [lblLeft + 0.01, 0.05, 0.14, 0.03], ...
                      'String', '', ...
                      'HorizontalAlignment', 'left', ...
                      'BackgroundColor', 'w');                  
uihdls.edit_fluency = edit_fluency;


% --- Buttons for labeling fluency --- %
uihdls.fluencyBtnLabel = uicontrol('Style', 'Text', ...
                                'Unit', 'Normalized', ...
                                'Position', [lblLeft + 0.11, 0.30, 0.04, 0.025], ...
                                'String', 'Fluency', ...
                                'HorizontalAlignment', 'left');     

uihdls.utterWords = splitstring(MAIN_UTTER);
for i1 = 1 : numel(uihdls.utterWords)
    uWord = uihdls.utterWords{i1};
    
    btnName = sprintf('bt_%s', uWord);
    uihdls.(btnName) = uicontrol('Style', 'pushbutton', ...
                                'Unit', 'Normalized', ...
                                'Position', [lblLeft + 0.11, 0.30 - 0.025 * i1, 0.04, 0.025], ...
                                'String', uWord, ...
                                'ForegroundColor', 'g', ...
                                'BackgroundColor', 'k', ...
                                'HorizontalAlignment', 'left');      
end
                  
bt_next = uicontrol('Style', 'pushbutton', ...
                     'Unit', 'Normalized', ...
                     'Position', [lblLeft + 0.10, 0.01, lblWidth, 0.03], ...
                     'String', 'Next', 'FontSize', 11);
uihdls.bt_next = bt_next;

%% Set upcall back functions
set(uihdls.hlist, 'Callback', {@list_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_playSigIn, 'Callback', {@playSig_cbk, dacacheFN, stateFN, uihdls, 'in'});
set(uihdls.bt_playSigOut, 'Callback', {@playSig_cbk, dacacheFN, stateFN, uihdls, 'out'});
set(uihdls.bt_reproc, 'Callback', {@reproc_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_relabel, 'Callback', {@relabel_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_relabel_focus, 'Callback', {@relabel_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.pm_rating, 'Callback', {@rating_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.pm_ostOkay, 'Callback', {@ostOkay_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.pm_asrOkay, 'Callback', {@asrOkay_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.edit_comments, 'Callback', {@comments_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.edit_fluency, 'Callback', {@fluency_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_auto_rmsThresh, 'Callback', {@auto_rmsThresh_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_auto_rmsThresh_all, 'Callback', {@auto_rmsThresh_all_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_auto_nLPC, 'Callback', {@auto_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.bt_auto_nLPC_all, 'Callback', {@auto_nLPC_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_nLPC_show_overall_best, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_set_overall_best, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_set_list_1st, 'Callback', {@set_nLPC_list_1st, dacacheFN, stateFN, uihdls});
set(uihdls.hmenu_nLPC_restore_user, 'Callback', {@restore_user_nLPC, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_rmsThresh_scan, 'Callback', {@rmsThresh_scan_cbk, dacacheFN, stateFN, uihdls});
% set(uihdls.bt_best_nLPC, 'Callback', {@best_nLPC_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.hmenu_comments_recover, 'Callback', {@recover_comments_from_file, dacacheFN, stateFN, uihdls});

set(uihdls.bt_next, 'Callback', {@next_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hreveal, 'Callback', {@reveal_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hShowComments, 'Callback', {@showComments_cbk, dacacheFN, stateFN, uihdls});


set(uihdls.lst_srt_nLPCs, 'Callback', {@lst_srt_nLPCs_cbk, dacacheFN, stateFN, uihdls});

set(uihdls.bt_reproc, 'Enable', 'off');
set(uihdls.bt_auto_rmsThresh, 'Enable', 'off');

set(uihdls.hzo, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hzi, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hpleft, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hpright, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});
set(uihdls.hzd, 'Callback', {@zoom_cbk, dacacheFN, stateFN, uihdls});

uihdls.utterWords = splitstring(MAIN_UTTER);
for i1 = 1 : numel(uihdls.utterWords)
    uWord = uihdls.utterWords{i1};
    btnName = sprintf('bt_%s', uWord);
    
    set(uihdls.(btnName), 'Callback', {@fluencyBtn_cbk, dacacheFN, stateFN, uihdls});
end

updateTrialList(state, uihdls);

return

function next_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(stateFN);
updateTrialList(state, uihdls, 'next', dacacheFN, stateFN); % gives state

return

function reveal_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
revButStr = get(uihdls.hreveal, 'String');
if isequal(revButStr, 'Reveal trial details')
    set(uihdls.hreveal, 'String', 'Hide trial details');
else
    set(uihdls.hreveal, 'String', 'Reveal trial details');
end

load(stateFN);  % gives state;
updateTrialList(state, uihdls)
return

function showComments_cbk(hObject, eventdata, dacacheFN, stateFN, uihdls)
showCommentsStr = get(uihdls.hShowComments, 'String');
if isequal(showCommentsStr, 'Show comments')
    set(uihdls.hShowComments, 'String', 'Hide comments');
else
    set(uihdls.hShowComments, 'String', 'Show comments');
end

load(stateFN);  % gives state;
updateTrialList(state, uihdls)
return

%%
function restore_user_nLPC(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % Gives pdata

flds = {'otherData', 'randData', 'sustData'};
bFoundMissing = 0;
a_nLPC_lst = [];
for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if length(pdata.(fld).rawDataFNs) == 0
        continue;
    end
    
    if ~isfield(pdata.(fld), 'srt_nLPCs')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).srt_nLPCs)
        if pdata.(fld).bDiscard(h2)
            continue;
        end

        if pdata.(fld).rating(h2) == 0
            continue;
        end

        if isempty(pdata.(fld).srt_nLPCs{h2})
            bFoundMissing = 1;
            break;
        else
            a_nLPC_lst = [a_nLPC_lst; pdata.(fld).srt_nLPCs{h2}];
        end
    end            
end

if bFoundMissing
    fprintf(1, 'WARNING: This action cannot be taken until auto nLPC has not been run on all trials.\n');
    return
end

if ~isfield(pdata, 'nLPC_status')
    fprintf(1, 'WARNING: User nLPCs have not been overwritten by best nLPC or list-first nLPCs yet.\nCannot perform this restoring at this moment.\n');
    return
else
    if isequal(pdata.nLPC_status, 'user')
        fprintf(1, 'WARNING: User nLPCs have not been overwritten by best nLPC or list-first nLPCs yet.\nCannot perform this restoring at this moment.\n');        
        return
    end
end

fields = {'randData', 'sustData'};
for i1 = 1 : numel(fields)
    fld = fields{i1};

    pdata.(fld).nLPC = pdata.(fld).user_nLPCs; % Restore user selections
end    

pdata.nLPC_status = 'user';
save(dacacheFN, 'pdata');
fprintf(1, 'Restored user nLPCs. \npdata saved to %s\n', dacacheFN);

set(uihdls.hlist, 'Enable', 'off');
lst_str = get(uihdls.hlist, 'String');
for i1 = 1 : numel(lst_str)
    set(uihdls.hlist, 'Value', i1);
    list_cbk(uihdls.hlist, [], dacacheFN, stateFN, uihdls);
    set(uihdls.hlist, 'Enable', 'off');
    drawnow;
end
set(uihdls.hlist, 'Enable', 'on');

fprintf(1, 'Restored user nLPCs. \npdata saved to %s\n', dacacheFN);
return

%%
function set_nLPC_list_1st(hObject, eventdata, dacacheFN, stateFN, uihdls)
load(dacacheFN); % Gives pdata

flds = {'otherData', 'randData', 'sustData'};
bFoundMissing = 0;
a_nLPC_lst = [];
for h1 = 1 : numel(flds)
    fld = flds{h1};
    
    if length(pdata.(fld).rawDataFNs) == 0
        continue;
    end
    
    if ~isfield(pdata.(fld), 'srt_nLPCs')
        bFoundMissing = 1;
        break;
    end

    for h2 = 1 : numel(pdata.(fld).srt_nLPCs)
        if pdata.(fld).bDiscard(h2)
            continue;
        end

        if pdata.(fld).rating(h2) == 0
            continue;
        end

        if isempty(pdata.(fld).srt_nLPCs{h2})
            bFoundMissing = 1;
            break;
        else
            a_nLPC_lst = [a_nLPC_lst; pdata.(fld).srt_nLPCs{h2}];
        end
    end            
end

if bFoundMissing
    fprintf(2, 'ERROR: Setting nLPC to 1st in list cannot be done until auto nLPC has not been run on all trials.\n');
    return
end
    
if ~isfield(pdata, 'nLPC_status')
    pdata.nLPC_status = 'user';
else
    if isequal(pdata.nLPC_status, 'list_1st')
        fprintf(1, 'Data are already set to list-first nLPCs. No changes will be made.\n');
        return;
    end
end

fields = {'randData', 'sustData'};
for i1 = 1 : numel(fields)
    fld = fields{i1};

    if isequal(pdata.nLPC_status, 'user')
        pdata.(fld).user_nLPCs = pdata.(fld).nLPC; % Backup user selections
    end

    for i2 = 1 : numel(pdata.(fld).nLPC)
        if ~isnan(pdata.(fld).nLPC)
            if isempty(pdata.(fld).srt_nLPCs{i2})
                continue;
            end
            pdata.(fld).nLPC(i2) = pdata.(fld).srt_nLPCs{i2}(1);
        end
    end
end    

pdata.nLPC_status = 'list_1st';

save(dacacheFN, 'pdata');
fprintf(1, 'Set all nLPCs to list-first. \npdata saved to %s\n', dacacheFN);

set(uihdls.hlist, 'Enable', 'off');
lst_str = get(uihdls.hlist, 'String');
for i1 = 1 : numel(lst_str)
    set(uihdls.hlist, 'Value', i1);
    list_cbk(uihdls.hlist, [], dacacheFN, stateFN, uihdls);
    set(uihdls.hlist, 'Enable', 'off');
    drawnow;
end
set(uihdls.hlist, 'Enable', 'on');
fprintf(1, 'Set all nLPCs to list-first. \npdata saved to %s\n', dacacheFN);

return