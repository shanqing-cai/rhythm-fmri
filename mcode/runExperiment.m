function varargout = runExperiment(configFN, varargin)
%% CONFIG
DEBUG=0;

colors.rhythm = [0, 0, 1];
colors.nonRhythm = [0, 0.4, 0];

sentMatFN = '../asr/vowel_indices.mat';
check_file(sentMatFN);

%% ---- Load config ---- %%
subject = read_subject_config(configFN);

if isequal(getHostName, 'smcg_w510')
%     subject.dataDir        	= 'E:\DATA\RHYTHM-FMRI\';
    subject.dataDir        	= 'G:\DATA\RHYTHM-FMRI\';
else
    subject.dataDir         = 'D:\CS_2004\PROJECTS\RHYTHM-FMRI\';
end

if subject.trigByScanner == 0
    check_file(subject.basePCF);
end

%% --- In case of fMRI experiment, force input of mouthMicDist ---
if subject.trigByScanner
    if ~isempty(fsic(varargin, 'mouthMicDist'))
        subject.mouthMicDist = varargin{fsic(varargin, 'mouthMicDist') + 1};
    else
        subject.mouthMicDist = input('Please input mouth-mic distance (cm): ');
    end
end

modify_expt_config(configFN, 'mouthMicDist', subject.mouthMicDist);
fprintf(1, 'INFO: Modified the mouthMicDist field config file: %s\n', configFN);

%% Check current working directory
cwd = pwd;
[~, cwd0] = fileparts(cwd);
if ~isequal(cwd0, 'mcode')
    error('Not currently under directry "mcode"');
end

%%
subject.date				=clock;

if (~isempty(findStringInCell(varargin, 'subject')))
    clear('subject');
    subject = varargin{findStringInCell(varargin, 'subject')+1};
end

% if (~isfield(subject,'pcrKnob'))
%     subject.pcrKnob=input('Phone / Ctrl Room knob = ');
% end
subject.pcrKnob=0;

%%
bNew=true;

dirname=fullfile(subject.dataDir,num2str(subject.name));

if (~isempty(findStringInCell(varargin,'dirname')))
    clear('dirname');
    dirname=varargin{findStringInCell(varargin,'dirname')+1};
end

if isdir(dirname)
    if ~isempty(fsic(varargin, 'confirmOverwrite'))
        bNew = true;
    else
        messg={sprintf('The specified directory "%s" already contains a previously recorded experiment', dirname)
            ''
            'Continue experiment, overwrite  or cancel ?'};
        button1 = questdlg(messg,'DIRECTORY NOT EMPTY','Continue','Overwrite','Cancel','Continue');
        switch button1
            case 'Overwrite'
                button2 = questdlg({sprintf('Are you sure you want to overwrite data in directory %s?', dirname)} ,'OVERWRITE EXPERIMENT ?');
                switch button2
                    case 'Yes',
                        rmdir(dirname,'s')
                    otherwise,
                        return
                end
            case 'Continue'
                bNew=false;

            otherwise,
                return

        end
    end
end

expt.subject=subject;
if bNew % set up new experiment
    mkdir(dirname)
    
    if expt.subject.trigByScanner    % fMRI experiment
        expt.allPhases={'pre', 'pract1', 'pract2', ...
                        'run1', 'inter1', 'run2', 'inter2', ...
                        'run3', 'inter3', 'run4', 'inter4', ...
                        'run5', 'inter5', 'run6', 'inter6'};       
        expt.recPhases={'pre', 'pract1', 'pract2', ...
                        'run1', 'inter1', 'run2', 'inter2', ...
                        'run3', 'inter3', 'run4', 'inter4', ...
                        'run5', 'inter5', 'run6', 'inter6'}; %SC The pahses during which the data are recorded
        expt.skippablePhases = {'pract1', 'pract2', 'inter1', 'inter2', 'inter3', 'inter4', 'inter5', 'inter6'};
                    
        
        expt.script.pre.nReps = expt.subject.NREPS_PRE;    %SC Numbers of repetitions in the stages   % !!1!!	
        expt.script.pract1.nReps = expt.subject.NREPS_PRACT1;
        expt.script.pract2.nReps = expt.subject.NREPS_PRACT2;
        
        expt.script.run1.nReps = expt.subject.NREPS_RUN;  %SC Default 10   %SC-Mod(09/26/2007)       % !!8!!
        expt.script.inter1.nReps = expt.subject.NREPS_INTER;
        expt.script.run2.nReps = expt.subject.NREPS_RUN;   %SC Default 15   %SC-Mod(09/26/2007)       % !!2!!
        expt.script.inter2.nReps = expt.subject.NREPS_INTER;
        expt.script.run3.nReps = expt.subject.NREPS_RUN;   %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!
        expt.script.inter3.nReps = expt.subject.NREPS_INTER;
        expt.script.run4.nReps = expt.subject.NREPS_RUN;    %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!
        expt.script.inter4.nReps = expt.subject.NREPS_INTER;
        expt.script.run5.nReps = expt.subject.NREPS_RUN;
        expt.script.inter5.nReps = expt.subject.NREPS_INTER;
        expt.script.run6.nReps = expt.subject.NREPS_RUN;
        expt.script.inter6.nReps = expt.subject.NREPS_INTER;
    else
%         expt.allPhases={'pract1', 'pract2', 'pre', 'run1', 'run2', 'run3', 'run4'};
%         expt.recPhases={'pract1', 'pract2', 'pre', 'run1', 'run2', 'run3', 'run4'}; %SC The pahses during which the data are recorded
        expt.allPhases={'pract1', 'pract2', 'pre', ...
                        'run1', 'inter1', 'run2', 'inter2', 'run3'};
        expt.recPhases={'pract1', 'pract2', 'pre', ...
                        'run1', 'inter1', 'run2', 'inter2', 'run3'}; %SC The pahses during which the data are recorded
        expt.skippablePhases = {'inter1', 'inter2'};
        
        expt.script.pract1.nReps = expt.subject.NREPS_PRACT1;
        expt.script.pract2.nReps = expt.subject.NREPS_PRACT2;
        expt.script.pre.nReps = expt.subject.NREPS_PRE;
        expt.script.run1.nReps = expt.subject.NREPS_RUN;
        expt.script.inter1.nReps = expt.subject.NREPS_INTER;
        expt.script.run2.nReps = expt.subject.NREPS_RUN;
        expt.script.inter2.nReps = expt.subject.NREPS_INTER;
        expt.script.run3.nReps = expt.subject.NREPS_RUN;
    end
    
%     expt.trialTypes=[1, 2, 3, 4];  % 1: non-rhythmic speech, 2: rhythmic speech, 3: non-rhythmic baseline, 4: rhythmic baseline. 
    expt.trialTypes = eval(expt.subject.trialTypes);
    expt.trialOrderRandReps=1;	%How many reps are randomized together
    
	expt.trialTypeDesc=cell(1,5);
	expt.trialTypeDesc{1}='Non-rhythmic speech';
	expt.trialTypeDesc{2}='Rhythmic speech';
	expt.trialTypeDesc{3}='Non-rhythmic baseline';
	expt.trialTypeDesc{4}='Rhythmic baseline';
    
    nSpeechTrialsPerRep = numel(find(expt.trialTypes <= 2));
    if expt.subject.trigByScanner
        nSents=(expt.script.pre.nReps+expt.script.run1.nReps+expt.script.run2.nReps+expt.script.run3.nReps+...
                expt.script.run4.nReps+expt.script.run5.nReps+expt.script.run6.nReps) * nSpeechTrialsPerRep;
    else
        nSents=(expt.script.pract1.nReps+expt.script.pract2.nReps+expt.script.pre.nReps+...
                expt.script.run1.nReps+expt.script.run2.nReps+expt.script.run3.nReps+...
                expt.script.inter1.nReps + expt.script.inter2.nReps) * nSpeechTrialsPerRep;
    end
    
    if expt.subject.trigByScanner == 0
        [expt.stimSents_all, expt.nSyls_all] = getRandSentences(nSents);
    else
        [expt.stimSents_all, expt.nSyls_all] = getRandSentences_fmri(expt);
    end
    
    sentCnt = 1;

	nPhases=length(expt.allPhases);
    sentCnt=1;
    for i1=1:nPhases
        t_nSents=expt.script.(expt.allPhases{i1}).nReps * nSpeechTrialsPerRep;
        expt.stimSents.(expt.allPhases{i1})=expt.stimSents_all(sentCnt:sentCnt+t_nSents-1);
        expt.stimSents_nSyls.(expt.allPhases{i1})=expt.nSyls_all(sentCnt:sentCnt+t_nSents-1);
        sentCnt=sentCnt+t_nSents;
    end
    
    if expt.subject.trigByScanner
        for k1 = 1 : numel(expt.allPhases)
            ph = expt.allPhases{k1};
            expt.script.(ph) = ...
                genPhaseScript_fmri(ph, expt.script.(ph).nReps,  ...
                                    expt.trialTypes, expt.stimSents.(ph), ...
                                    expt.stimSents_nSyls.(ph), expt.trialOrderRandReps, expt.subject.trigByScanner);
        end
    else
        for k1 = 1 : numel(expt.allPhases)
            ph = expt.allPhases{k1};
            expt.script.(ph) = ...
                genPhaseScript_behav(ph, expt.script.(ph).nReps, expt.subject);
        end
    end
    
    p = getTSMDefaultParams(subject.sex,'closedLoopGain', expt.subject.closedLoopGain,...
        'trialLen', expt.subject.trialLen, ...
        'mouthMicDist', expt.subject.mouthMicDist);
    
    state.phase=1;
    state.rep=1;
    state.params=p;
    rmsPeaks=[];    
    
    save(fullfile(dirname,'expt.mat'),'expt');
    save(fullfile(dirname,'state.mat'),'state');
    
    % -- Create data structure for holding IVI (inter-vowel interval)
    % related data -- %
    timingDat = init_timingDat(expt.script);
    
    recMeanSylDurs.nonRhythm=[];
    recMeanSylDurs.rhythm=[];
    recMeanPeakRMS.nonRhythm=[];
    recMeanPeakRMS.rhythm=[];
    recCV_IVI.nonRhythm = [];
    recCV_IVI.rhythm = [];

    adaptMeanSylDur=expt.subject.paceStim.meanSylDur;
else % load expt
    load(fullfile(dirname,'state.mat'));
    load(fullfile(dirname,'expt.mat'));            
    p=state.params;
%     nPeaks=length(expt.trainWords);
%     if state.phase > 1
%         rmsPeaks=ones(length(expt.trainWords),1)*p.rmsMeanPeak;    %SC ***Bug!!***
%     end
    rmsPeaks=[];
    subject=expt.subject;
    
    if isfile(fullfile(dirname, 'timingDat.mat'))
        load(fullfile(dirname, 'timingDat.mat'));
    else
        timingDat = init_timingDat(expt.script);
        
        recMeanSylDurs.nonRhythm=[];
        recMeanSylDurs.rhythm=[];
        recMeanPeakRMS.nonRhythm=[];
        recMeanPeakRMS.rhythm=[];
        recCV_IVI.nonRhythm = [];
        recCV_IVI.rhythm = [];

        adaptMeanSylDur=expt.subject.paceStim.meanSylDur;
    end
end
trackMeanSylDurs = [repmat(expt.subject.paceStim.meanSylDur,1,4)];

% Make a copy of the configuration file
dconfig = dir(fullfile(dirname, 'config*.txt'));
configBackupFN = fullfile(dirname, sprintf('config_%.2d.txt', length(dconfig) + 1));
copyfile(configFN, configBackupFN);

% --- Optional N-R color reversal --- %
if expt.subject.colorReverseNR == 1
    clrTmp = colors.rhythm;
    colors.rhythm = colors.nonRhythm;
    colors.nonRhythm = clrTmp;
end

%% Determine the log file name
dlogs = dir(fullfile(dirname, 'log_*.txt'));
logFN = fullfile(dirname, sprintf('log_%.2d.txt', numel(dlogs) + 1));

%% initialize algorithm
MexIO('init',p);      %SC Set the initial (default) parameters

msglog(logFN, ['runExperiment started at ', datestr(clock)]);

if ((p.frameShift-round(p.frameShift)~=0) || (p.frameShift>p.frameLen))
    uiwait(errordlg(['Frameshift = ' num2str(p.frameShift) ' is a bad value. Set nWin and frameLen appropriately. Frameshift must be an integer & Frameshift <= Framelen'],'!! Error !!'))
    return
else

    msglog(logFN, ' ');
    msglog(logFN, ' ');
    TransShiftMex(0);           %SC Gives input/output device info, and serves as an initialization.

    msglog(logFN, 'Settings :')
    msglog(logFN, sprintf('DMA Buffer    = %i samples',p.frameLen)); %SC Buffer length after downsampling
    msglog(logFN, sprintf('Samplerate    = %4.2f kHz',p.sr/1000));   %SC sampling rate after downsampling
    msglog(logFN, sprintf('Analysis win  = %4.2f msec',p.bufLen/p.sr*1000));
    msglog(logFN, sprintf('LPC  window   = %4.2f msec',p.anaLen/p.sr*1000));

    msglog(logFN, sprintf('Process delay = %4.2f msec',p.nDelay*p.frameLen/p.sr*1000));
    msglog(logFN, sprintf('Process/sec   = %4.2f', p.sr/p.frameShift));

end

TransShiftMex(2);

%% Load the multi-talker babble noise
[x,fs_mtb]=wavread('mtbabble48k.wav');
lenMTB=round(2.5*fs_mtb);

% gainMTB_fb=dBSPL2WaveAmp(subject.lvNoise,1000)/sqrt(2)/calcMaskNoiseRMS;
% gainMTB_fb3=dBSPL2WaveAmp(subject.lvNoise3,1000,subject.pcrKnob)/sqrt(2)/calcMaskNoiseRMS;
% x_mtb=cell(1,3);
% x_mtb{1}=x(1:lenMTB);               x_mtb{1}=x_mtb{1}-mean(x_mtb{1});
% x_mtb{2}=x(lenMTB+1:lenMTB*2);      x_mtb{2}=x_mtb{2}-mean(x_mtb{2});
% x_mtb{3}=x(lenMTB*2+1:lenMTB*3);    x_mtb{3}=x_mtb{3}-mean(x_mtb{3});
% TransShiftMex(3,'datapb',x_mtb{1});

%% expt
figIdDat = makeFigDataMon;
% if expt.subject.trigByScanner == 0
figUFBDat = makeFigUserFB;
guidata(figUFBDat.fid, figUFBDat);
% end

if expt.subject.trigByScanner == 1
    set(figUFBDat.fid, 'visible', 'off')
end


% wordList=expt.words;

allPhases=expt.allPhases;
recPhases=expt.recPhases;
% nWords=length(wordList);

% if expt.subject.trigByScanner == 0
    hgui = UIRecorder('figIdDat', figIdDat, 'figUFBDat', figUFBDat);
% else
%     hgui = UIRecorder('figIdDat', figIdDat);
% end

% if (expt.subject.designNum==2)
%     expt.script=addFaceInfo(expt.script,hgui.skin.dFaces);
%     expt.dFaces=hgui.skin.dFaces;
% end
hgui.logFN = logFN;

hgui.interfaceMode = subject.interfaceMode;
hgui.simDataDir = subject.simDataDir;
if ~isempty(hgui.simDataDir)
    if ~isdir(hgui.simDataDir)
        error('Input sim data directory does not exist: %s', subject.simDataDir);
    end
end

hgui.stcsData = load(sentMatFN);

hgui.colors = colors;

hgui.pcrKnob=subject.pcrKnob;
hgui.ITI=expt.subject.ITI;
hgui.trigByScanner=expt.subject.trigByScanner;
hgui.TA=expt.subject.TA;
hgui.dBRange=expt.subject.dBRange1;
hgui.trialLen=expt.subject.trialLen;

hgui.modePromptDur = expt.subject.modePromptDur;
hgui.type1Prompt = expt.subject.type1Prompt;
hgui.type2Prompt = expt.subject.type2Prompt;
hgui.type3Prompt = expt.subject.type3Prompt;
hgui.type4Prompt = expt.subject.type4Prompt;

% hgui.skin.faceOrder=randperm(length(hgui.skin.dFaces));
hgui.skin.facePnt=1;

hgui.trialTypes = expt.subject.trialTypes;

hgui.manualTrigPhases = expt.subject.manualTrigPhases;

hgui.meanSylDur=expt.subject.paceStim.meanSylDur;
hgui.minSylDur = expt.subject.minSylDur;
hgui.maxSylDur = expt.subject.maxSylDur;
hgui.sylDurRange_R = expt.subject.sylDurRange_R;

hgui.minVwlLevel = expt.subject.(['minVwlLevel_', lower(expt.subject.sex)]);
hgui.maxVwlLevel = expt.subject.(['maxVwlLevel_', lower(expt.subject.sex)]);

hgui.showRhythmHint = expt.subject.showRhythmHint;
hgui.nonInformativeFixationCross = expt.subject.nonInformativeFixationCross;

hgui.audRhythmAlways = expt.subject.audRhythmAlways;

hgui.showRhythmicityFB_phases = expt.subject.showRhythmicityFB_phases;
hgui.showRateFB_phases = expt.subject.showRateFB_phases;
hgui.showIntFB_phases = expt.subject.showIntFB_phases;
% hgui.showRhythmicityWarn_phases = expt.subject.showRhythmicityWarn_phases;
hgui.showRateWarn_phases = expt.subject.showRateWarn_phases;
hgui.showIntWarn_phases = expt.subject.showIntWarn_phases;

hgui.showRateWarn_run_R_sameType = expt.subject.showRateWarn_run_R_sameType;

hgui.rateErrRepeat_phases = expt.subject.rateErrRepeat_phases;
hgui.intErrRepeat_phases = expt.subject.intErrRepeat_phases;

hgui.colorAcceptRanges = expt.subject.colorAcceptRanges;

hgui.showRhythmicityFB_onlyRhythm = expt.subject.showRhythmicityFB_onlyRhythm;

hgui.toneDur=expt.subject.paceStim.toneDur;
hgui.toneFreq=expt.subject.paceStim.toneFreq;
hgui.toneAmp=expt.subject.paceStim.toneAmp;
hgui.toneRamp=expt.subject.paceStim.toneRamp;
hgui.TPaceStim=expt.subject.TPaceStim;
hgui.TVisStim=expt.subject.TVisStim;

hgui.vumeterMode=expt.subject.vumeterMode;

hgui.rmsTransTarg_spl=getSPLTarg(expt.subject.mouthMicDist);
load('micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));

msglog(logFN, ' ');
msglog(logFN, ['Mouth-microphone distance = ',num2str(expt.subject.mouthMicDist),' cm']);
msglog(logFN, ['hgui.rmsTransTarg_spl = ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
msglog(logFN, ' ');

hgui.vocaLen=round(300*p.sr/(p.frameLen*1000)); % 300 ms, 225 frames
hgui.lenRange=round(250*p.sr/(p.frameLen*1000));  % single-sided tolerance range: 0.4*250 = 100 ms
msglog(logFN, ['Vowel duration range: [',num2str(300-0.4*250),',',num2str(300+0.4*250),'] ms.']);

hgui.debug = DEBUG;
hgui.trigKey = expt.subject.trigKey;

hgui.FmtShiftStat0 = expt.subject.FmtShiftStat0;
hgui.FmtShiftStat1 = expt.subject.FmtShiftStat1;

guidata(hgui.UIrecorder, hgui);

if (isempty(findStringInCell(varargin,'twoScreens')))
% 	set(hgui.UIrecorder,...
% 		'position', [0    5.0000  250.6667   65.8750],...
% 		'toolbar','none');  %SC Set the position of the expt window, partially for the use of multiple monitors.
else
% 	if (expt.subject.trigByScanner==1)
		ms=get(0,'MonitorPosition');
        if size(ms, 1) < 2
            error('It appears that only one monitor is connected to the system.');
        end
        
		set(hgui.UIrecorder,'Position',[ms(2,1),ms(1,4)-ms(2,4),ms(2,3)-ms(2,1)+1,ms(2,4)+20],'toolbar','none','doublebuffer','on','renderer','painters');
        set(hgui.UIrecorder, 'NumberTitle', 'off', 'Name', 'UI Recorder');
		pos_win=get(hgui.UIrecorder,'Position');
		pos_strh=get(hgui.strh,'Position');
		pos_axes_pic=get(hgui.axes_pic,'Position');
		pos_rms_axes=get(hgui.rms_axes,'Position');
		pos_speed_axes=get(hgui.speed_axes,'Position');
		pos_rms_label=get(hgui.rms_label,'Position');
		pos_rms_too_soft=get(hgui.rms_too_soft,'Position');
		pos_rms_too_loud=get(hgui.rms_too_loud,'Position');
		pos_speed_label=get(hgui.speed_label,'Position');
		pos_speed_too_slow=get(hgui.speed_too_slow,'Position');
		pos_speed_too_fast=get(hgui.speed_too_fast,'Position');
		set(hgui.strh,'Position',[(pos_win(3)-pos_strh(3))/2+5,(pos_win(4)-pos_strh(4))/2-75,pos_strh(3),pos_strh(4)]);                        
		set(hgui.axes_pic,'Position',[(pos_win(3)-pos_axes_pic(3))/2,(pos_win(4)-pos_axes_pic(4))/2,pos_axes_pic(3),pos_axes_pic(4)]);
		set(hgui.rms_axes,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_axes(2),pos_rms_axes(3),pos_rms_axes(4)]);
		set(hgui.rms_label,'Position',[(pos_win(3)-pos_rms_label(3))/2,pos_rms_label(2),pos_rms_label(3),pos_rms_label(4)]);
		set(hgui.rms_too_soft,'Position',[(pos_win(3)-pos_rms_axes(3))/2,pos_rms_too_soft(2),pos_rms_too_soft(3),pos_rms_too_soft(4)]);
		set(hgui.rms_too_loud,'Position',[(pos_win(3)-pos_rms_axes(3))/2+pos_rms_axes(3)-pos_rms_too_loud(3),pos_rms_too_loud(2),pos_rms_too_loud(3),pos_rms_too_loud(4)]);
		set(hgui.speed_axes,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_axes(2),pos_speed_axes(3),pos_speed_axes(4)]);		
		set(hgui.speed_label,'Position',[(pos_win(3)-pos_speed_label(3))/2,pos_speed_label(2),pos_speed_label(3),pos_speed_label(4)]);
		set(hgui.speed_too_slow,'Position',[(pos_win(3)-pos_speed_axes(3))/2,pos_speed_too_slow(2),pos_speed_too_slow(3),pos_speed_too_slow(4)]);
		set(hgui.speed_too_fast,'Position',[(pos_win(3)-pos_speed_axes(3))/2+pos_speed_axes(3)-pos_speed_too_fast(3),pos_speed_too_fast(2),pos_speed_too_fast(3),pos_speed_too_fast(4)]);
        set(hgui.msgh, 'FontSize', 20);
        
%         if hgui.trigByScanner == 0
        fpos = get(figUFBDat.fid, 'Position');
        set(figUFBDat.fid, 'Position', ...
            [ms(2, 1) + (ms(2, 3) - ms(2, 1) - fpos(3)) / 2, ...
             ms(1, 4) - ms(2, 4) + (ms(2, 4) - fpos(4)), ...
             fpos(3), fpos(4)]);
%         end
% 	else
% 		set(hgui.UIrecorder,'Position',[-1400,180,1254,857],'toolbar','none');
% 	end
	
end

if (subject.showProgress)
	set(hgui.progress_axes,'visible','on');
	set(hgui.progress_imgh,'visible','on');
	progress_meter=0.5*ones(1,100,3);
	progress_mask=zeros(1,100,3);
	set(hgui.progress_imgh,'Cdata',progress_meter.*progress_mask);        
else
	set(hgui.progress_axes,'visible','off');
	set(hgui.progress_imgh,'visible','off');
end



%% 

rProgress=0;
startPhase=state.phase; %SC For the purpose of resumed experiments
startRep=state.rep;     %SC For the purpose of resumed experiments

% --- Update timingDat.trialCnt --- %
timingDat.trialCnt = 1;
aPhases = fields(expt.script);
for i1 = 1 : startPhase - 1
    t_phase = aPhases{i1};
    timingDat.trialCnt = timingDat.trialCnt + expt.script.(t_phase).nTrials;
end

trialTypes = {'N', 'R'};

% trialCnt = 1;
for n=startPhase:length(allPhases)
    state.phase=n;
    state.rep=1;
    thisphase=allPhases{1,n};
    
    subdirname = fullfile(dirname,thisphase);
    if ~isdir(subdirname)
        mkdir(subdirname);
    end
    
    % -- Option to skip phases -- %
    if ~isempty(fsic(expt.skippablePhases, thisphase))
        inputOK = 0;
        while ~inputOK
            optDo = input(sprintf('Perform the optional phase {%s} (yes | no): ', thisphase), 's');
            if ~(isequal(lower(optDo), 'yes') || isequal(lower(optDo), 'no'))
                inputOK = 0;
                fprintf(2, 'Please input only yes or no\n');
            else
                inputOK = 1;
            end                                
        end

        if isequal(lower(optDo), 'yes')
            msglog(logFN, sprintf('User chose to include optional phase {%s}', thisphase));
        else
            msglog(logFN, sprintf('User chose to skip optional phase {%s}', thisphase));
        end

        if isequal(lower(optDo), 'no')
            continue;
        end
    end
    
    if expt.subject.trigByScanner == 1
        % -- Set the visibility of the user visual feedback window -- %
        if ~isempty(fsic(expt.skippablePhases, thisphase))
            set(figUFBDat.fid, 'visible', 'on');
        else
            set(figUFBDat.fid, 'visible', 'off');
        end
    end
    
    hgui.phase=thisphase;
    
    subjectMessageDialog(thisphase, getMsgStr(expt.subject.trigByScanner, thisphase));
    
    % === Get perturbation parameters === %
    if expt.subject.trigByScanner == 0 && ...
            length(thisphase) > 3 && isequal(thisphase(1 : 3), 'run')               
        sent = strrep(expt.subject.pertSent, '_', ' ');
        
        inputDirs = {};
        for i2 = 1 : n - 1
            if ~isequal(allPhases{i2}, 'pract1')
                inputDirs{end + 1} = fullfile(dirname, allPhases{i2});
            end
        end
        % TODO: Incremental updates
        
        [ostMat, timeWarpCfg, rmsThresh, maxInterval_5_10] = ...
            get_ost_pcf_params_rhy(inputDirs, sent, 0, ...
                                   'warpOnsetTime', expt.subject.warpOnsetTime, ...
                                   'decelWarpRate', expt.subject.decelWarpRate, ...
                                   'accelWarpRate', expt.subject.accelWarpRate, ...
                                   'F1ShiftRatio', expt.subject.F1ShiftRatio, ...
                                   'FmtShiftStat0', expt.subject.FmtShiftStat0, ...
                                   'FmtShiftStat1', expt.subject.FmtShiftStat1);
        
        % == Write ostMat to .ost file == %
        msglog(logFN, sprintf('\n'));
        ost_fns = struct;
        for i1 = 1 : numel(trialTypes)
            tt = trialTypes{i1};
            
            ost_fns.(tt) = fullfile(subdirname, sprintf('%s.ost', tt));
            write_ost_to_file(ostMat.(tt), expt.subject.rmsSlopeWin, ...
                              maxInterval_5_10.(tt), ost_fns.(tt));
            check_file(ost_fns.(tt));
            
            msglog(logFN, ...
                   sprintf('Generated online state tracking file for trialType %s: %s', ...
                           tt, ost_fns.(tt)));
        end
        msglog(logFN, sprintf('\n'));
        
        % == Write timeWarpCfg to .pcf file == %
        pcf_fmt_fns = struct;   % For formant perturbation
        pcf_twarp_fns = struct; % For time warping
        for i1 = 1 : numel(trialTypes)
            tt = trialTypes{i1}; 
            
            pcf_fmt_fns.(tt) = fullfile(subdirname, sprintf('fmt_%s.pcf', tt));
            format_pcf(expt.subject.basePCF, pcf_fmt_fns.(tt), ...
                       1234, ...        % Use a very large, practically unreachable stat number
                       timeWarpCfg.timeWarp_tBegin.(tt), ...
                       timeWarpCfg.timeWarp_rate1.(tt), ...
                       timeWarpCfg.timeWarp_dur1.(tt), ...
                       timeWarpCfg.timeWarp_durHold.(tt), ...
                       timeWarpCfg.timeWarp_rate2.(tt), ...
                       expt.subject.F1ShiftRatio, expt.subject.FmtShiftStat0, expt.subject.FmtShiftStat1);
            check_file(pcf_fmt_fns.(tt));
            
            msglog(logFN, sprintf('Generated formant perturbation configuration file for trialType %s: %s', ...
                                  tt, pcf_fmt_fns.(tt)));
            
            pcf_twarp_fns.(tt) = fullfile(subdirname, sprintf('twarp_%s.pcf', tt));
            format_pcf(expt.subject.basePCF, pcf_twarp_fns.(tt), ...
                       timeWarpCfg.timeWarp_initStat.(tt), ...
                       timeWarpCfg.timeWarp_tBegin.(tt), ...
                       timeWarpCfg.timeWarp_rate1.(tt), ...
                       timeWarpCfg.timeWarp_dur1.(tt), ...
                       timeWarpCfg.timeWarp_durHold.(tt), ...
                       timeWarpCfg.timeWarp_rate2.(tt), ...
                       0, expt.subject.FmtShiftStat0, expt.subject.FmtShiftStat1);
            check_file(pcf_twarp_fns.(tt));
            
            msglog(logFN, sprintf('Generated time-warping perturbation configuration file for trialType %s: %s', ...
                                  tt, pcf_twarp_fns.(tt)));
        end
        
        pcf_noPert_fn = fullfile(subdirname, 'noPert.pcf');
        format_pcf(expt.subject.basePCF, pcf_noPert_fn, ...
                   timeWarpCfg.timeWarp_initStat.(tt), ...
                   timeWarpCfg.timeWarp_tBegin.(tt), ...
                   1.0, ...
                   timeWarpCfg.timeWarp_dur1.(tt), ...
                   timeWarpCfg.timeWarp_durHold.(tt), ...
                   1.0, ...
                   0, expt.subject.FmtShiftStat0, expt.subject.FmtShiftStat1);
        check_file(pcf_noPert_fn);
        
        msglog(logFN, sprintf('Generated no-perturbation configuration file for trialType: %s', ...
                              pcf_noPert_fn));
    end
    msglog(logFN, sprintf('\n'));
    
    % Adjust the number of reps
    msglog(logFN, ['--- Coming up: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
        '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
    
    if expt.subject.trigByScanner
        nRepsNew = input('(Enter to skip) nRepsNew = ', 's');
    else
        nRepsNew = '';
    end
    
    nRepsNew=str2num(nRepsNew);
    if (~isempty(nRepsNew) && ~ischar(nRepsNew) && nRepsNew~=expt.script.(thisphase).nReps)
        expt.script.(thisphase).nReps=nRepsNew;
        expt.script.(thisphase)=genPhaseScript(thisphase, expt.script.(thisphase).nReps,...
            expt.trialTypes,expt.trainWords,expt.testWords,expt.pseudoWords,...
            expt.trialOrderRandReps,expt.subject.designNum);
        msglog(logFN, ['Changed: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
            '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
        save(fullfile(dirname,'expt.mat'),'expt');
        msglog(logFN, ['Saved ', fullfile(dirname,'expt.mat')]);
    end
    % Adjust the number of reps
    
    nReps=expt.script.(thisphase).nReps;
    if ~isequal(thisphase,'stay')
        phaseTrialCnt=1;
    end

    expt.script.(thisphase).startTime=clock;

	set(hgui.rms_axes,'visible','off');
    set(hgui.rms_imgh,'visible','off');
    set(hgui.rms_label,'visible','off');
	set(hgui.rms_too_soft,'visible','off');
	set(hgui.rms_too_loud,'visible','off');	
    set(hgui.speed_axes,'visible','off');
    set(hgui.speed_imgh,'visible','off');
    set(hgui.speed_label,'visible','off');
	set(hgui.speed_too_slow,'visible','off');
    set(hgui.speed_too_fast,'visible','off');
    
    hgui.bSpeedRepeat=0;
    hgui.bRmsRepeat=0;
	
	if (subject.showPlayButton==0)
		set(hgui.play,'visible','off');
    end
    
    set(hgui.play, 'cdata', hgui.skin.play, 'userdata', 0);
    
    if isequal(thisphase,'pre')
%         case 'pre'
        p.bDetect=0;
        p.bShift = 0;   %SC No shift in the pre phase

        hgui.bRmsRepeat=0;
        hgui.bSpeedRepeat=0;			

% 			set(hgui.rms_axes,'visible','off');
% 		    set(hgui.rms_imgh,'visible','off');
% 		    set(hgui.rms_label,'visible','off');
% 			set(hgui.rms_too_soft,'visible','off');
% 			set(hgui.rms_too_loud,'visible','off');	
% 		    set(hgui.speed_axes,'visible','off');
% 		    set(hgui.speed_imgh,'visible','off');
% 		    set(hgui.speed_label,'visible','off');
% 			set(hgui.speed_too_slow,'visible','off');
% 		    set(hgui.speed_too_fast,'visible','off');

        if (~isempty(rmsPeaks))
            p.rmsMeanPeak=mean(rmsPeaks);
            p.rmsThresh=p.rmsMeanPeak/4;       %SC !! Adaptive RMS threshold setting. Always updating 
        end

        hgui.showTextCue=0;
        
    elseif isequal(thisphase(1:3),'run')
        % SC(2008/06/10) Manually determine the optimum tracking params
        % Warning: for consistency, don't change nDelay			
        set(hgui.msgh,'visible','on');
%             set(hgui.msgh_imgh,'CData',CDataMessage.ftparampicking,'visible','on');
        drawnow;

        set(hgui.msgh,'string',{'Please stand by...'},'visible','on');        
        
        if (hgui.debug==0)
%             [vowelF0Mean,vowelF0SD]=getVowelPitches(dirname);
%             disp(['Vowel meanF0 = ',num2str(vowelF0Mean),' Hz: stdF0 = ',num2str(vowelF0SD),' Hz']);
%             disp(['Recommended cepsWinWidth = ',num2str(round(p.sr/vowelF0Mean*0.54))]);
%             [vowelF1Mean,vowelF2Mean]=getVowelMeanF1F2(dirname);
%             disp(['Vowel meanF1 = ',num2str(vowelF1Mean),' Hz; meanF2 = ',num2str(vowelF2Mean),' Hz']);
%             optimFTParams=compareFormTrackParams(dirname);
%             p.nLPC=optimFTParams.nLPC;
%             p.nDelay=optimFTParams.nDelay;
%             p.bufLen=(2*p.nDelay-1)*(p.frameLen);
%             p.anaLen=p.frameShift+2*(p.nDelay-1)*p.frameLen;
%             p.avgLen=optimFTParams.avgLen;
%             p.bCepsLift=optimFTParams.bCepsLift;
%             p.cepsWinWidth=optimFTParams.cepsWinWidth;
%             p.fn1=optimFTParams.fn1;
%             p.fn2=optimFTParams.fn2;       
%             p.aFact=optimFTParams.aFact;
%             p.bFact=optimFTParams.bFact;
%             p.gFact=optimFTParams.gFact;
        end
        % ~SC(2008/06/10) Manually determine the optimum tracking

        set(hgui.msgh,'string',{''},'visible','on'); 
        p.rmsMeanPeak=mean(rmsPeaks);
        p.rmsThresh=p.rmsMeanPeak / 4;    %SC !! Adaptive RMS threshold setting. Always updating
% 			if (p.rmsThresh>0.015)
% 				p.rmsThresh=0.015;
% 				disp('********* Warning: rms too high! Limited at 0.015. *********');
% 			end
        p.bDetect=0;
        p.bShift=0;
        set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
        hgui.showTextCue=1;

    end
    drawnow
    

%     set(hgui.msgh,'string',getMsgStr(thisphase),'visible','on');
    % TODO %

    MexIO('init',p);  %SC Inject p to TransShiftMex

    for i0 = startRep:nReps    %SC Loop for the reps in the phase
        repString=['rep',num2str(i0)];
        state.rep=i0;
        state.params=p;
        save(fullfile(dirname,'state.mat'),'state');
        
        nTrials=length(expt.script.(thisphase).(repString).trialOrder);

        subsubdirname=fullfile(subdirname,repString);
        mkdir(subsubdirname);
		
		% --- Perturbation field ---
		p.pertF2=linspace(p.F2Min,p.F2Max,p.pertFieldN);
%         switch (thisphase)
%             case 'ramp'
% 				p.pertAmp=i0/(expt.script.ramp.nReps+1)*subject.shiftRatio*ones(1,p.pertFieldN);
%             case 'stay'
%                 p.pertAmp=subject.shiftRatio*ones(1,p.pertFieldN);				
%             otherwise,
        p.pertAmp=zeros(1,p.pertFieldN);	
% 		end
% 		if isequal(subject.shiftDirection,'F1Up')
% 			p.pertPhi=0*ones(1,p.pertFieldN);
% 		elseif isequal(subject.shiftDirection,'F1Down')
% 			p.pertPhi=pi*ones(1,p.pertFieldN);
% 		end
		MexIO('init',p);
		% --- ~Perturbation field ---

        for k = 1 : nTrials
            thisTrial = expt.script.(thisphase).(repString).trialOrder(k); % 0: silent; 1: no noise; 2: noise only;
            
            if expt.subject.trigByScanner == 0
                thisPert = expt.script.(thisphase).(repString).pertType(k); % 0: noPert; 1: F1; 2: Decel; 4: other
            else
                thisPert = -1; % Not applicable
            end
            
            thisWord=expt.script.(thisphase).(repString).word{k};     %SC Retrieve the word from the randomly shuffled list
            nSyls=expt.script.(thisphase).(repString).nSyls(k);

			hgui.trialType = thisTrial;
            hgui.thisPert = thisPert;
			hgui.word = thisWord;            
            hgui.subject = expt.subject;
            hgui.nSyls = nSyls;
            hgui.saveDataFN = fullfile(subsubdirname, ['trial-', num2str(k), '-', num2str(thisTrial)]);
            hgui.asrDir = fullfile(subsubdirname, ['trial-', num2str(k), '-', num2str(thisTrial), '_asr']);
            
            % -- Get the next trial type -- %
            if k < nTrials
                hgui.nextTrialType = expt.script.(thisphase).(repString).trialOrder(k + 1);
            elseif i0 < nReps
                hgui.nextTrialType = expt.script.(thisphase).(['rep', num2str(i0 + 1)]).trialOrder(1);
            else
                hgui.nextTrialType = NaN;
            end
            
            if expt.subject.bAdaptRate
                hgui.meanSylDur=adaptMeanSylDur;
            end

%             if (hgui.trialType==2 | hgui.trialType==3)	% Speech with masking noise or passively listening to masking noise
%                 TransShiftMex(3,'datapb',gainMTB_fb*x_mtb{3-mod(k,3)},0);
% 			end
               
			msglog(logFN, '');
            if (ischar(thisWord))
    			msglog(logFN, ['[', datestr(clock), ']: ', thisphase,' - ',repString, ...
                      ', k = ',num2str(k),': trialType = ',num2str(hgui.trialType), ...
                      '; pertType = ', num2str(hgui.thisPert), ...
                      ' - ',thisWord, ' (trialLen = ',num2str(hgui.trialLen), ' s)']);
            else
                msglog(logFN, ['[', datestr(clock), ']: ', thisphase,' - ',repString, ...
                     ', k = ',num2str(k),': trialType = ',num2str(hgui.trialType), ...
                     '; pertType = ', num2str(hgui.thisPert), ...
                     ' - Pseudoword-', num2str(thisWord), ' (trialLen = ',num2str(hgui.trialLen), ' s)']);
            end
            
            % Count down    
            if ~(isequal(thisphase,'start') || isequal(thisphase,'ramp') || isequal(thisphase,'stay') || isequal(thisphase,'end'))
                msglog(logFN, ['Left: ',num2str(expt.script.(thisphase).nTrials-phaseTrialCnt+1),'/',num2str(expt.script.(thisphase).nTrials)]);
            else 
                if ~(isequal(thisphase,'ramp') || isequal(thisphase,'stay'))
                    msglog(logFN, ['Left: ',num2str(expt.script.(thisphase).nTrials-phaseTrialCnt+1),'/',num2str(expt.script.(thisphase).nTrials),...
                        ', ',num2str((expt.script.(thisphase).nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
                else
                    msglog(logFN, ['Left: ',num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1),'/',...
                        num2str(expt.script.ramp.nTrials+expt.script.stay.nTrials),...
                        ', ',num2str((expt.script.ramp.nTrials+expt.script.stay.nTrials-phaseTrialCnt+1)*hgui.ITI),' sec']);
                end
            end
            % ~Count down
            
            bPrompt = 0;
            if (hgui.trialType >= 2)   %SC The distinction between train and test words                             
                TransShiftMex(3, 'bdetect', 0, bPrompt);
                TransShiftMex(3, 'bshift', 0, bPrompt);
			else
                TransShiftMex(3, 'bdetect', p.bDetect, bPrompt);
                TransShiftMex(3, 'bshift', p.bShift, bPrompt);                
            end
            
            if (thisTrial==5)
				hgui.skin.facePnt=expt.script.(thisphase).(repString).face(k);
            end
            
            if hgui.interfaceMode == 2 % winAudio + sim 
                simDir = fullfile(hgui.simDataDir, thisphase, sprintf('rep%d', i0));
                if ~isdir(simDir)
                    error('Cannot find subdirectory in sim data directory: %s', simDir);
                end
                
                dsimfns = dir(fullfile(simDir, sprintf('trial-*-%d.mat', hgui.trialType)));
                if isempty(dsimfns)
                    error('Cannot find sim data file for phase %s, rep %d, trial %d', ...
                          thisphase, i0, k);
                end
                hgui.simDataFN = fullfile(simDir, dsimfns(1).name);
            end
            
            % == Perturbation related configurations == %
            if expt.subject.trigByScanner == 0 && ...
               length(thisphase) > 3 && isequal(thisphase(1 : 3), 'run')
                if thisTrial == 1
                    TransShiftMex(3, 'rmsthr', rmsThresh.N);
                elseif thisTrial == 2                    
                    TransShiftMex(3, 'rmsthr', rmsThresh.R);
                end
                
                tt = trialTypes{thisTrial};
                if thisPert == 1  % F1 shift
                    TransShiftMex(3, 'bbypassfmt', 0, 0);
                    TransShiftMex(3, 'bshift', 1, 0);
                    TransShiftMex(3, 'bpitchshift', 0, 0);
                    
                    TransShiftMex(8, ost_fns.(tt), 0);
                    TransShiftMex(9, pcf_fmt_fns.(tt), 0);
                    
                    % --- Take care of the feedback intensity mismatch issue --- %
                    TransShiftMex(3, 'scale', p.dScale * 1.4);
                    
                elseif thisPert == 2 % Time warping
                    TransShiftMex(3, 'bbypassfmt', 1, 0);
                    TransShiftMex(3, 'bshift', 0, 0);
                    TransShiftMex(3, 'bpitchshift', 1, 0);
                    
                    TransShiftMex(8, ost_fns.(tt), 0);
                    TransShiftMex(9, pcf_twarp_fns.(tt), 0);
                    
                    TransShiftMex(3, 'scale', p.dScale);
                    
                else % noPert
                    TransShiftMex(3, 'bbypassfmt', 1, 0);
                    TransShiftMex(3, 'bshift', 0, 0);
                    TransShiftMex(3, 'bpitchshift', 1, 0);
                    
                    TransShiftMex(8, ost_fns.(tt), 0);
                    TransShiftMex(9, pcf_noPert_fn, 0);
                    
                    TransShiftMex(3, 'scale', p.dScale);
                    
                end
            else
                TransShiftMex(3, 'bbypassfmt', 0, 0);
                TransShiftMex(3, 'bshift', 0, 0);
                TransShiftMex(3, 'bpitchshift', 0, 0);
                
                TransShiftMex(3, 'scale', p.dScale * 1.4);
            end
            TransShiftMex(3, 'nfb', 1, 0);
            % == ~Perturbation related configurations == %
            
            
            % == Adjust trial length (For behavioral sessions only) == %
            if expt.subject.trigByScanner == 0
                suffix = thisphase;
                if length(suffix) > 3 && isequal(suffix(1 : 3), 'run')
                    suffix = suffix(1 : 3);
                end
                
                if length(thisphase) > 5 && isequal(thisphase(1 : 5), 'inter')
                    hgui.trialLen = expt.subject.trialLen_pre;
                else
                    hgui.trialLen = expt.subject.(['trialLen_', suffix]);
                end
                
                TransShiftMex(3, 'triallen', hgui.trialLen, 0);                
            end
            
            b_fMRI_practInter = (expt.subject.trigByScanner == 1 && (thisTrial == 1 || thisTrial == 2) ...
                                 && length(thisphase) > 5 && (isequal(thisphase(1 : 5), 'pract') || isequal(thisphase(1 : 5), 'inter')));
            % == Loop trial until trial is performed correctly or no
            % repetition requirement is in place == %
            bRepeat = 1;
            repeatCnt = 1;
            while bRepeat
                MexIO('reset');

                UIRecorder('singleTrial', hgui.play, 1, hgui);

                data = get(hgui.UIrecorder,'UserData');           %SC Retrieve the data
                hgui1 = guidata(hgui.UIrecorder);
                
                % == Save timing data == %
                tDatIdx = get_trial_index(expt.script, thisphase, i0, k);
                if expt.subject.trigByScanner == 0 || b_fMRI_practInter % Behavioral sessions or fMRI practice / inter trials 
                    bToSaveTimingDat = timingDat.trialType(tDatIdx) <= 2;
                else % fMRI sessions
                    tDatIdx = tDatIdx - 1;
                    if tDatIdx <= 0
                        bToSaveTimingDat = 0;
                    else
                        bToSaveTimingDat = (timingDat.trialType(tDatIdx) <= 2);
                    end
                end
                
                if bToSaveTimingDat
                    if isfield(hgui1, 't_mean_ivi')
                        msglog(logFN, sprintf('Mean IVI = %f s', hgui1.t_mean_ivi))
                        
                        if expt.subject.trigByScanner == 1
                            bSaveDat = hgui1.bASRGood;
                        else
                            bSaveDat = data.bASRGood;
                        end
                        
                        if bSaveDat
                            timingDat.mean_ivi(tDatIdx) = hgui1.t_mean_ivi;
                            timingDat.cv_ivi(tDatIdx) = hgui1.t_cv_ivis;
                        else
                            fprintf(2, 'WARNING: ASR result from the last speech trials appears to contain errors. Not using the ASR results.\n');
                        end

                        if timingDat.trialType(tDatIdx) == 1
                            if bSaveDat
                                trackMeanSylDurs(end + 1) = hgui1.t_mean_ivi;
                            end
                            if length(trackMeanSylDurs) > 5 && ...
                               ~isnan(nanmean(trackMeanSylDurs(end - 5 + 1 : end))) && ...
                               (~isequal(thisphase, 'pract1') && ~isequal(thisphase, 'pract2') && ~isequal(thisphase, 'pre'))
                                adaptMeanSylDur = nanmean(trackMeanSylDurs(end - 5 + 1 : end));
                            end
                            recMeanSylDurs.nonRhythm(end+1) =  hgui1.t_mean_ivi;
                            recMeanPeakRMS.nonRhythm(end+1) = hgui1.t_mean_vwl_lv;
                            recCV_IVI.nonRhythm(end + 1) = hgui1.t_cv_ivis;
                        else
                            recMeanSylDurs.rhythm(end+1) = hgui1.t_mean_ivi;
                            recMeanPeakRMS.rhythm(end+1)  = hgui1.t_mean_vwl_lv;
                            recCV_IVI.rhythm(end + 1) = hgui1.t_cv_ivis;
                        end

                        save(fullfile(dirname, 'timingDat.mat'), ...
                             'timingDat', 'trackMeanSylDurs', 'adaptMeanSylDur', ...
                             'recMeanSylDurs', 'recMeanPeakRMS', 'recCV_IVI');
                    end
                end

                data.timeStamp=clock;

%                 if (thisTrial==1 || thisTrial==2)
                    set(0,'CurrentFigure', figIdDat(1));
                    set(gcf,'CurrentAxes',figIdDat(5));
                    cla;
                    plot(recMeanSylDurs.nonRhythm, '.-', 'Color', colors.nonRhythm);
                    hold on;
                    plot(recMeanSylDurs.rhythm, '.-', 'Color', colors.rhythm); 
%                     set(gca,'XLim',[0, 100],'YLim',[0,1]);
                    set(gca,'YLim',[0, 0.75]);
                    legend({'nonRhythm','rhythm'});
    %                 xlabel('Block #');
                    ylabel('Mean syllable duration (sec)');

                    set(gcf,'CurrentAxes',figIdDat(6));
                    cla;
%                     plot(20*log10(recMeanPeakRMS.nonRhythm), '.-', 'Color', colors.nonRhythm);
                    plot(recMeanPeakRMS.nonRhythm, '.-', 'Color', colors.nonRhythm);
                    hold on;
%                     plot(20*log10(recMeanPeakRMS.rhythm), '.-', 'Color', colors.rhythm);
                    plot(recMeanPeakRMS.rhythm, '.-', 'Color', colors.rhythm);
%                     set(gca,'XLim',[0,100]);
                    xlabel('Block #');
                    ylabel('Log Mean peak RMS (a.u.)');

                    set(gcf, 'CurrentAxes', figIdDat(8));
                    cla;
                    plot(recCV_IVI.nonRhythm, '.-', 'Color', colors.nonRhythm); hold on;
                    plot(recCV_IVI.rhythm, '.-', 'Color', colors.rhythm);
%                     set(gca,'XLim',[0,100]);
%                     legend({'nonRhythm','rhythm'});
                    ylabel('CV of IVIs');
%                 end

    %             if (thisTrial==1)
    %                 if ~isempty(data.rms)
    %                     switch (thisphase)  %SC Record the RMS peaks in the bout
    %                         case 'pre'
    %                             rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];
    %                         case 'pract1',
    %                             rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
    %                         case 'pract2',
    %                             rmsPeaks=[rmsPeaks ; max(data.rms(:,1))];                    
    %                         otherwise,
    %                     end
    %                 end
    %             end

    %             if (isequal(thisphase,'pract1'))
    %                 if (thisTrial==1 || thisTrial==2)
    %                     if (isfield(data,'vowelLevel') && ~isempty(data.vowelLevel) && ~isnan(data.vowelLevel) && ~isinf(data.vowelLevel))
    %                         subjProdLevel=[subjProdLevel,data.vowelLevel];
    %                     end
    %                 end
    %             end
                
                if expt.subject.trigByScanner && ~b_fMRI_practInter
                    bRepeat = 0;   % -- Override -- %
                    dataFN = fullfile(subsubdirname, ...
                                      ['trial-', num2str(k), '-', num2str(thisTrial), '.mat']);
                else 
                    bRepeat = 0;
                    if isequal(expt.subject.intErrRepeat_phases, 'all') || ...
                        ~isempty(fsic(strsplit(expt.subject.intErrRepeat_phases, ','), thisphase))
                        bRepeat = bRepeat | (data.bRmsGood == 0);
                    end

                    if isequal(expt.subject.rateErrRepeat_phases, 'all') || ...
                        ~isempty(fsic(strsplit(expt.subject.rateErrRepeat_phases, ','), thisphase))
                        bRepeat = bRepeat | (data.bSpeedGood == 0);
                    end

                    if bRepeat
                        dataFN = fullfile(subsubdirname, ...
                                          ['trial-', num2str(k), '-', num2str(thisTrial), '_bad', num2str(repeatCnt), '.mat']);

                        % -- Rename the ASR directory -- %
                        asrDir1 = strrep(hgui.asrDir, '_asr', sprintf('_bad%d_asr', repeatCnt));
                        if isdir(hgui.asrDir)
                            movefile(hgui.asrDir, asrDir1);
                            check_dir(asrDir1);
                            msglog(logFN, sprintf('Move directory: %s --> %s', hgui.asrDir, asrDir1));
                        end

                        repeatCnt = repeatCnt + 1;
                    else
                        dataFN = fullfile(subsubdirname, ...
                                          ['trial-', num2str(k), '-', num2str(thisTrial), '.mat']);
                    end

                end
                
                save(dataFN, 'data');
                check_file(dataFN);
                msglog(logFN, ['Saved ', dataFN,'.']);
                msglog(logFN, ' ');
            end
            
%             if thisTrial <= 2
                timingDat.trialCnt = timingDat.trialCnt + 1;
%             end
            phaseTrialCnt=phaseTrialCnt+1;

            % Calculate and show progress
            if (subject.showProgress)
                rProgress=calcExpProgress(expt,thisphase,i0,k,rProgress);
				if (~isnan(rProgress))
	                progress_mask=zeros(size(progress_meter));
		            progress_mask(:,1:round(rProgress*100),:)=1;            
			        set(hgui.progress_imgh,'Cdata',progress_meter.*progress_mask);
				end
            end
            
        end
    end
        
    startRep=1;
end
set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
set(hgui.msgh,'string',...
	{'Congratulations!';...
	'You have finished this experiment'},'visible','on');
% set(hgui.msgh_imgh,'CData',CDataMessage.finish,'visible','on');
pause(3);
close(hgui.UIrecorder)
pause(2);
% saveExperiment(dirname);

save(fullfile(dirname,'expt.mat'),'expt');
save(fullfile(dirname,'state.mat'),'state');

return

