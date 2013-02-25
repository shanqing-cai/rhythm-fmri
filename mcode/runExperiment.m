function varargout = runExperiment(configFN, varargin)
%% CONFIG
DEBUG=0;

colors.rhythm = [0, 0, 1];
colors.nonRhythm = [0, 0.4, 0];

%% ---- Modify -----
subject = read_subject_config(configFN);

if isequal(getHostName, 'smcg_w510')
    subject.dataDir        	= 'E:\DATA\RHYTHM-FMRI\';
else
    subject.dataDir         = 'D:\CS_2004\PROJECTS\RHYTHM-FMRI\';
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
    messg={'The specified directory already contains a previously recorded experiment'
        ''
        'Continue experiment, overwrite  or cancel ?'};
    button1 = questdlg(messg,'DIRECTORY NOT EMPTY','Continue','Overwrite','Cancel','Continue');
    switch button1
        case 'Overwrite'
            button2 = questdlg({'Are you sure you want to overwrite experiment'} ,'OVERWRITE EXPERIMENT ?');
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

% Make a copy of the configuration file
dconfig = dir(fullfile(dirname, 'config*.txt'));
configBackupFN = fullfile(dirname, sprintf('config_%.2d.txt', length(dconfig) + 1));
copyfile(configFN, configBackupFN);

if bNew % set up new experiment
    mkdir(dirname)
    
	expt.subject=subject;
    expt.allPhases={'pre','run1','run2','run3','run4','run5','run6'};
    expt.recPhases={'pre','run1','run2','run3','run4','run5','run6'}; %SC The pahses during which the data are recorded
    
    expt.trialTypes=[1,2,3,4];  % 1: non-rhythmic speech, 2: rhythmic speech, 3: non-rhythmic baseline, 4: rhythmic baseline. 
    expt.trialOrderRandReps=1;	%How many reps are randomized together
    
    expt.script.pre.nReps=2;    %SC Numbers of repetitions in the stages   % !!1!!	
    expt.script.run1.nReps=8;  %SC Default 10   %SC-Mod(09/26/2007)       % !!8!!
    expt.script.run2.nReps=8;   %SC Default 15   %SC-Mod(09/26/2007)       % !!2!!
    expt.script.run3.nReps=8;   %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!
    expt.script.run4.nReps=8;    %SC Default 20   %SC-Mod(09/26/2007)       % !!8!!
    expt.script.run5.nReps=8;
    expt.script.run6.nReps=8;

	expt.trialTypeDesc=cell(1,5);
	expt.trialTypeDesc{1}='Non-rhythmic speech';
	expt.trialTypeDesc{2}='Rhythmic speech';
	expt.trialTypeDesc{3}='Non-rhythmic baseline';
	expt.trialTypeDesc{4}='Rhythmic baseline';
    
    nSents=(expt.script.pre.nReps+expt.script.run1.nReps+expt.script.run2.nReps+expt.script.run3.nReps+...
           expt.script.run4.nReps+expt.script.run5.nReps+expt.script.run6.nReps)*2;
    [expt.stimSents_all,expt.nSyls_all]=getRandSentences(nSents);
    sentCnt=1;

	nPhases=length(expt.allPhases);
    sentCnt=1;
    for i1=1:nPhases
        t_nSents=expt.script.(expt.allPhases{i1}).nReps*2;
        expt.stimSents.(expt.allPhases{i1})=expt.stimSents_all(sentCnt:sentCnt+t_nSents-1);
        expt.stimSents_nSyls.(expt.allPhases{i1})=expt.nSyls_all(sentCnt:sentCnt+t_nSents-1);
        sentCnt=sentCnt+t_nSents;
    end
    
    expt.script.pre  = genPhaseScript('pre',  expt.script.pre.nReps,  expt.trialTypes, expt.stimSents.pre,  expt.stimSents_nSyls.pre,  expt.trialOrderRandReps);
    expt.script.run1 = genPhaseScript('run1', expt.script.run1.nReps, expt.trialTypes, expt.stimSents.run1, expt.stimSents_nSyls.run1, expt.trialOrderRandReps);
    expt.script.run2 = genPhaseScript('run2', expt.script.run2.nReps, expt.trialTypes, expt.stimSents.run2, expt.stimSents_nSyls.run2, expt.trialOrderRandReps);
    expt.script.run3 = genPhaseScript('run3', expt.script.run3.nReps, expt.trialTypes, expt.stimSents.run3, expt.stimSents_nSyls.run3, expt.trialOrderRandReps);
    expt.script.run4 = genPhaseScript('run4', expt.script.run4.nReps, expt.trialTypes, expt.stimSents.run4, expt.stimSents_nSyls.run4, expt.trialOrderRandReps);
    expt.script.run5 = genPhaseScript('run5', expt.script.run5.nReps, expt.trialTypes, expt.stimSents.run5, expt.stimSents_nSyls.run5, expt.trialOrderRandReps);
    expt.script.run6 = genPhaseScript('run6', expt.script.run6.nReps, expt.trialTypes, expt.stimSents.run6, expt.stimSents_nSyls.run6, expt.trialOrderRandReps);
    
    p=getTSMDefaultParams(subject.sex,'closedLoopGain',expt.subject.closedLoopGain,...
        'trialLen',expt.subject.trialLen,...
        'mouthMicDist',expt.subject.mouthMicDist);
    
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

    msglog(logFN, '\nSettings :')
    msglog(logFN, sprintf('DMA Buffer    = %i samples',p.frameLen)); %SC Buffer length after downsampling
    msglog(logFN, sprintf('Samplerate    = %4.2f kHz',p.sr/1000));   %SC sampling rate after downsampling
    msglog(logFN, sprintf('Analysis win  = %4.2f msec',p.bufLen/p.sr*1000));
    msglog(logFN, sprintf('LPC  window   = %4.2f msec',p.anaLen/p.sr*1000));

    msglog(logFN, sprintf('Process delay = %4.2f msec',p.nDelay*p.frameLen/p.sr*1000));
    msglog(logFN, sprintf('Process/sec   = %4.2f', p.sr/p.frameShift));

end

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
      
figIdDat=makeFigDataMon;

% wordList=expt.words;

allPhases=expt.allPhases;
recPhases=expt.recPhases;
% nWords=length(wordList);

hgui=UIRecorder('figIdDat', figIdDat);

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

hgui.pcrKnob=subject.pcrKnob;
hgui.ITI=expt.subject.ITI;
hgui.trigByScanner=expt.subject.trigByScanner;
hgui.TA=expt.subject.TA;
hgui.dBRange=expt.subject.dBRange1;
hgui.trialLen=expt.subject.trialLen;
% hgui.skin.faceOrder=randperm(length(hgui.skin.dFaces));
hgui.skin.facePnt=1;

hgui.meanSylDur=expt.subject.paceStim.meanSylDur;
hgui.minSylDur = expt.subject.minSylDur;
hgui.maxSylDur = expt.subject.maxSylDur;

hgui.toneDur=expt.subject.paceStim.toneDur;
hgui.toneFreq=expt.subject.paceStim.toneFreq;
hgui.toneAmp=expt.subject.paceStim.toneAmp;
hgui.toneRamp=expt.subject.paceStim.toneRamp;
hgui.TPaceStim=expt.subject.TPaceStim;
hgui.TVisStim=expt.subject.TVisStim;

hgui.vumeterMode=expt.subject.vumeterMode;

hgui.rmsTransTarg_spl=getSPLTarg(expt.subject.mouthMicDist);
load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
hgui.rmsTransTarg=micRMS_100dBA / (10^((100-hgui.rmsTransTarg_spl)/20));

msglog(logFN, ' ');
msglog(logFN, ['Mouth-microphone distance = ',num2str(expt.subject.mouthMicDist),' cm']);
msglog(logFN, ['hgui.rmsTransTarg_spl = ',num2str(hgui.rmsTransTarg_spl),' dBA SPL']);
msglog(logFN, ' ');

hgui.vocaLen=round(300*p.sr/(p.frameLen*1000)); % 300 ms, 225 frames
hgui.lenRange=round(250*p.sr/(p.frameLen*1000));  % single-sided tolerance range: 0.4*250 = 100 ms
msglog(logFN, ['Vowel duration range: [',num2str(300-0.4*250),',',num2str(300+0.4*250),'] ms.']);

hgui.debug=DEBUG;
% hgui.trigKey='equal';

if (isempty(findStringInCell(varargin,'twoScreens')))
% 	set(hgui.UIrecorder,...
% 		'position', [0    5.0000  250.6667   65.8750],...
% 		'toolbar','none');  %SC Set the position of the expt window, partially for the use of multiple monitors.
else
% 	if (expt.subject.trigByScanner==1)
		ms=get(0,'MonitorPosition');
		set(hgui.UIrecorder,'Position',[ms(2,1),ms(1,4)-ms(2,4),ms(2,3)-ms(2,1)+1,ms(2,4)+20],'toolbar','none','doublebuffer','on','renderer','painters');
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

% trialCnt = 1;
for n=startPhase:length(allPhases)
    state.phase=n;
    state.rep=1;
    thisphase=allPhases{1,n};
    subdirname=fullfile(dirname,thisphase);
    if ~isdir(subdirname)
        mkdir(subdirname);
    end
    
    hgui.phase=thisphase;
    
    % Adjust the number of reps
    msglog(logFN, ['--- Coming up: ',thisphase,'. nReps = ',num2str(expt.script.(thisphase).nReps),...
        '; nTrials = ',num2str(expt.script.(thisphase).nTrials),' ---']);
    nRepsNew=input('(Enter to skip) nRepsNew = ','s');
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
    
    if isequal(thisphase,'pre')
%         case 'pre'
        set(hgui.play,'cdata',hgui.skin.play,'userdata',0);
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
    

    set(hgui.msgh,'string',getMsgStr(thisphase),'visible','on');        

    MexIO('init',p);  %SC Inject p to TransShiftMex

    for i0=startRep:nReps    %SC Loop for the reps in the phase
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
            thisTrial=expt.script.(thisphase).(repString).trialOrder(k); % 0: silent; 1: no noise; 2: noise only; 			
            thisWord=expt.script.(thisphase).(repString).word{k};     %SC Retrieve the word from the randomly shuffled list
            nSyls=expt.script.(thisphase).(repString).nSyls(k);

			hgui.trialType = thisTrial;
			hgui.word = thisWord;
            hgui.subject = expt.subject;
            hgui.nSyls = nSyls;
            hgui.saveDataFN = fullfile(subsubdirname, ['trial-', num2str(k), '-', num2str(thisTrial)]);
            hgui.asrDir = fullfile(subsubdirname, ['trial-', num2str(k), '-', num2str(thisTrial), '_asr']);
            
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
                      ' - ',thisWord]);
            else
                msglog(logFN, ['[', datestr(clock), ']: ', thisphase,' - ',repString, ...
                     ', k = ',num2str(k),': trialType = ',num2str(hgui.trialType), ...
                     ' - Pseudoword-',num2str(thisWord)]);
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
            if (hgui.trialType>=2)   %SC The distinction between train and test words                             
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
            
            UIRecorder('singleTrial', hgui.play, 1, hgui);

            data=get(hgui.UIrecorder,'UserData');           %SC Retrieve the data
            hgui1 = guidata(hgui.UIrecorder);
            
            if timingDat.trialCnt > 1 && (timingDat.trialType(timingDat.trialCnt - 1) <= 2)
                if isfield(hgui1, 't_mean_ivi')
                    msglog(logFN, sprintf('Mean IVI = %f s', hgui1.t_mean_ivi))

                    timingDat.mean_ivi(timingDat.trialCnt - 1) = hgui1.t_mean_ivi;
                    timingDat.cv_ivi(timingDat.trialCnt - 1) = hgui1.t_cv_ivis;

                    if timingDat.trialType(timingDat.trialCnt - 1) == 1
                        trackMeanSylDurs(end + 1) = hgui1.t_mean_ivi;
                        adaptMeanSylDur = nanmean(trackMeanSylDurs(end - 5 + 1 : end));
                        recMeanSylDurs.nonRhythm(end+1) =  hgui1.t_mean_ivi;
                        recMeanPeakRMS.nonRhythm(end+1) = data.meanPeakRMS;
                        recCV_IVI.nonRhythm(end + 1) = hgui1.t_cv_ivis;
                    else
                        recMeanSylDurs.rhythm(end+1) = hgui1.t_mean_ivi;
                        recMeanPeakRMS.rhythm(end+1)  = data.meanPeakRMS;
                        recCV_IVI.rhythm(end + 1) = hgui1.t_cv_ivis;
                    end

                    save(fullfile(dirname, 'timingDat.mat'), ...
                         'timingDat', 'trackMeanSylDurs', 'adaptMeanSylDur', ...
                         'recMeanSylDurs', 'recMeanPeakRMS', 'recCV_IVI');
                end
            end
            
            timingDat.trialCnt = timingDat.trialCnt + 1;
            
            data.timeStamp=clock;
%             data.subject=expt.subject;
%             data.params.name=thisWord;
%             data.params.trialType=thisTrial;
            
%             if (thisTrial==1)
%                 trackMeanSylDurs=[trackMeanSylDurs,data.meanSylDur];
%                 adaptMeanSylDur = nanmean(trackMeanSylDurs(end - 5 + 1 : end));
%                 recMeanSylDurs.nonRhythm(end+1)=data.meanSylDur;
%                 recMeanPeakRMS.nonRhythm(end+1)=data.meanPeakRMS;
%             elseif (thisTrial==2)
%                 recMeanSylDurs.rhythm(end+1)=data.meanSylDur;
%                 recMeanPeakRMS.rhythm(end+1)=data.meanPeakRMS;
%             end
            
            if (thisTrial==1 || thisTrial==2)
                set(0,'CurrentFigure', figIdDat(1));
                set(gcf,'CurrentAxes',figIdDat(5));
                cla;
                plot(recMeanSylDurs.nonRhythm, '.-', 'Color', colors.nonRhythm);
                hold on;
                plot(recMeanSylDurs.rhythm, '.-', 'Color', colors.rhythm); 
                set(gca,'XLim',[0,100],'YLim',[0,1]);
                legend({'nonRhythm','rhythm'});
%                 xlabel('Block #');
                ylabel('Mean syllable duration (sec)');
                
                set(gcf,'CurrentAxes',figIdDat(6));
                cla;
                plot(20*log10(recMeanPeakRMS.nonRhythm), '.-', 'Color', colors.nonRhythm);
                hold on;
                plot(20*log10(recMeanPeakRMS.rhythm), '.-', 'Color', colors.rhythm);
                set(gca,'XLim',[0,100]);
                xlabel('Block #');
                ylabel('Log Mean peak RMS (a.u.)');
                
                set(gcf, 'CurrentAxes', figIdDat(8));
                cla;
                plot(recCV_IVI.nonRhythm, '.-', 'Color', colors.nonRhythm); hold on;
                plot(recCV_IVI.rhythm, '.-', 'Color', colors.rhythm);
                set(gca,'XLim',[0,100]);
                legend({'nonRhythm','rhythm'});
                ylabel('CV of IVIs');
            end
            
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
			
            save(fullfile(subsubdirname, ['trial-', num2str(k), '-', num2str(thisTrial)]), 'data');
            msglog(logFN, ['Saved ',fullfile(subsubdirname,['trial-',num2str(k),'-',num2str(thisTrial)]),'.']);
            msglog(logFN, ' ');
            
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
	'You have finished the expt.'},'visible','on');
% set(hgui.msgh_imgh,'CData',CDataMessage.finish,'visible','on');
pause(3);
close(hgui.UIrecorder)
pause(2);
% saveExperiment(dirname);

save(fullfile(dirname,'expt.mat'),'expt');
save(fullfile(dirname,'state.mat'),'state');

return

%% 
function phaseScript=genPhaseScript(stage,nReps,trialTypes,trainWords,nSyls,randReps)
	phaseScript=struct();
	phaseScript.nReps=nReps;
    phaseScript.nTrials=0;

    symbols='@#$%^&*()_+=-<>/\|[]{}';
    nSymbols=length(symbols);

    twCnt=1;
    ptwCnt=1;
    for n=1:nReps
        bt=[trialTypes];
        trainWordsUsed=trainWords;
%         pseudoWordsUsed=pseudoWords(randperm(length(trainWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));            
       
        bt=bt(randperm(length(bt)));
        oneRep=struct;
        oneRep.trialOrder=[];
        oneRep.word=cell(1,0);
        oneRep.nSyls=[];
        cntTW=1;
        for m=1:length(bt)
            oneRep.trialOrder=[oneRep.trialOrder,bt(m)];
            if (bt(m)==1 || bt(m)==2)					                
                oneRep.word{length(oneRep.word)+1}=trainWordsUsed{twCnt};
                oneRep.nSyls(end+1)=nSyls(twCnt);
                twCnt=twCnt+1;
            elseif (bt(m)==3 || bt(m)==4)
                t_sent=trainWordsUsed{ptwCnt};                
                for i1=1:length(t_sent)
                    if ~isequal(t_sent(i1),' ')
                        t_sent(i1)=symbols(round(rand*(nSymbols-1))+1);
                    end
                end
                oneRep.word{length(oneRep.word)+1}=t_sent;
                oneRep.nSyls(end+1)=nSyls(ptwCnt);
%                 oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW);
%                 oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
%                 cntTW=cntTW+2;
                ptwCnt=ptwCnt+1;
            end
        end

        phaseScript.(['rep',num2str(n)])=oneRep;
        phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);

        if (isequal(stage(1:3),'run') && n==nReps)
            phaseScript.nTrials=phaseScript.nTrials+1;
            
            idx0=find(phaseScript.(['rep',num2str(n)]).trialOrder==3,1);            
            phaseScript.(['rep',num2str(n)]).trialOrder(end+1)=phaseScript.(['rep',num2str(n)]).trialOrder(idx0);
            phaseScript.(['rep',num2str(n)]).word{end+1}=phaseScript.(['rep',num2str(n)]).word{idx0};
            phaseScript.(['rep',num2str(n)]).nSyls(end+1)=phaseScript.(['rep',num2str(n)]).nSyls(idx0);
        end
    end

return

% function script1=addFaceInfo(script0,dFaces)
%     script1=script0;
%     
%     IDs.male=[];
%     IDs.female=[];
%     
%     for n=1:length(dFaces.d)
%         if isequal(dFaces.sex{n},'M')
%             if isempty(find(IDs.male==dFaces.subjID(n)))
%                 IDs.male=[IDs.male,dFaces.subjID(n)];
%             end
%         elseif isequal(dFaces.sex{n},'F')
%             if isempty(find(IDs.female==dFaces.subjID(n)))
%                 IDs.female=[IDs.female,dFaces.subjID(n)];
%             end
%         end
%     end    
%     IDs.male=sort(IDs.male);
%     IDs.female=sort(IDs.female);
%     
%     IDs.male=IDs.male(randperm(length(IDs.male)));
%     IDs.female=IDs.female(randperm(length(IDs.female)));    
%     
%     stages=fields(script1);
%     for n=1:length(stages);
%         stg=stages{n};
%         if ~isempty(find(script1.(stg).rep1.trialOrder==5))
%             nPersons=script1.(stg).nReps;
%             tPersons=[];
%             for k=1:nPersons
%                 if (mod(k,2)==0)    % female
%                     tPersons=[tPersons,IDs.female(1)];
%                     IDs.female=IDs.female(2:end);
%                 else    % male
%                     tPersons=[tPersons,IDs.male(1)];
%                     IDs.male=IDs.male(2:end);
%                 end
%             end
%             idxFaces=[];
%             for k=1:length(tPersons)
%                 idxFaces0=find(dFaces.subjID==tPersons(k));
%                 idxFaces0=idxFaces0(randperm(length(idxFaces0)));
%                 idxFaces=[idxFaces,idxFaces0(1:4)];
%             end
%             idxFaces=idxFaces(randperm(length(idxFaces)));
%             cnt=1;
%             for k=1:script1.(stg).nReps
%                 repString=['rep',num2str(k)];
%                 script1.(stg).(repString).face=zeros(size(script1.(stg).(repString).trialOrder));
%                 idx=find(script1.(stg).(repString).trialOrder==5)
%                 script1.(stg).(repString).face(idx)=idxFaces(cnt:cnt+length(idx)-1);
%                 cnt=cnt+length(idx);
%             end
%         end
%     end
% return