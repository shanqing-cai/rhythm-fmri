function varargout = UIRecorder(varargin)
% UIRECORDER M-file for UIrecorder.fig
%      UIRECORDER, by itself, creates a new UIRECORDER or raises the existing
%      singleton*.
%
%      H = UIRECORDER returns the handle to a new UIRECORDER or the handle to
%      the existing singleton*.
%
%      UIRECORDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIRECORDER.M with the given input arguments.
%
%      UIRECORDER('Property','Value',...) creates a new UIRECORDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UIrecorder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UIrecorder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UIrecorder

% Last Modified by GUIDE v2.5 22-Jun-2008 15:54:34

%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UIrecorder_OpeningFcn, ...
    'gui_OutputFcn',  @UIrecorder_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% SCai: data displaying function
if (~isempty(findStringInCell(varargin,'figIdDat'))) 
    figIdDat=varargin{findStringInCell(varargin,'figIdDat')+1};
end

%% --- Executes just before UIrecorder is made visible.
function UIrecorder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
if (length(varargin)>=2)
    figIdDat=varargin{2};
else
    figIdDat=[];
end

if (~isempty(findStringInCell(varargin,'trad')))
	isTrad=1;
else 
	isTrad=0;
end

handles.trigByScanner=0;
handles.TA=2.5;
handles.phase='';
handles.trigKey='add';	% To change
handles.trialLen=2.2;
handles.debug=0;
handles.vumeterMode=NaN;  % 1: 10 ticks; 2: 3 ticks;

handles.timeCreated=clock;

if isequal(getHostName, 'smcg_w510');
    handles.msgImgDir = 'E:\speechres\adapt\triphcode\uimg';
    handles.utterImgDir = 'E:\speechres\adapt\triphcode\utterimg';
else
    handles.msgImgDir = 'D:\speechres\adapt\triphcode\uimg';
    handles.utterImgDir = 'D:\speechres\adapt\triphcode\utterimg';
end

if (isTrad)
	handles.msgImgDir=fullfile(handles.msgImgDir,'trad');
	handles.utterImgDir=fullfile(handles.utterImgDir,'trad');
end

set(hObject,'visible','off');
set(handles.UIrecorder,'interruptible','on','busyaction','queue')
%--------------------------------------------------------------------------
%SC Construct the volume/speed indicator template
vumeter  = permute(jet(100),[1,3,2]);
color1=vumeter(20,:,:);
color2=vumeter(53,:,:);
color3=vumeter(90,:,:);
for n=1:100
    if n<=30
        vumeter(n,:,:)=color1;
    elseif n<=70
        vumeter(n,:,:)=color2;
    else
        vumeter(n,:,:)=color3;
    end
end
vumeter2=vumeter;
vumeter(10:10:90,:)=0;    %SC-Commented(12/11/2007)
vumeter0=nan(size(vumeter,2),size(vumeter,1),size(vumeter,3));
vumeter0(:,:,1)=transpose(vumeter(:,:,1));
vumeter0(:,:,2)=transpose(vumeter(:,:,2));
vumeter0(:,:,3)=transpose(vumeter(:,:,3));
vumeter=vumeter0;
% vubounds=[1,30,70,100];%SC The boundaries are at 29 and 69.

vumeter2([30,70],:)=0;
vumeter02=nan(size(vumeter2,2),size(vumeter2,1),size(vumeter2,3));
vumeter02(:,:,1)=transpose(vumeter2(:,:,1));
vumeter02(:,:,2)=transpose(vumeter2(:,:,2));
vumeter02(:,:,3)=transpose(vumeter2(:,:,3));
vumeter2=vumeter02;
%--------------------------------------------------------------------------
%SC Construct the progress indicator template
progressmeter=1*ones(1,100,3);
%SC ~Construct the progress indicator template
%--------------------------------------------------------------------------
if (handles.vumeterMode==1)
    handles.rms_imgh = image(vumeter,'parent',handles.rms_axes);
    handles.speed_imgh=image(vumeter,'parent',handles.speed_axes);
    set(handles.rms_imgh,'CData',zeros(size(vumeter)));
    set(handles.speed_imgh,'CData',zeros(size(vumeter)));
else
    handles.rms_imgh = image(vumeter,'parent',handles.rms_axes);
    handles.speed_imgh=image(vumeter,'parent',handles.speed_axes);
    set(handles.rms_imgh,'CData',zeros(size(vumeter2)));
    set(handles.speed_imgh,'CData',zeros(size(vumeter2))); 
end

if (~isempty(findStringInCell(varargin,'showVuMeter')))    
    set(handles.rms_imgh,'CData',vumeter0);
end



% if (~isempty(findStringInCell(varargin,'showVuMeter')))
%     set(handles.speed_imgh,'CData',vumeter0);
% end

% set(handles.phrase_axes,'Box','off');
% set(handles.axes_msgh,'Box','off');
set(handles.axes_pic,'Box','off');


handles.progress_imgh = image(progressmeter,'parent',handles.progress_axes);

set(handles.rms_label,'string','Volume');
set(handles.speed_label,'string','Speed');

handles.pcrKnob=NaN;
% handles.trialType=4;
% handles.word='Ready...';

% handles.bAuto=1;

handles.time1=[];
handles.time2=[];

% set(handles.auto_btn,'Value',get(handles.auto_btn,'Max'));
% skin=struct('pause', imread(fullfile(pwd,'graphics','skin-pause.jpg')),...
%     'play', imread(fullfile(pwd,'graphics','skin-play.jpg')),...
%     'good', imread(fullfile(pwd,'graphics','choice-yes.gif')),...
%     'bad', imread(fullfile(pwd,'graphics','choice-cancel.gif')),...
% 	'fixation',imread(fullfile(pwd,'graphics','fixation.bmp')),...
%     'fixation_paced',imread(fullfile(pwd,'graphics','fixation-paced.bmp')),...
%     'fixation_nonpaced',imread(fullfile(pwd,'graphics','fixation-nonpaced.bmp')),...
% 	'vumeter',  vumeter,...
%     'vumeter2',  vumeter2,...
% 	'dFaces',getDFaces('./graphics/faces/face*.bmp'),...
% 	'dPseudowords',dir('./graphics/pseudochars/pseudoword-*.bmp'),...
%     'dWords',dir('./graphics/words/word*.bmp'),...
% 	'faceOrder',[],...
% 	'facePnt',1);
skin=struct('pause', imread(fullfile(pwd,'graphics','skin-pause.jpg')),...
    'play', imread(fullfile(pwd,'graphics','skin-play.jpg')),...
    'good', imread(fullfile(pwd,'graphics','choice-yes.gif')),...
    'bad', imread(fullfile(pwd,'graphics','choice-cancel.gif')),...
	'fixation',imread(fullfile(pwd,'graphics','fixation.bmp')),...
    'fixation_paced',imread(fullfile(pwd,'graphics','fixation-paced.bmp')),...
    'fixation_nonpaced',imread(fullfile(pwd,'graphics','fixation-nonpaced.bmp')),...
	'vumeter',  vumeter,...
    'vumeter2',  vumeter2,...
	'dPseudowords',dir('./graphics/pseudochars/pseudoword-*.bmp'),...
    'dWords',dir('./graphics/words/word*.bmp'));
handles.skin=skin;

handles.pic_imgh=image(handles.skin.fixation,'parent',handles.axes_pic);
set(handles.pic_imgh,'visible','off');

% set(handles.prev,'cdata',skin.prev);
% set(handles.next,'cdata',skin.next);
set(handles.play,'cdata',skin.play);
set(handles.play,'UserData',0);
set(handles.play,'Value',get(handles.play,'Min'));
set(handles.rms_axes,'xtick',[],'ytick',[]);
axis(handles.rms_axes, 'xy');

set(handles.speed_axes,'xtick',[],'ytick',[])
axis(handles.speed_axes, 'xy');

set(handles.axes_pic,'xtick',[],'ytick',[],'box','off','visible','off');

set(handles.progress_axes,'xtick',[],'ytick',[]);
axis(handles.progress_axes, 'xy');

handles.figIdDat=figIdDat;

handles.dataOut=[];
handles.bRmsRepeat=0;
handles.bSpeedRepeat=0;
handles.vocaLen=NaN;    %SC-Mod(2008/01/05) Old value: 300 
handles.lenRange=NaN;   %SC(2008/01/05)

handles.ITI=6;			%SC(2009/02/05) Inter-trial interval

handles.showTextCue=0;  %SC(2008/01/06)

handles.dBRange=NaN;
handles.rmsTransTarg_spl=NaN;
% load calibMic;  % gets micGain: wav rms/ Pa rms (Pa^-1)
load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
handles.rmsTransTarg=micRMS_100dBA / (10^((100-handles.rmsTransTarg_spl)/20));

handles.nextMessage=imread(fullfile(handles.msgImgDir,'message_pre2.bmp'));

set(handles.UIrecorder,'keyPressFcn',@key_Callback);
set(handles.strh,'keyPressFcn',@key_Callback);
set(handles.play,'keyPressFcn',@key_Callback);
% set(handles.rec_slider,'keyPressFcn',@key_Callback);
set(handles.msgh,'keyPressFcn',@key_Callback);

set(handles.strh,'string','HELLO','visible','on');

set(hObject,'visible','on');

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes UIrecorder wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% --- Outputs from this function are returned to the command line.
function varargout = UIrecorder_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = handles;

%% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play

% set(handles.button_next,'visible','off');

% CDataMessageBlank=zeros(750,720,3);
% CDataMessageBlank(:,:,1)=64/255*ones(750,720);
% CDataMessageBlank(:,:,3)=64/255*ones(750,720);

if(get(handles.play,'userdata')==0) % currently in pause mode
    set(handles.play,'cdata',handles.skin.pause,'userdata',1); % now in play mode
    set(handles.msgh,'string','');
	handles.trialType=-1;
	handles.word='Ready...';	

    singleTrial(handles.play,[],handles)   %%SC
else % currently in play mode
    set(handles.play,'cdata',handles.skin.play,'userdata',0); % now in pause mode
    TransShiftMex(2) %%SC stop TransShiftMex
    set(handles.msgh,'string','Press play to continue...');
end

function key_Callback(src, evnt)
hgui=guidata(src);
timeNow=clock;
eTime=etime(timeNow,hgui.timeCreated);

if (isequal(evnt.Key,hgui.trigKey) || isequal(evnt.Key,'a'))
% 	set(hgui.UIrecorder,'UserData','go');
    disp(['--> Trigger at ',num2str(eTime),' sec <--']);
	uiresume(hgui.UIrecorder);    
else
% 	set(hgui.UIrecorder,'UserData','nogo');
end

return

% %% --- Executes on button press in prev.
% function prev_Callback(hObject, eventdata, handles)
% % hObject    handle to prev (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% singleTrial(handles.prev,[],handles);

%% --- Executes on button press in next.
% function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% --- Single trial callback function.
function singleTrial(hObject, eventdata, handles, varargin)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.next,'enable','off','visible','off')
% set(handles.prev,'enable','off','visible','off')
record(handles);

%% SUBFUNCTIONS
%% --------------------------------------------------------------------------
function [dataOut, bRmsGood, bSpeedGood, t_rms] = checkData(data, handles, paceStimParams)
% This function is called after stoprec is executed. It checks the data and
% displays the rms and transition length . If the rms and speed are in range
%, the data is stored in handles.dataOut.

t0=1/data.params.sr;
taxis=0:t0:t0*(length(data.signalIn)-1);

dBrange=handles.dBRange; %SC-Mod(2008/01/05). One side tolerance: 0.4*dBrange 
rmsval=50;
speedval=0;
dataOut=data;
dataOut.paceStimParams=paceStimParams;

% if (handles.bMon==1)
% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
% [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));
[sentDur,meanPeakRMS,threshRMS,idx1,idx2] = analyzeSentence(data);
if ~isnan(handles.nSyls)
    meanSylDur=sentDur/handles.nSyls;
else
    meanSylDur=sentDur/8.5;
end

dataOut.meanSylDur=meanSylDur;
dataOut.meanPeakRMS=meanPeakRMS;

if ~isempty(data)
% 	vocaLen=handles.vocaLen;
% 	lenRange=handles.lenRange;
% 	rmsTransTarg=handles.rmsTransTarg;

% 	t1=k1*data.params.frameLen;
% 	t2=k2*data.params.frameLen;

% 	vocaLenNow=k2-k1+1;   %SC-Mod(2008/04/06)

	%SC-Mod(2008/04/06): Look at the rms during the transition, instead of
	%   during the entire vocal part.
% 	if (isnan(t1) | isnan(t2) | isempty(t1) | isempty(t2) | t1>=t2)
% 		rmsTrans=0;
% 		rmsBGNoise=0;
%         if ~isempty(data.signalIn)
%             rmsBGNoise=calcAWeightedRMS(data.signalIn(1:round(0.2*data.params.sr)),data.params.sr);
%         end
% 	else
% 		rmsTrans=sqrt(mean(data.signalIn(t1:t2).^2));
% 		rmsBGNoise=sqrt(mean(data.signalIn(1:round(0.2*data.params.sr)).^2));
% 		rmsTrans=calcAWeightedRMS(data.signalIn(t1:t2),data.params.sr);
		rmsBGNoise=calcAWeightedRMS(data.signalIn(1:round(0.2*data.params.sr)),data.params.sr);
% 	end

% 	rmsval   = round(100/dBrange*max(0,min(dBrange,dBrange/2+10*log10(rmsTrans/rmsTransTarg))));   %SC-Mod(2007/12/29)
% 	speedval = round(100/lenRange*max(0,min(lenRange,lenRange/2+(vocaLen-vocaLenNow)/2)));   
end

bRmsGood=1;
bSpeedGood=1;
%--------------------------------------------------------------------------
%SC Set the volume/speed indicator
% if (handles.vumeterMode==1)
%     vumeter=handles.skin.vumeter;
% elseif (handles.vumeterMode==2);
%     vumeter=handles.skin.vumeter2;
% end
% vumeter0=vumeter*0.5;
% vubounds=handles.skin.vubounds; %SC(12/11/2007)


% mask=zeros(size(vumeter));
% mask0=zeros(size(vumeter));

% rmsval1=floor(rmsval/10)*10;
% if (rmsval1+10>size(vumeter,2))
% 	rmsval1=size(vumeter,2)-10;
% end
% if (rmsval1<0)
% 	rmsval1=0;
% end
% if (handles.trialType==3 || handles.trialType==4 || handles.trialType==5)
% 	rmsval1=40+rand*10;
% 	if (rmsval1>45) rmsval1=50;
% 	else rmsval1=40;
% 	end
% end

% if (handles.vumeterMode==1)
%     mask(:,rmsval1+1:rmsval1+10,:) = 1;   %SC-Commented(12/11/2007)
%     mask0=1-mask;
% elseif (handles.vumeterMode==2)
%     if (rmsval1<30)
%         mask(:,1:30,:) = 1;
%     elseif (rmsval1>=30 && rmsval1<70)
%         mask(:,31:70,:) = 1;
%     else
%         mask(:,70:100,:) = 1;                
%     end
%     mask0=1-mask;    
% end

% set(handles.rms_imgh,'Cdata',vumeter.*mask+vumeter0.*mask0);
% 
% mask=zeros(size(vumeter));
% mask0=zeros(size(vumeter));
% 
% speedval1=floor(speedval/10)*10;
% if (speedval1+10>size(vumeter,2))
% 	speedval1=size(vumeter,2)-10;
% end
% if (speedval1<0)
% 	speedval1=0;
% end
% if (handles.trialType==3 | handles.trialType==4 | handles.trialType==5)
% 	speedval1=40+rand*10;
% 	if (speedval1>45) speedval1=50;
% 	else speedval1=40;
% 	end	
% end

% if (handles.vumeterMode==1)
%     mask(:,speedval1+1:speedval1+10,:) = 1;         %SC(2008/01/05)
%     mask0=1-mask;
% elseif (handles.vumeterMode==2)
%     if (speedval1<30)
%         mask(:,1:30,:) = 1;
%     elseif (speedval1>=30 && speedval1<70)
%         mask(:,31:70,:) = 1;
%     else
%         mask(:,70:100,:) = 1;          
%     end
%     mask0=1-mask;
% end

% if (handles.trialType==1 | handles.trialType==2)
% set(handles.speed_imgh,'Cdata',vumeter.*mask+vumeter0.*mask0);
    
% end
%SC ~Set the volume/speed indicator
%--------------------------------------------------------------------------
% drawnow

% msg1='';    msg2='';
% instr1='';  instr2='';
% if (rmsval < 70 && rmsval > 30)
% 	bRmsGood=1;
% else
% 	bRmsGood=0;
% 	if (rmsval >= 70)   %SC(2008/01/05)
% 		msg1='Softer';
% 		instr2='Loud';
% 	else % Then rmsval <= 30
% 		msg1='Louder';
% 		instr2='Soft';
% 	end
% end

% if (speedval < 70 && speedval > 30) %SC-Mod(2008/01/05) Used to be speedval > 20
% 	bSpeedGood=1;
% else
% 	bSpeedGood=0;
% 	if (speedval >= 70) %SC (2008/01/05)
% 		msg2='Slower';
% 		instr1='Fast';
% 	else
% 		msg2='Faster';
% 		instr1='Slow';
% 	end
% end

% if (handles.trialType==3 | handles.trialType==4 | handles.trialType==5)
% 	bRmsGood=1;
% 	bSpeedGood=1;
% end

%SC(2008/01/05)
% if (~bRmsGood | ~bSpeedGood)
% 	if (~bRmsGood & ~bSpeedGood)
% %         msgc=['Please speak ',msg1,' and ',msg2,'.'];
% 		msgc=[msg1,' and ',lower(msg2),' please!'];
% 	elseif (~bRmsGood)
% %         msgc=['Please speak ',msg1,'.'];
% 		msgc=[msg1,' please!'];        
% 	elseif (~bSpeedGood)
% 		msgc=['Please speak ',msg2,'.'];
% 		msgc=[msg2,' please!'];
% 	end
% 
% 	if (handles.showTextCue)
% 		set(handles.msgh,'string',{'';msgc});
% % 		pause(1);
% 	end
% end

% SCai: update the data monitor window
set(0,'CurrentFigure',handles.figIdDat(1));
set(gcf,'CurrentAxes',handles.figIdDat(2));
cla;
plot(taxis,data.signalIn);      hold on;
set(gca,'XLim',[taxis(1);taxis(end)]);
set(gca,'YLim',[-1,1]);
ylabel('Wave In');

set(gcf,'CurrentAxes',handles.figIdDat(3));
cla;
taxis=0:t0:t0*(length(data.signalOut)-1);
plot(taxis,data.signalOut*data.params.dScale);     hold on;
set(gca,'XLim',[taxis(1);taxis(end)]);
set(gca,'YLim',[-1,1]);
xlabel('Time (sec)');
ylabel('Wave Out');

% [i1,i2,f1,f2,iv1,iv2]=getFmtPlotBounds(data.fmts(:,1),data.fmts(:,2));
% [k1,k2]=detectVowel(data.fmts(:,1),data.fmts(:,2),iv1,iv2,'eh','rms',data.rms(:,1));
% if (~isnan(i1) && ~isnan(i2) && ~isempty(i1) && ~isempty(i2) && k2 >= k1)
% 	t1=k1*data.params.frameLen;
% 	t2=k2*data.params.frameLen;
% 	tv1=min(find(data.fmts(:,1)>0));
% 	tv2=max(find(data.fmts(:,1)>0));

% 	idx1=max([1,tv1-50]);
% 	idx2=min([tv2+50,length(data.signalIn)]);

%     wavInGain=0.13827;  % w/Pa
% 	p0=20e-6;           % Pa
% 	tRMSIn=sqrt(mean((data.signalIn(t1:t2)).^2));

% 	set(gcf,'CurrentAxes',handles.figIdDat(2));
% 	xs=get(gca,'XLim'); ys=get(gca,'YLim');
	
% 	if (~isnan(t1) && ~isnan(t2) && t1>0 && t2>0 && t2>t1)
% 		plot([taxis(t1),taxis(t1)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
% 		plot([taxis(t2),taxis(t2)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
% 		tRMSOut=calcAWeightedRMS(data.signalOut(t1:t2),data.params.sr);
% 	else 
% 		tRMSOut=0;
% 	end

%         load calibMic;  % gets micGain: wav rms/ Pa rms (Pa^-1)
	load('../../signals/leveltest/micRMS_100dBA.mat');  % Gives micRMS_100dBA: the rms the microphone should read when the sound is at 100 dBA SPL
%     text(xs(1)+0.05*range(xs),ys(2)-0.1*range(ys),['RMS(In)=',num2str(tRMSIn)]);
    xs=get(gca,'XLim'); ys=get(gca,'YLim');
	text(xs(1)+0.05*range(xs),ys(2)-0.21*range(ys),...
		['meanPeak soundLevel=',num2str(100+20*log10((meanPeakRMS/micRMS_100dBA))),' dBA SPL'],'FontSize',10);
    text(xs(1)+0.05*range(xs),ys(1)+0.24*range(ys),...
        ['Sentence duration=',num2str(sentDur),' sec'],'FontSize',10);
    text(xs(1)+0.05*range(xs),ys(1)+0.18*range(ys),sprintf('nSyls=%d',handles.nSyls),'FontSize',10);
    text(xs(1)+0.05*range(xs),ys(1)+0.12*range(ys),sprintf('meanSylDur=%f sec',meanSylDur),'FontSize',10);
    
    dataOut.meanPeakLevel=100+20*log10((meanPeakRMS/micRMS_100dBA));

	text(xs(1)+0.05*range(xs),ys(2)-0.27*range(ys),...
		['BGNoiseLevel=',num2str(100+20*log10((rmsBGNoise/micRMS_100dBA))),' dBA SPL'],'FontSize',10);

	text(xs(1)+0.05*range(xs),ys(2)-0.33*range(ys),...
		['SNR=',num2str(20*log10(meanPeakRMS/rmsBGNoise))],'FontSize',10);

%         load calibOutput;   
	% gives 'freq' and 'voltGains', measured at 'shanqing' M-audio configuration 
	% and -1.65 Phone volume knob.
%         mvg=mean(voltGains) * sqrt(2);    % mean voltage gain (V_rms / wavAmp_rms)
	
% 	soundLvOut=20*log10(tRMSOut*data.params.dScale/(dBSPL2WaveAmp(0,1000)/sqrt(2)));  % dBA SPL

	set(gcf,'CurrentAxes',handles.figIdDat(3));
	xs=get(gca,'XLim'); ys=get(gca,'YLim');
% 	if (~isnan(t1) && ~isnan(t2) && t1>0 && t2>0 && t2>t1)	
% 		plot([taxis(t1),taxis(t1)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
% 		plot([taxis(t2),taxis(t2)],[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
%     end
	text(xs(1)+0.05*range(xs),ys(2)-0.15*range(ys),...
		['dScale=',num2str(data.params.dScale)],'FontSize',10);    
% 	text(xs(1)+0.05*range(xs),ys(2)-0.2*range(ys),...
% 		['soundLevel=',num2str(soundLvOut),' dBA SPL'],'FontSize',10);


	set(gcf,'CurrentAxes',handles.figIdDat(4));
	cla;
    t_rms=data.rms(:,1);
    frameDur=data.params.frameLen/data.params.sr;
    taxis1=0:frameDur:(frameDur*(length(t_rms)-1));
    plot(taxis1, t_rms); hold on;
    xlabel('Time (sec)');
    ylabel('Signal RMS');
    plot([taxis1(1),taxis1(end)],repmat(threshRMS,1,2),'k-');
    ys=get(gca,'YLim');
    for i1=1:length(idx1)
        plot(repmat(taxis1(idx1(i1)),1,2),ys,'k--');
    end
    for i1=1:length(idx2)
        plot(repmat(taxis1(idx2(i1)),1,2),ys,'k-');
    end    
    hold on;
    set(gca,'XLim',[taxis1(1),taxis1(end)]);
    
% 	if (data.params.frameLen*idx1>=1 & data.params.frameLen*idx2<=length(taxis) & idx1>=1 & idx2 <= size(data.fmts,1))
% 		plot(taxis(data.params.frameLen*(idx1:idx2)),data.fmts(idx1:idx2,1),'k-','LineWidth',1.5);   hold on;
% 		plot(taxis(data.params.frameLen*(idx1:idx2)),data.fmts(idx1:idx2,2),'k-','LineWidth',1.5); 
% 		plot(taxis(data.params.frameLen*(idx1:idx2)),data.sfmts(idx1:idx2,1),'b-','LineWidth',1.5);
% 		plot(taxis(data.params.frameLen*(idx1:idx2)),data.sfmts(idx1:idx2,2),'b-','LineWidth',1.5);
% 		set(gca,'XLim',taxis(data.params.frameLen*([idx1,idx2])));
% 		set(gca,'YLim',[0,3000]);
% 		xs=get(gca,'XLim'); ys=get(gca,'YLim');
% 		plot(taxis(data.params.frameLen*([k1,k1])),[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;
% 		plot(taxis(data.params.frameLen*([k2,k2])),[ys(1),ys(2)],'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);  hold on;    
% 		xlabel('Time (sec)');
% 		ylabel('Formant freqs (Hz)');
% 	else
% 		cla;
% 	end

% 	set(gcf,'CurrentAxes',handles.figIdDat(5));
% 	cla;
% 	if (~isnan(k1) && ~isnan(k2) && k1>0 && k2>0 && t2>t1)
% 		plot(data.fmts(k1:k2,1),data.fmts(k1:k2,2),'b-','LineWidth',1.5);   hold on;
% 		plot(data.sfmts(k1:k2,1),data.sfmts(k1:k2,2),'b-','LineWidth',1.5);   hold off;
% 	end
% 	set(gca,'XLim',[0,2000]);
% 	set(gca,'YLim',[0,3000]);
% 	grid on;
% 	xlabel('F1 (Hz)');
% 	ylabel('F2 (Hz)');

% 	set(gcf,'CurrentAxes',handles.figIdDat(6));
% 	cla;

% 	set(gcf,'CurrentAxes',handles.figIdDat(7));
% 	cla;
% 	if (~isnan(k1) && ~isnan(k2) && k1>0 && k2>0 && t2>t1)
% 		plot(hz2mel(data.fmts(k1:k2,1)),hz2mel(data.fmts(k1:k2,2)),'b-','LineWidth',1.5);   hold on;
% 		plot(hz2mel(data.sfmts(k1:k2,1)),hz2mel(data.sfmts(k1:k2,2)),'b-','LineWidth',1.5);   hold off;
% 	end
% 	set(gca,'XLim',[0,1000]);
% 	set(gca,'YLim',[0,2000]);
% 	grid on;	
% 	xlabel('F1 (mel)');
% 	ylabel('F2 (mel)');    

% else
% 	pause(1e-3);
% end
% ~SCai: update the data monitor window
% --------------------------------------------------------------------------

dataOut.subject = handles.subject;
data.params.trialType= handles.trialType;
dataOut.params.name = handles.word;

handles.time2=clock;
if (handles.trigByScanner==0)
	timeToPause=handles.ITI-etime(handles.time2,handles.time1);
	if (timeToPause>0)
		pause(timeToPause);
    end
else
    pause(0.25);
end



return

% set(handles.msgh,'string','');
% if (handles.vumeterMode==1)
%     vumeter=handles.skin.vumeter;
% elseif (handles.vumeterMode==2)
%     vumeter=handles.skin.vumeter2;
% end
% mask=0.5*ones(size(vumeter));
% set(handles.rms_imgh,'Cdata',vumeter.*mask);
% set(handles.speed_imgh,'Cdata',vumeter.*mask);

%% --------------------------------------------------------------------------
function startRec(obj,event,handles)
% startRec displays a string and starts a timer object (handles.rect_timer)
% who's TimerFcn Callback (@stopRec) is called after a timeout period, set by
% the Value of the slider (handles.rec_slider)
% the StartFcn / StopFcn of the timer object starts/stops the recording
% (@soundDevice('start')/('stop')
% CDataPhrase=imread('utterimg/phrase.bmp');

fprintf('startRec\n')

handles.dataOut=[];
str=get(handles.strh,'string');
set(handles.rms_imgh,'Cdata',zeros(size(get(handles.rms_imgh,'Cdata'))));
set(handles.speed_imgh,'Cdata',zeros(size(get(handles.speed_imgh,'Cdata'))));

% set(handles.phrase_axes,'CData',CDataPhrase,'visible','on');

% clockStart=clock;
% clockStart(6) = clockStart(6) + get(handles.rec_slider,'Value');
% startat(handles.rec_timer, clockStart);
set(handles.strh,'string',str);

%% --------------------------------------------------------------------------
function stopRec(obj,event,handles)
% this function stops the recording if the timer object (handles.rec_timer)
% is still running. After the recording is stoped, checkData is executed,
% and the next and previous buttons are activated
fprintf('stopRec\n')
% if(strcmp(get(handles.rec_timer,'Running'),'on'))
%     stop(handles.rec_timer)
% end

handles.dataOut=checkData(getData,handles);


guidata(handles.UIrecorder, handles);

% set(handles.next,'enable','on')
% set(handles.prev,'enable','on')
% set(handles.next,'visible','on')
% set(handles.prev,'visible','on')
% if(get(handles.auto_btn,'Value')==get(handles.auto_btn,'Max'))
%     next_Callback(handles.UIrecorder,[],handles)
% end

%% --------------------------------------------------------------------------
% function soundDevice(obj,event,handles,action)
% % interface to teh external sounddevice
% switch(action)
%     case 'init'
%         TransShiftMex(0)
%     case 'start'
%         TransShiftMex(1)
%     case 'stop'
%         TransShiftMex(2)
%     otherwise,
% end

%% --------------------------------------------------------------------------
function dataOut= getData
% gets the data
dataOut=MexIO('getData');

%% --- Executes on button press in auto_btn.
% function auto_btn_Callback(hObject, eventdata, handles)
% hObject    handle to auto_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_btn

%%
function record(handles)
handles1 = guidata(handles.UIrecorder);
if isfield(handles1, 'lastData');
    handles.lastData = handles1.lastData;
    handles.lastSaveDataFN = handles1.lastSaveDataFN;
    handles.lastAsrDir = handles1.lastAsrDir;
end

if (handles.trigByScanner==1)
    set(handles.play,'userdata',1);
    uiwait(handles.UIrecorder);
else
    waitfor(handles.play,'userdata',1);
end

if (handles.debug==0)
    
    go=get(handles.UIrecorder,'UserData');
    if isequal(go,'nogo')
        return
    end

    handles.dataOut=[];
    guidata(handles.UIrecorder,handles);

    set(handles.strh,'visible','off');
    set(handles.msgh,'visible','off');

    handles.time1=clock;

    if (handles.trigByScanner==1)
        if (~isequal(handles.phase,'pre') && ~isequal(handles.phase,'pract1') && ~isequal(handles.phase,'pract2'))
            pause(handles.TA);
        else
            pause(0.25);
        end
    else
        pause(0.25);
    end
    
    % pause(0.1);
    if (handles.trialType==1 || handles.trialType==3)
        set(handles.pic_imgh,'cdata',handles.skin.fixation_nonpaced,'visible','on');
    elseif (handles.trialType==2 || handles.trialType==4)
        set(handles.pic_imgh,'cdata',handles.skin.fixation_paced,'visible','on');
    else
        set(handles.pic_imgh,'cdata',handles.skin.fixation,'visible','on');
    end
    
    drawnow;
    
    if (handles.trialType==1 || handles.trialType==3)   % Non-rhythmic
        paceStimParams=playRandToneSeq(8,'toneDur',handles.toneDur,'toneFreq',handles.toneFreq,...
            'toneAmp',handles.toneAmp,'toneRamp',handles.toneRamp,'totDur',handles.TPaceStim,...
            'seqDur',handles.meanSylDur*(8-1));
    else
        paceStimParams=playRandToneSeq(8,'uniform',...
            handles.toneDur,'toneFreq',handles.toneFreq,'toneAmp',handles.toneAmp,...
            'toneRamp',handles.toneRamp,'totDur',handles.TPaceStim,...
            'seqDur',handles.meanSylDur*(8-1));
    end
        
    
%     pause(0.18);
    set(handles.pic_imgh,'visible','off');
    drawnow;

    if (isequal(handles.word,'Ready...') || handles.trialType==-1)
        return
    end
    
    
    if (handles.trialType==1 || handles.trialType==2 || handles.trialType==3 || handles.trialType==4)
        set(handles.strh, 'string', handles.word,'visible','on');
        drawnow;
    end    

    if handles.interfaceMode == 0
        tic;
        TransShiftMex(3, 'fb', 1);

        fprintf('meanSylDur = %f sec\n',handles.meanSylDur);
        
        TransShiftMex(1);
        audapterStartupTime = toc;
%         pause(handles.trialLen);  % Changed 2008/06/18 to make the pause longer +1.5 --> +2.0        
    elseif handles.interfaceMode == 1 % winAudio for I/O (TODO)
        error('%s: interfaceMode == 1 has not been implemented yet', mfilname);
    elseif handles.interfaceMode == 2 % winAudio + sim 
        tic;
        data = load(handles.simDataFN);
        [dataOut, bRmsGood, bSpeedGood, t_rms] = checkData(data.data, handles, paceStimParams);
        dataLoadTime = toc;                
    else
        error('%s: unrecognized interfaceMode: %d', mfilename, handles.interfaceMode)
    end
    
    % --- Perform ASR on the latest speech data --- %    
    if isfield(handles, 'lastData')
        tic;
        asrPAlign = run_julian(handles.lastData, 'outDir', handles.lastAsrDir);

        set(0, 'CurrentFigure', handles.figIdDat(1));
        set(gcf, 'CurrentAxes', handles.figIdDat(7));
        cla;
        show_spectrogram(handles.lastData.signalIn, handles.lastData.params.sr, 'noFig');
        set(gca, 'XTick', []);

        plot_phn_align(asrPAlign);        
        title(['Last speech trial: "', handles.lastData.params.name, '"']);
        
        xs = get(gca, 'XLim');
        set(gcf, 'CurrentAxes', handles.figIdDat(4))
        set(gca, 'XLim', xs);

        % --- Calculate vowel timing --- %
        % TODO: code to deal with production error
        vidx = get_vowel_indices(handles.lastData.params.name);
        vts = get_vowel_t(asrPAlign, vidx);
        t_ivis = diff(vts);
        t_mean_ivi = mean(t_ivis);
        t_cv_ivis = std(t_ivis) / mean(t_ivis);
%             vts = get_vowel_t(asrPAlign, vidx, 'peakRMS', t_rms, handles.lastData.params.frameLen / handles.lastData.params.sr);

        handles = rmfield(handles, {'lastData', 'lastSaveDataFN', 'lastAsrDir'});
        handles.t_nVowels = length(vts);
        handles.t_mean_ivi = t_mean_ivi;
        handles.t_cv_ivis = t_cv_ivis;
        handles.t_ivis = t_ivis;
        
        asrTime = toc;
        fprintf(1, 'INFO: ASR on the last speech trial took %.2f seconds.\n', asrTime);
    else
        asrTime = 0;
    end   
    % --- ~Perform ASR on the latest speech data --- %
    
    
    
    if handles.interfaceMode == 0
        pause(max([0, handles.trialLen - audapterStartupTime - asrTime]));
    elseif handles.interfaceMode == 2
        pause(max([0, handles.trialLen - dataLoadTime - asrTime]));
    end
    
    if handles.interfaceMode == 0
        if get(handles.play, 'userdata')==0 % in pause mode
            record(handles); % re-do recording
        end
        TransShiftMex(2);

        [dataOut, bRmsGood, bSpeedGood, t_rms] = checkData(getData, handles, paceStimParams);
    end    
    
    % --- Store speech data of the current trial for later ASR --- %
    if handles.trialType == 1 || handles.trialType == 2
        handles.lastData = dataOut;
        handles.lastSaveDataFN = handles.saveDataFN;
        handles.lastAsrDir = handles.asrDir;
    end
    
    
    % --- ~Store speech data of the current trial for later ASR --- %
    
    

    bRmsRepeat=handles.bRmsRepeat;
    bSpeedRepeat=handles.bSpeedRepeat;
    if (handles.trialType==3 || handles.trialType==4)
        if (bRmsRepeat==1)
            bRmsRepeat=0;
        end
        if (bSpeedRepeat==1)
            bSpeedRepeat=0;
        end
    end


    if ((~bRmsGood && bRmsRepeat) || (~bSpeedGood && bSpeedRepeat))
        record(handles)    
    else
    % data is saved as UserData in the fig handle (wicht is the signal for
    % the host function to launch the next single trial
    set(handles.UIrecorder,'UserData',dataOut)
    end
    % end
    if (handles.trigByScanner==0)
        pause(0.25);
    end
    set(handles.strh,'visible','off');
    % set(handles.pic_imgh,'cdata',handles.skin.fixation);
    set(handles.pic_imgh,'visible','off');
else
    dataOut=struct;
    dataOut.signalIn=[]; dataOut.signalOut=[];
    dataOut.rms=[];
    set(handles.UIrecorder,'UserData',dataOut);
end

guidata(handles.UIrecorder, handles);

return

% --- Executes on button press in button_next.
% function button_next_Callback(hObject, eventdata, handles)
% % hObject    handle to button_next (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% if (isfield(handles,'nextMessage'))
%     set(handles.msgh_imgh,'CData',handles.nextMessage,'visible','on');
% end
% set(handles.button_next,'visible','off');
% 
% set(handles.play,'visible','on');

% function dFaces=getDFaces(fileMask)
%     dFaces=struct;
%     dFaces.d=dir(fileMask);
%     dFaces.sex=cell(1,length(dFaces.d));
%     dFaces.subjID=nan(1,length(dFaces.d));
%     
%     for n=1:length(dFaces.d)
%         fn=dFaces.d(n).name;
%         idx=strfind(fn,'.bmp');
%         dFaces.sex{n}=dFaces.d(n).name(idx-1);
%         idx1=strfind(fn,'-s');
%         idx2=strfind(fn,['-',dFaces.sex{n}]);
%         dFaces.subjID(n)=str2num(dFaces.d(n).name(idx1+2:idx2-1));
%     end
% return