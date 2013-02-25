function [dur,meanPeakRMS,threshRMS,idx1,idx2]=analyzeSentence(in,varargin)
%% Config
silWin=150e-3;
minSegDur=200e-3;
peakNum=6;

bPlot=0;
if ~isempty(fsic(varargin,'plot'))
    bPlot=1;
end

%%
if isstruct(in)
    data=in;
else
    if isfile(in)
        load(in);   % gives data
    end
end

%%
t_rms=data.rms(:,1);
t_rms=mva(t_rms,33);
frameDur=data.params.frameLen/data.params.sr;
taxis0=0:frameDur:frameDur*(length(t_rms)-1);

silRMS=t_rms(taxis0<silWin);

threshRMS=mean(silRMS)+5*std(silRMS);

segAbvRMS=nan(0,2);

if t_rms(1)>threshRMS
    idx=min(find(t_rms<=threshRMS));
    t_rms=t_rms(idx:end);    
end

stat=0;
idx1=find(t_rms(1:end-1)<=threshRMS & t_rms(2:end)>threshRMS);
idx2=find(t_rms(1:end-1)>threshRMS & t_rms(2:end)<=threshRMS);

if length(idx1)>length(idx2)
    idx1=idx1(1:length(idx2));
end

segLens=(idx2-idx1)*frameDur;
idx1=idx1(segLens>minSegDur);
idx2=idx2(segLens>minSegDur);
segLens=segLens(segLens>minSegDur);

% dur=sum(segLes);
if ~isempty(idx1) && ~isempty(idx2)
    dur=(idx2(end)-idx1(1))*frameDur;
else
    dur=NaN;
end

peakRMSs=[];
for i1=1:length(idx1)
    rms_seg=t_rms(idx1(i1):idx2(i1));
    d_rms_seg=diff(rms_seg);
    idx_peaks=find(d_rms_seg(1:end-1)>0 & d_rms_seg(2:end)<0);
    peakRMSs=[peakRMSs,rms_seg(idx_peaks+1)];    
end
peakRMSs=sort(peakRMSs);
if length(peakRMSs)>=peakNum
    peakRMSs=peakRMSs(end-peakNum+1:end);    
end

meanPeakRMS=rms(peakRMSs);

%% Visualization
if bPlot
    figure;
    plot(taxis0,t_rms); hold on;
    plot([taxis0(1),taxis0(end)],repmat(threshRMS,1,2),'k-');
    set(gca,'XLim',[taxis0(1),taxis0(end)]);
    ys=get(gca,'YLim');
    for i1=1:length(idx1)
        plot(repmat(taxis0(idx1(i1)),1,2),ys,'k-');
        plot(repmat(taxis0(idx2(i1)),1,2),ys,'k-');
    end
    plot([taxis0(1),taxis0(end)],repmat(meanPeakRMS,1,2),'r-');
end
return