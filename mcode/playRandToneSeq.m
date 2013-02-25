function p = playRandToneSeq(nTones, varargin)
%%
randSeqFN=['randIntervals_',num2str(nTones-1),'.mat'];

seqDur=2.5;

maxNTones=64;
toneDur=25e-3;  % sec
toneFreq=6e3;   % Hz
toneAmp=0.2;    % Wave peak
toneRamp=5e-3;  % sec

if ~isempty(fsic(varargin,'toneDur'))
    toneDur=varargin{fsic(varargin,'toneDur')+1};   
end
if ~isempty(fsic(varargin,'toneFreq'))
    toneFreq=varargin{fsic(varargin,'toneFreq')+1};   
end
if ~isempty(fsic(varargin,'toneAmp'))
    toneAmp=varargin{fsic(varargin,'toneAmp')+1};   
end
if ~isempty(fsic(varargin,'toneRamp'))
    toneRamp=varargin{fsic(varargin,'toneRamp')+1};   
end

totDur=NaN;
if ~isempty(fsic(varargin,'totDur'))
    totDur=varargin{fsic(varargin,'totDur')+1};   
end
if ~isempty(fsic(varargin,'seqDur'))
    seqDur=varargin{fsic(varargin,'seqDur')+1}; 
end

bUniform=0;
if ~isempty(fsic(varargin,'uniform'))
    bUniform=1;
end

%%
load(randSeqFN);    % gives seqs
nSeqs=size(seqs,1);
% nTones=size(seqs,2)+1;
k=round(rand*(nSeqs-1))+1;

t_int=seqs(k,:);
t_int=t_int/sum(t_int); % Just for saveguard

if bUniform==1
    t_int=repmat(1/(nTones-1),1,nTones-1);
end

p.tsgNTones=nTones;
p.tsgToneDur=[repmat(toneDur,1,nTones),zeros(1,maxNTones-nTones)];
p.tsgToneFreq=[repmat(toneFreq,1,nTones),zeros(1,maxNTones-nTones)];
p.tsgToneAmp=[repmat(toneAmp,1,nTones),zeros(1,maxNTones-nTones)];
p.tsgToneRamp=[repmat(toneRamp,1,nTones),zeros(1,maxNTones-nTones)];
p.tsgInt=[t_int*seqDur,toneDur+0.05,zeros(1,maxNTones-nTones-1)];

%%
bPrompt = 0;
TransShiftMex(3, 'tsgntones', p.tsgNTones, bPrompt);
TransShiftMex(3, 'tsgtonedur', p.tsgToneDur, bPrompt);
TransShiftMex(3, 'tsgtonefreq', p.tsgToneFreq, bPrompt);
TransShiftMex(3, 'tsgtoneamp', p.tsgToneAmp, bPrompt);
TransShiftMex(3, 'tsgtoneramp', p.tsgToneRamp, bPrompt);
TransShiftMex(3, 'tsgint', p.tsgInt, bPrompt);
TransShiftMex(3, 'wgtime', 0, bPrompt);

if isnan(totDur)
    totDur=sum(p.tsgInt);    
end
p.totDur=totDur;
p.seqDur=seqDur;

TransShiftMex(13);
pause(totDur);
TransShiftMex(2);
return