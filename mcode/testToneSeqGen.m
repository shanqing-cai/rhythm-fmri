function testToneSeqGen()
%% Parameters
maxNTones=64;
p.tsgNTones=8;
p.tsgToneDur=[repmat(25e-3,1,p.tsgNTones),zeros(1,maxNTones-p.tsgNTones)];
p.tsgToneFreq=[repmat(6000,1,p.tsgNTones),zeros(1,maxNTones-p.tsgNTones)];
p.tsgToneAmp=[repmat(0.2,1,p.tsgNTones),zeros(1,maxNTones-p.tsgNTones)];
p.tsgToneRamp=[repmat(5e-3,1,p.tsgNTones),zeros(1,maxNTones-p.tsgNTones)];
p.tsgInt=[repmat(0.333,1,p.tsgNTones),zeros(1,maxNTones-p.tsgNTones)];

%%
TransShiftMex(3,'tsgntones',p.tsgNTones);
TransShiftMex(3,'tsgtonedur',p.tsgToneDur);
TransShiftMex(3,'tsgtonefreq',p.tsgToneFreq);
TransShiftMex(3,'tsgtoneamp',p.tsgToneAmp);
TransShiftMex(3,'tsgtoneramp',p.tsgToneRamp);
TransShiftMex(3,'tsgint',p.tsgInt);
TransShiftMex(3,'wgtime',0);

totDur=sum(p.tsgInt);

TransShiftMex(13);
pause(totDur);
TransShiftMex(2);
return