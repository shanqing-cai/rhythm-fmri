function testFeedback(configFN)
subject = read_subject_config(configFN);

subject.mouthMicDist = input('Please input mouth-mic distance (cm): ');

p = getTSMDefaultParams(subject.sex, ...
                        'closedLoopGain', subject.closedLoopGain,...
                        'trialLen', subject.trialLen, ...
                        'mouthMicDist', subject.mouthMicDist);
MexIO('init', p);

TransShiftMex(3, 'fb', 1);
MexIO('reset');

TransShiftMex(1);
foo = input('Hit Enter to stop ...');
TransShiftMex(2);
return

