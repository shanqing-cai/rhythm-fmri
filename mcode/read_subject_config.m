function subject = read_subject_config(configFN)
subject = read_struct_from_text(configFN);

subject.paceStim = struct('meanSylDur', subject.paceStim_meanSylDur, ...
                          'toneDur', subject.paceStim_toneDur, ...
                          'toneFreq', subject.paceStim_toneFreq, ...
                          'toneAmp', subject.paceStim_toneAmp, ...
                          'toneRamp', subject.paceStim_toneRamp);
subject = rmfield(subject, {'paceStim_meanSylDur', 'paceStim_toneDur', ...
                            'paceStim_toneFreq', 'paceStim_toneAmp', 'paceStim_toneRamp'});
                       
if subject.trigByScanner == 1
	subject.showProgress		=0;
	subject.showPlayButton      =0;
else
	subject.showProgress		=1;
	subject.showPlayButton      =1;
end

subject.ITI					=subject.TA + subject.TPaceStim + subject.TVisStim;
return