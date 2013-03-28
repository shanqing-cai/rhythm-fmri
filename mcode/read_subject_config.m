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

subject.showRhythmicityFB_phases = proc_phase_items(subject.showRhythmicityFB_phases);
subject.showRateFB_phases = proc_phase_items(subject.showRateFB_phases);
subject.showIntFB_phases = proc_phase_items(subject.showIntFB_phases);

% subject.showRhythmicityWarn_phases = proc_phase_items(subject.showRhythmicityWarn_phases);
subject.showRateWarn_phases = proc_phase_items(subject.showRateWarn_phases);
subject.showIntWarn_phases = proc_phase_items(subject.showIntWarn_phases);

subject.ITI					=subject.TA + subject.TPaceStim + subject.TVisStim;
return

function cellout = proc_phase_items(strin)
if isequal(lower(strin), 'none')
    cellout = {};
else
    cellout = strsplit(strin, ',');
end
return