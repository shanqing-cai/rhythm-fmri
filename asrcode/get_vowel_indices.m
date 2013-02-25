function vi = get_vowel_indices(sent)
%% CONFIG
vifn = 'E:/speechres/rhythm-fmri/asr/vowel_indices.mat';

%%
if ~isfile(vifn)
    error('Vowel index file missing: %s', vifn);
end

vidat = load(vifn);
if isempty(fsic(vidat.stcs, sent))
    vi = [];
else
    vi = vidat.vowelIndices{fsic(vidat.stcs, sent)};
end


return