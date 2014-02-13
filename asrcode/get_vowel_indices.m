function vi = get_vowel_indices(sent)
%% CONFIG
vifn = '../asr/vowel_indices.mat';

%%
if ~isfile(vifn)
    error('Vowel index file missing: %s', vifn);
end

vidat = load(vifn);
if isempty(fsic(vidat.stcs, sent))
    fprintf(2, 'WARNING: sentence "%s" is not found in vowel_indices file: %s\n', ...
            sent, vifn);
    vi = [];
else
    vi = vidat.vowelIndices{fsic(vidat.stcs, sent)};
end


return