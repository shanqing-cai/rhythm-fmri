function selStcs = select_IEEE_sents(nv, nics, outFN)
% 
% Inputs: nv: number of vowels
%         nics: [1x2] vector, the number of intervening consonants between
%         adjacent vowels will be >= nics[1] and <= nics[2]
%% CONFIG
ASR_DATA_DIR = 'E:/speechres/rhythm-fmri/asr';
vowelIndicesFN = fullfile(ASR_DATA_DIR, 'vowel_indices.mat');

%%
vidat = load(vowelIndicesFN);

ns = numel(vidat.stcs);
nVowels = nan(1, ns);
minIntCons = nan(1, ns);
maxIntCons = nan(1, ns);

for i1 = 1 : ns
    nVowels(i1) = length(vidat.vowelIndices{i1});
    minIntCons(i1) = min(diff(vidat.vowelIndices{i1})) - 1;
    maxIntCons(i1) = max(diff(vidat.vowelIndices{i1})) - 1;
end

idxs = find(nVowels == nv & minIntCons >= nics(1) & maxIntCons <= nics(2));

selStcs = vidat.stcs(idxs);

%% Save to output file
save(outFN, 'selStcs');
if isfile(outFN)
    fprintf(1, 'INFO: Saved results to file %s\n', outFN);
end
return