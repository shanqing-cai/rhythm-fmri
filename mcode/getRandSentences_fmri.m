function [sents, nSyls] = getRandSentences_fmri(expt)
%% Ad hoc!
NSYLS_DEFAULT = 8;

%%
load('IEEE_8vwls_nics1-3.mat');  % Gives selStcs

%%
all_phases = fields(expt.script);

sents = {};
nSyls = [];

nSpeechN = numel(find(expt.trialTypes == 1));
nSpeechR = numel(find(expt.trialTypes == 2));
if nSpeechN ~= nSpeechR
    error('Unequal number of non-rhythmic and rhytmic trials');
end

nSpeech = nSpeechN * 2;
for i1 = 1 : numel(all_phases)
    ph = all_phases{i1};
    nst = expt.script.(ph).nReps * nSpeech;
    if isequal(ph, 'pre') || isequal(ph, 'pract1') || isequal(ph, 'pract2') || (length(ph) > 5 && isequal(ph(1 : 5), 'inter'))
        idxrp = [];
        while length(idxrp) < nst;
            idxrp = [idxrp, randperm(length(selStcs))];
        end
        
        sents = [sents, selStcs(idxrp(1 :  nst))];
        nSyls = [nSyls, repmat(NSYLS_DEFAULT, size(sents))];
    elseif length(ph) > 3 && isequal(ph(1 : 3), 'run')
        if mod(nst / 2, length(selStcs)) ~= 0
            error('The number of rhythmic (or natural) speech trials in phase %s of the functional run is not an integer multiple of the number of sentences in the pool (%d)', ...
                  ph, length(selStcs));              
        end
        
        k = nst / length(selStcs) / 2;
        idxrp = [];
        for i1 = 1 : k
            idxrp = [idxrp, randperm(length(selStcs))];
        end
        
        t_sents = selStcs(idxrp);
        t_sents = reshape(t_sents, [expt.script.(ph).nReps, nSpeechN]);
        t_sents = repmat(t_sents, 1, 2);
        t_sents = reshape(t_sents', [1, nst]);
        t_nsyls = zeros(size(t_sents));
        
        sents = [sents, t_sents];
        nSyls = [nSyls, t_nsyls];
    end
end

return