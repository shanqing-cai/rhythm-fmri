function [sents, nSyls] = getRandSentences(nSent)
%%
% load('sentencePool.mat');   % gives 'orig_sents' and 'orig_nSyls'
load('IEEE_8vwls_nics1-3.mat'); 

%% Discard the sentences with flapped /t/ in A.E.
bKeep = ones(1, length(selStcs));
for i1 = 1 : numel(bKeep)
    bKeep(i1) = isempty(strfind(lower(selStcs{i1}), 'dirty')) & ...
                isempty(strfind(lower(selStcs{i1}), 'thirty')) & ...
                isempty(strfind(lower(selStcs{i1}), 'kitten'));
end
selStcs = selStcs(find(bKeep));
%%

orig_sents = selStcs;
orig_nSyls = 8 * ones(size(orig_sents));

%% 
if nSent > length(orig_sents)
    orig_sents = repmat(orig_sents, 1, ceil(nSent / length(orig_sents)));
    orig_nSyls = repmat(orig_nSyls, 1, ceil(nSent / length(orig_nSyls)));
end

%% Random permutation and selection

if nSent>length(orig_sents)
    fprintf('%s: Error: nSent larget than the number of available sentneces.\n',mfilename);
    sents=[];
    return
end
idxperm=randperm(length(orig_sents));
reshuf_sents=orig_sents(idxperm);
reshuf_nSyls=orig_nSyls(idxperm);
sents=reshuf_sents(1:nSent);
nSyls=reshuf_nSyls(1:nSent);

return