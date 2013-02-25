function [sents,nSyls]=getRandSentences(nSent)
%%
load('sentencePool.mat');   % gives 'orig_sents' and 'orig_nSyls'

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