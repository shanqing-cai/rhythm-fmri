function genRandIntervals
%% Config
nInt=7;
nSeq=500;
saveFN='randIntervals_7.mat';

bSave=1;
%%
seqs=nan(0,nInt);
for i1=1:nSeq
    t_int=rand(1,nInt)*0.75+0.25;
    t_int=t_int/sum(t_int);
    seqs=[seqs;t_int];
end

%% visualization
figure;
for i1=1:size(seqs,1)
    t_int=seqs(i1,:);
    c_int=cumsum(t_int);
    for i2=1:nInt        
        plot(repmat(c_int(i2),1,2),[i1-1,i1],'k-');
        hold on;
    end
end

%% Saving
if bSave
    save(saveFN,'seqs');
    fprintf('seqs saved in %s.\n',saveFN);
end
return