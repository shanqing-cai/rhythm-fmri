function phaseScript = ...
    genPhaseScript_fmri(stage, nReps, trialTypes, ...
                        trainWords, nSyls, randReps, trigByScanner)
	phaseScript=struct();
	phaseScript.nReps=nReps;
    phaseScript.nTrials=0;

    symbols='@#$%^&*()_+=-<>/\|[]{}';
    nSymbols=length(symbols);

    twCnt=1;
    ptwCnt=1;
    for n=1:nReps
        bt=[trialTypes];
        trainWordsUsed=trainWords;
%         pseudoWordsUsed=pseudoWords(randperm(length(trainWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));            
       
        bt=bt(randperm(length(bt)));
        oneRep=struct;
        oneRep.trialOrder=[];
        oneRep.word=cell(1,0);
        oneRep.nSyls=[];
        cntTW=1;
        for m=1:length(bt)
            oneRep.trialOrder=[oneRep.trialOrder,bt(m)];
            if (bt(m)==1 || bt(m)==2)					                
                oneRep.word{length(oneRep.word)+1}=trainWordsUsed{twCnt};
                oneRep.nSyls(end+1)=nSyls(twCnt);
                twCnt=twCnt+1;
            elseif (bt(m)==3 || bt(m)==4)
                t_sent=trainWordsUsed{ptwCnt};                
                for i1=1:length(t_sent)
                    if ~isequal(t_sent(i1),' ')
                        t_sent(i1)=symbols(round(rand*(nSymbols-1))+1);
                    end
                end
                oneRep.word{length(oneRep.word)+1}=t_sent;
                oneRep.nSyls(end+1)=nSyls(ptwCnt);
%                 oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW);
%                 oneRep.word{length(oneRep.word)+1}=pseudoWordsUsed(cntTW+1);
%                 cntTW=cntTW+2;
                ptwCnt=ptwCnt+1;
            end
        end
        
%         if ~trigByScanner % -- Behavioral session: remove the no-speech trials -- %
%             idxKeep = find(oneRep.trialOrder <= 2);
%             oneRep.trialOrder = oneRep.trialOrder(idxKeep);
%             oneRep.word = oneRep.word(idxKeep);
%             oneRep.nSyls = oneRep.nSyls(idxKeep);
%         end
% 
        phaseScript.(['rep',num2str(n)])=oneRep;
        phaseScript.nTrials=phaseScript.nTrials+length(oneRep.trialOrder);

        if (isequal(stage(1:3),'run') && n==nReps && trigByScanner == 1)
            phaseScript.nTrials=phaseScript.nTrials+1;
            
            idx0=find(phaseScript.(['rep',num2str(n)]).trialOrder==3,1);            
            phaseScript.(['rep',num2str(n)]).trialOrder(end+1)=phaseScript.(['rep',num2str(n)]).trialOrder(idx0);
            phaseScript.(['rep',num2str(n)]).word{end+1}=phaseScript.(['rep',num2str(n)]).word{idx0};
            phaseScript.(['rep',num2str(n)]).nSyls(end+1)=phaseScript.(['rep',num2str(n)]).nSyls(idx0);
        end
    end

return