function phaseScript = ...
    genPhaseScript_fmri(stage, nReps, trialTypes, ...
                        trainWords, nSyls, randReps, trigByScanner)
	phaseScript=struct();
	phaseScript.nReps=nReps;
    phaseScript.nTrials=0;

    symbols='@#$%^&*()_+=-<>/\|[]{}';
    nSymbols=length(symbols);
    
    % --- Ensuring that the sentence sets are identical between the N and R
    % conditions ---
    nSpeechN = numel(find(trialTypes == 1));
    nSpeechR = numel(find(trialTypes == 2));    
    if nSpeechN ~= nSpeechR
        error('Unequal number of non-rhythmic and rhytmic trials');
    end
    nSpeech = nSpeechN * 2;
    
    trainWordsUsed = {};
    nSylsUsed = [];
    for i1 = 1 : nReps
        trainWordsUsed = [trainWordsUsed, ...
                          repmat(trainWords((i1 - 1) * nSpeech + 1 : (i1 - 1) * nSpeech + nSpeechN), 1, 2)];
        nSylsUsed = [nSylsUsed, repmat(nSyls((i1 - 1) * nSpeech + 1 : (i1 - 1) * nSpeech + nSpeechN), 1, 2)];
    end

    twCnt=1;
    ptwCnt=1;
    for n=1:nReps
        bt=[trialTypes];
        trainWordsRep = trainWordsUsed((n - 1) * nSpeech + 1 : n * nSpeech);
        nSylsRep = nSylsUsed((n - 1) * nSpeech + 1 : n * nSpeech);
        
%         pseudoWordsUsed=pseudoWords(randperm(length(trainWords)));
%             testWordsUsed2=testWords(randperm(length(testWords)));

        idxrp = randperm(length(bt));        
        bt = bt(idxrp);
        
        oneRep=struct;
        oneRep.trialOrder=[];
        oneRep.word=cell(1,0);
        oneRep.nSyls=[];
        
        nCnt = 1;
        rCnt = nSpeechN;
        bCnt = 1;
        
        for m=1:length(bt)
            oneRep.trialOrder=[oneRep.trialOrder,bt(m)];
            if (bt(m) == 1 || bt(m) == 2)
                if bt(m) == 1 % -- Non-rhythmic -- %
                    oneRep.word{length(oneRep.word) + 1} = trainWordsRep{nCnt};
                    oneRep.nSyls(end+1) = nSyls(nCnt);
                    
                    nCnt = nCnt + 1;
                else % -- Rhythmic -- %
                    oneRep.word{length(oneRep.word) + 1} = trainWordsRep{rCnt + nSpeechN};
                    oneRep.nSyls(end+1) = nSyls(rCnt + nSpeechN);
                    
                    rCnt = rCnt - 1;
                end
                
            elseif (bt(m) == 3 || bt(m) == 4)
                t_sent = trainWordsRep{bCnt};
                
                for i1=1:length(t_sent)
                    if ~isequal(t_sent(i1),' ')
                        t_sent(i1) = symbols(round(rand * (nSymbols - 1)) + 1);
                    end
                end
                oneRep.word{length(oneRep.word)+1} = t_sent;
                oneRep.nSyls(end+1) = nSyls(bCnt);
                
                bCnt = bCnt + 1;
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