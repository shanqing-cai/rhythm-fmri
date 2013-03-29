function phaseScript = genPhaseScript_behav(phase, nReps, sInfo)
load('IEEE_8vwls_nics1-3.mat');  % Gives selStcs

phasesNoPert = {'pre', 'pract1', 'pract2'};

%% Check the sanity of trial number settings
if sInfo.nTrialsF1UpPerRep + sInfo.nTrialsDecelPerRep + sInfo.nTrialsOtherPerRep >= sInfo.nTrialsPerRep
    error('sInfo.nTrialsPerRep is too small');
end
nTrialsNonePerRep = sInfo.nTrialsPerRep - (sInfo.nTrialsF1UpPerRep + sInfo.nTrialsDecelPerRep + sInfo.nTrialsOtherPerRep);

if mod(sInfo.nTrialsPerRep, 2) ~= 0
    error('nTrialsPerRep is not an even interger');
end
if mod(sInfo.nTrialsF1UpPerRep, 2) ~= 0
    error('nTrialsF1UpPerRep is not an even interger');
end
if mod(sInfo.nTrialsDecelPerRep, 2) ~= 0
    error('nTrialsDecelPerRep is not an even interger');
end
if mod(sInfo.nTrialsOtherPerRep, 2) ~= 0
    error('nTrialsOtherPerRep is not an even interger');
end

%% Locate the pertSent in selStcs
idxPertSent = -1;
for i1 = 1 : numel(selStcs)
    if isequal(lower(strrep(sInfo.pertSent, '_', ' ')), lower(selStcs{i1}))
        idxPertSent = i1; 
        break;
    end
end

if idxPertSent < 0
    error('Failed to located pert sentence: %s', sInfo.pertSent);
end

rndStcs = selStcs(setxor(1 : length(selStcs), idxPertSent));
nRnd = length(rndStcs);

%%
% -- Pert type: 0 - None
% --            1 - F1Up
% --            2 - Decel
% --            4 - Other 

phaseScript = struct('nReps', nReps, 'nTrials', 0);
for i1 = 1 : nReps
    thisRep = struct;
    thisRep.trialOrder = [1 * ones(1, sInfo.nTrialsPerRep / 2), 2 * ones(1, sInfo.nTrialsPerRep / 2 - 1)];
    thisRep.pertType = [0 * ones(1, nTrialsNonePerRep / 2), ...
                        1 * ones(1, sInfo.nTrialsF1UpPerRep / 2), ...
                        2 * ones(1, sInfo.nTrialsDecelPerRep / 2), ...
                        4 * ones(1, sInfo.nTrialsOtherPerRep / 2)];
    thisRep.pertType = repmat(thisRep.pertType, 1, 2);
    thisRep.pertType = thisRep.pertType(1 : end - 1);
    
    bGood = 0;
    while ~bGood
        idxrp = randperm(length(thisRep.trialOrder));
        t_pertType = thisRep.pertType(idxrp);
        idxIsPert = find(t_pertType >= 1 & t_pertType <= 2);
        bGood = isempty(find(diff(idxIsPert) == 1, 1));
    end
    
    thisRep.trialOrder = [thisRep.trialOrder(idxrp), 2];
    thisRep.pertType = [thisRep.pertType(idxrp), 4];
    
    thisRep.word = cell(1, length(thisRep.trialOrder));
    thisRep.nSyls = 8 * ones(1, length(thisRep.trialOrder));
    for i2 = 1 : numel(thisRep.trialOrder)
        if thisRep.pertType(i2) < 4
            thisRep.word{i2} = selStcs{idxPertSent};
        else
            iRand = randperm(nRnd);
            thisRep.word{i2} = rndStcs{iRand(1)};
        end
    end
    
    if ~isempty(fsic(phasesNoPert, phase))
         thisRep.pertType(thisRep.pertType >= 1 & thisRep.pertType <= 2) = 0;
    end
    
    phaseScript.(['rep', num2str(i1)]) = thisRep;
    phaseScript.nTrials = phaseScript.nTrials + length(thisRep.trialOrder);
    
    
end
return