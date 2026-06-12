function pSig = timePermClusterMulti(signalCell, clusterargs)
    arguments
        signalCell  % Cell array of signals (conditions x trials x time series)
        clusterargs.statstype = 'f-stat';
        clusterargs.nPerm = 1000; 
        clusterargs.numTail = 1;
        clusterargs.pThresh = 0.05;
    end
    
    pThresh = clusterargs.pThresh;
    nPerm = clusterargs.nPerm;
    
    % Ensure all conditions have the same number of time points
    nConditions = numel(signalCell);
    timeLength = size(signalCell{1}, 2);
    assert(all(cellfun(@(x) size(x, 2), signalCell) == timeLength), 'Conditions have different time lengths.');
    
    % Compute F-statistic across conditions at each time point
    statMetric = statfun(signalCell, type=clusterargs.statstype);
    
    % Combine signals from all conditions for permutation testing
    combSignals = vertcat(signalCell{:});
    nComb = size(combSignals, 1);
    fprintf(1, '\n>> Calculating permutation statistic...\n');
    permStat = zeros(nPerm, timeLength);
    
    parfor iPerm = 1:nPerm
        permIds = randperm(nComb);
        permSignals = cell(nConditions, 1);
        startIdx = 1;
        
        % Split permuted signals back into conditions
        for iCond = 1:nConditions
            nTrials = size(signalCell{iCond}, 1);
            permSignals{iCond} = combSignals(permIds(startIdx:startIdx + nTrials - 1), :);
            startIdx = startIdx + nTrials;
        end
        
        permStat(iPerm, :) = statfun(permSignals, type=clusterargs.statstype);
        fprintf(1, '.');
    end
    fprintf(1, '] done\n');
    
    % Compute p-values based on permutation statistics
    pStatOne = sum(permStat > statMetric, 1) ./ nPerm;
    pStatTwo = sum(abs(permStat) > abs(statMetric), 1) ./ nPerm;
    pStatOne(isnan(statMetric)) = nan;
    pStatTwo(isnan(statMetric)) = nan;
    pStatOne = adjustPVals(pStatOne, nPerm);
    pStatTwo = adjustPVals(pStatTwo, nPerm);
    
    fprintf(1, '\n>> Calculating permutation p-values...\n');
    for iPerm = 1:nPerm
        permRecur = setdiff(1:nPerm, iPerm);
        pStatPermOne(iPerm,:) = sum(permStat(permRecur,:) > permStat(iPerm,:), 1) ./ (nPerm - 1);
        pStatPermOne(iPerm,isnan(permStat(iPerm,:)))=nan;
        pStatPermOne(iPerm,:) = adjustPVals(pStatPermOne(iPerm,:), nPerm - 1);
        pStatPermTwo(iPerm,:) = sum(abs(permStat(permRecur,:)) > abs(permStat(iPerm,:)), 1) ./ (nPerm - 1);
        pStatPermTwo(iPerm,isnan(permStat(iPerm,:)))=nan;
        pStatPermTwo(iPerm,:) = adjustPVals(pStatPermTwo(iPerm,:), nPerm - 1);
        fprintf(1,'.');
    end
    fprintf(1, '] done\n');
    
    if(clusterargs.numTail == 1)
        fprintf(1,'\n>> Performing one-sided cluster correction');
        [pValsRaw, actClust] = timePermClusterAfterPermPValues(pStatOne, pStatPermOne, pThresh);
    else
        fprintf(1,'\n>> Performing two-sided cluster correction');
        [pValsRaw, actClust] = timePermClusterAfterPermPValues(pStatTwo, pStatPermTwo, pThresh);
    end
    
    h_sig_05 = zeros(1,length(pValsRaw));
    h_sig_01 = zeros(1,length(pValsRaw));
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            h_sig_05(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
            if actClust.Size{iClust}>actClust.perm99
                h_sig_01(actClust.Start{iClust}: ...
                    actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                    = 1;
            end
        end
    end
    
    pSig.pVals = pValsRaw;
    pSig.clust = actClust;
    pSig.h_sig_05 = h_sig_05;
    pSig.h_sig_01 = h_sig_01;

end

function statMetric = statfun(signalCell, statargs)
arguments
    signalCell % Cell array of signals (conditions x trials x time series)
    statargs.type {mustBeMember(statargs.type, {'mean-sub', 't-stat', 'f-stat'})} = 'f-stat'
end

switch statargs.type
    case 'mean-sub'
        % Compute mean difference between first two conditions as an example (not applicable for >2 conditions)
        assert(numel(signalCell) == 2, 'Mean subtraction only valid for two conditions.');
        statMetric = mean(signalCell{1}, 1) - mean(signalCell{2}, 1);

    case 't-stat'
        % Compute t-statistic between first two conditions as an example (not applicable for >2 conditions)
        assert(numel(signalCell) == 2, 'T-test only valid for two conditions.');
        [~, ~, ~, stats] = ttest2(signalCell{1}, signalCell{2});
        statMetric = stats.tstat;

    case 'f-stat'
        % Compute F-statistic across all conditions at each time point
        nConditions = numel(signalCell);
        dataMatrix = vertcat(signalCell{:}); % Combine all trials into one matrix (trials x time series)
        groupLabels = cell2mat(arrayfun(@(x) repmat(x, size(signalCell{x}, 1), 1), (1:nConditions)', 'UniformOutput', false)); % Group labels
        
        %statMetric = zeros(1, size(dataMatrix, 2));
        statMetric = computeFStatistic(dataMatrix, groupLabels);
        % parfor tIdx=1:size(dataMatrix,2)
        %     [~, tbl] = anova1(dataMatrix(:, tIdx), groupLabels,'off');
        %     statMetric(tIdx) = tbl{2,5}; % Extract F-statistic from ANOVA table
        % end
end

end

function [adjustedPVals] = adjustPVals(pVals, nPerm)
adjustedPVals = max(min(pVals, 1 - 1/nPerm), 1/nPerm);

end
