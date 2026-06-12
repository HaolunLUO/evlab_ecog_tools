function statMetric = computeFStatistic(dataMatrix, groupLabels)
    % Compute F-statistic for each time point across trials
    %
    % Parameters:
    %   dataMatrix: trials x time points matrix of data
    %   groupLabels: vector specifying the group/condition for each trial
    %
    % Returns:
    %   statMetric: 1 x time points vector of F-statistics

    % Validate inputs
    assert(size(dataMatrix, 1) == numel(groupLabels), ...
        'Number of rows in dataMatrix must match the length of groupLabels.');

    % Unique group labels
    uniqueGroups = unique(groupLabels);

    % Number of groups and time points
    numGroups = numel(uniqueGroups);
    numTimePoints = size(dataMatrix, 2);

    % Preallocate arrays for efficiency
    groupMeans = zeros(numGroups, numTimePoints);
    groupSizes = zeros(numGroups, 1);

    % Compute group means and sizes
    for g = 1:numGroups
        groupIndices = groupLabels == uniqueGroups(g);
        groupMeans(g, :) = mean(dataMatrix(groupIndices, :), 1); % Mean for each group at each time point
        groupSizes(g) = sum(groupIndices); % Number of trials in each group
    end

    % Overall mean across all groups
    overallMean = mean(dataMatrix, 1);

    % Between-group variance (SSB)
    SSB = sum(groupSizes .* (groupMeans - overallMean).^2, 1);

    % Within-group variance (SSW)
    SSW = zeros(1, numTimePoints);
    for g = 1:numGroups
        groupIndices = groupLabels == uniqueGroups(g);
        SSW = SSW + sum((dataMatrix(groupIndices, :) - groupMeans(g, :)).^2, 1);
    end

    % Degrees of freedom
    dfBetween = numGroups - 1;
    dfWithin = size(dataMatrix, 1) - numGroups;

    % Mean squares
    MSB = SSB / dfBetween; % Mean square between groups
    MSW = SSW / dfWithin;  % Mean square within groups

    % F-statistic
    statMetric = MSB ./ MSW; % F-statistic for each time point
end
