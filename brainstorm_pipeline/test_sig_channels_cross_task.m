function test_sig_channels_cross_task()
% TEST_SIG_CHANNELS_CROSS_TASK
%   Test whether channels significant in a source task respond to another
%   task's primary condition contrast (e.g. langloc S vs N -> WM Hard vs Easy).
%
%   For each subject with both tasks processed:
%     1) Load source-task sn_obj and read sig channels (s_vs_n or wordwise).
%     2) Load target-task sn_obj for the same subject.
%     3) Match channels by label (case/punctuation insensitive).
%     4) Compute Spearman rho and mean HG difference for the target contrast.
%     5) Bar-plot averaged z-scored high-gamma power per target condition
%        (grand average across source-sig matched channels).
%
%   Requires finished crunched files from complete_mit_pipeline_brainstorm:
%     <workingDir>/<Subject>_<TaskType>_crunched.mat
%
%   Saves CSV tables, group stats, per-task plots, and combined multi-task
%   figures (combined_all_conditions_*.png) under outputDir.

%% USER SETTINGS
repoRoot     = fileparts(fileparts(mfilename('fullpath')));
workingDir   = 'F:\seeg\luohong\analysisEV';

sourceTask   = 'MITSWJNTask';   % task used to define sig channels
targetTasks  = {'vWM','Math','WM','MSIT','vMSIT','MITSWJNTask'};  % tasks whose conditions are tested

subjects     = {'Subject12'};              % {} -> auto-discover from source task files
signalType   = 'bipolar';       % 'bipolar' | 'unipolar'
sigSource    = 's_vs_n';        % 's_vs_n' | 'wordwise'

% Override target contrast (leave '' to use get_mit_task_config defaults)
targetS_condition = struct();   % e.g. targetS_condition.WM = 'Hard'
targetN_condition = struct();   % e.g. targetN_condition.WM = 'Easy'

compareGroup = 'source_sig';    % 'source_sig' | 'all_valid'
doPlots      = true;
outputDir    = fullfile(workingDir, 'output', 'cross_task_sig_test', ...
    sprintf('%s_to_%s', sourceTask, strjoin(targetTasks, '_')));

auxSearchDirs = {workingDir, fullfile(workingDir, 'output', sourceTask)};

%% SETUP
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

sourceCfg = get_mit_task_config(sourceTask);
subjects = subjects_to_cellstr(subjects);
if isempty(subjects)
    subjects = discover_task_subjects(workingDir, sourceTask);
end
if isempty(subjects)
    error('No subjects with %s crunched files in %s', sourceTask, workingDir);
end

fprintf('\n=== CROSS-TASK SIG CHANNEL TEST ===\n');
fprintf('Source task:  %s (%s, %s)\n', sourceTask, signalType, sigSource);
fprintf('Target tasks: %s\n', strjoin(targetTasks, ', '));
fprintf('Subjects:     %s\n', strjoin(subjects, ', '));
fprintf('Output:       %s\n', outputDir);

allRows = table();
subjectSummaries = table();
statsOut = struct();
groupBarRows = table();

for s = 1:numel(subjects)
    subj = subjects{s};
    fprintf('\n--- %s ---\n', subj);

    try
        [snSrc, metaSrc] = load_sn_obj_for_comparison(subj, workingDir, ...
            repoRoot, sourceTask, auxSearchDirs, 'brainstorm');
    catch ME
        warning('Skipping %s (source): %s', subj, ME.message);
        continue;
    end

    try
        sigMaskSrc = get_lang_sig_mask(snSrc, signalType, sigSource);
    catch ME
        warning('Skipping %s: no sig mask on source (%s)', subj, ME.message);
        continue;
    end

    srcLabels = channel_labels(snSrc, signalType);
    srcContrast = compute_task_contrast(snSrc, sourceCfg, signalType, srcLabels);
    nSigSrc = sum(sigMaskSrc);

    fprintf('  Source %s: %d / %d sig channels | %s\n', ...
        sourceTask, nSigSrc, numel(sigMaskSrc), metaSrc.file);

    for t = 1:numel(targetTasks)
        tgtTask = targetTasks{t};
        tgtCfg = get_mit_task_config(tgtTask);
        [condA, condB] = resolve_target_conditions(tgtCfg, targetS_condition, targetN_condition, tgtTask);

        try
            [snTgt, metaTgt] = load_sn_obj_for_comparison(subj, workingDir, ...
                repoRoot, tgtTask, auxSearchDirs, 'brainstorm');
        catch ME
            warning('  %s: no target file (%s) — %s', subj, tgtTask, ME.message);
            continue;
        end

        tgtLabels = channel_labels(snTgt, signalType);
        [matchIdx, matchedMask] = match_channels_by_label(srcLabels, tgtLabels, signalType);

        tgtContrast = compute_task_contrast(snTgt, tgtCfg, signalType, tgtLabels, ...
            'condA', condA, 'condB', condB);

        nMatched = sum(matchedMask);
        nMatchedSig = sum(matchedMask & sigMaskSrc);
        fprintf('  -> %s: matched %d channels (%d source-sig) | contrast %s vs %s\n', ...
            tgtTask, nMatched, nMatchedSig, condA, condB);

        for ci = 1:numel(srcLabels)
            if ~matchedMask(ci)
                continue;
            end
            ti = matchIdx(ci);
            row = table();
            row.subject = string(subj);
            row.source_task = string(sourceTask);
            row.target_task = string(tgtTask);
            row.channel = string(srcLabels{ci});
            row.source_ch_idx = ci;
            row.target_ch_idx = ti;
            row.source_sig = logical(sigMaskSrc(ci));
            row.source_rho = srcContrast.rho(ci);
            row.source_mean_diff = srcContrast.mean_diff(ci);
            row.target_rho = tgtContrast.rho(ti);
            row.target_mean_diff = tgtContrast.mean_diff(ti);
            row.target_cond_A = string(condA);
            row.target_cond_B = string(condB);
            row.n_trials_A = tgtContrast.n_trials_A(ti);
            row.n_trials_B = tgtContrast.n_trials_B(ti);
            allRows = [allRows; row]; %#ok<AGROW>
        end

        subSum = summarize_subject_contrast(allRows, subj, tgtTask, sigMaskSrc, matchedMask);
        subjectSummaries = [subjectSummaries; subSum]; %#ok<AGROW>

        key = matlab.lang.makeValidName(sprintf('%s_%s', subj, tgtTask));
        statsOut.(key) = group_compare_channels(allRows, subj, tgtTask, compareGroup);

        tgtConds = available_task_conditions(snTgt, tgtCfg);
        sigChanTgt = matchIdx(matchedMask & sigMaskSrc);
        nonsigChanTgt = matchIdx(matchedMask & ~sigMaskSrc);
        [hgSig, hgSigSem] = mean_hg_by_condition(snTgt, tgtCfg, signalType, sigChanTgt, tgtConds);
        [hgNon, hgNonSem] = mean_hg_by_condition(snTgt, tgtCfg, signalType, nonsigChanTgt, tgtConds);

        if doPlots
            subjPlotDir = fullfile(outputDir, subj);
            if ~exist(subjPlotDir, 'dir')
                mkdir(subjPlotDir);
            end
            plot_condition_hg_bars(tgtConds, hgSig, hgSigSem, ...
                sprintf('%s | %s\nsource-sig channels (n=%d)', subj, tgtTask, numel(sigChanTgt)), ...
                fullfile(subjPlotDir, sprintf('%s_%s_source_sig_hg_bar.png', subj, tgtTask)));
            if ~isempty(nonsigChanTgt)
                plot_condition_hg_bars(tgtConds, hgNon, hgNonSem, ...
                    sprintf('%s | %s\nsource-nonsig channels (n=%d)', subj, tgtTask, numel(nonsigChanTgt)), ...
                    fullfile(subjPlotDir, sprintf('%s_%s_source_nonsig_hg_bar.png', subj, tgtTask)));
            end
            plot_grouped_condition_bars(tgtConds, hgSig, hgSigSem, hgNon, hgNonSem, ...
                sprintf('%s | %s\nsource-sig vs source-nonsig', subj, tgtTask), ...
                fullfile(subjPlotDir, sprintf('%s_%s_sig_vs_nonsig_hg_bar.png', subj, tgtTask)));
        end

        for c = 1:numel(tgtConds)
            groupBarRows = [groupBarRows; table( ...
                string(subj), string(tgtTask), string(tgtConds{c}), ...
                "source_sig", numel(sigChanTgt), hgSig(c), hgSigSem(c), ...
                'VariableNames', {'subject', 'target_task', 'condition', ...
                'channel_group', 'n_channels', 'mean_hg', 'sem_hg'})]; %#ok<AGROW>
            if ~isempty(nonsigChanTgt)
                groupBarRows = [groupBarRows; table( ...
                    string(subj), string(tgtTask), string(tgtConds{c}), ...
                    "source_nonsig", numel(nonsigChanTgt), hgNon(c), hgNonSem(c), ...
                    'VariableNames', {'subject', 'target_task', 'condition', ...
                    'channel_group', 'n_channels', 'mean_hg', 'sem_hg'})]; %#ok<AGROW>
            end
        end
    end
end

if isempty(allRows)
    error('No cross-task channel rows collected. Check crunched files and subject list.');
end

writetable(allRows, fullfile(outputDir, 'channel_cross_task_metrics.csv'));
if ~isempty(subjectSummaries)
    writetable(subjectSummaries, fullfile(outputDir, 'subject_summary.csv'));
end
if ~isempty(groupBarRows)
    writetable(groupBarRows, fullfile(outputDir, 'condition_hg_bar_data.csv'));
end

groupStats = struct();
for t = 1:numel(targetTasks)
    tgtTask = targetTasks{t};
    groupStats.(matlab.lang.makeValidName(tgtTask)) = ...
        group_compare_channels(allRows, '', tgtTask, compareGroup);
end
save(fullfile(outputDir, 'cross_task_stats.mat'), 'statsOut', 'groupStats', 'allRows', ...
    'groupBarRows', 'sourceTask', 'targetTasks', 'signalType', 'sigSource', 'compareGroup');

print_group_summary(groupStats, targetTasks);

if doPlots && ~isempty(groupBarRows)
    plot_group_average_bars(groupBarRows, targetTasks, outputDir);
    plot_combined_cross_task_bars(groupBarRows, targetTasks, outputDir, ...
        'sourceTask', sourceTask, 'channelGroup', 'source_sig', ...
        'subjects', subjects, 'filePrefix', 'combined');
    plot_combined_cross_task_bars(groupBarRows, targetTasks, outputDir, ...
        'sourceTask', sourceTask, 'channelGroup', 'both', ...
        'subjects', subjects, 'filePrefix', 'combined');
end

fprintf('\nSaved results to: %s\n', outputDir);
end


function labels = channel_labels(sn_obj, signalType)
if strcmpi(signalType, 'bipolar')
    labels = sn_obj.bip_ch_label(:);
else
    labels = sn_obj.elec_ch_label(:);
end
end


function [condA, condB] = resolve_target_conditions(tgtCfg, targetS_condition, targetN_condition, tgtTask)
condA = tgtCfg.S_condition;
condB = tgtCfg.N_condition;
if isstruct(targetS_condition) && isfield(targetS_condition, tgtTask) ...
        && ~isempty(targetS_condition.(tgtTask))
    condA = targetS_condition.(tgtTask);
end
if isstruct(targetN_condition) && isfield(targetN_condition, tgtTask) ...
        && ~isempty(targetN_condition.(tgtTask))
    condB = targetN_condition.(tgtTask);
end
end


function contrast = compute_task_contrast(sn_obj, taskCfg, signalType, labels, varargin)
p = inputParser();
addParameter(p, 'condA', taskCfg.S_condition);
addParameter(p, 'condB', taskCfg.N_condition);
parse(p, varargin{:});
condA = p.Results.condA;
condB = p.Results.condB;

nChan = numel(labels);
contrast = struct();
contrast.rho = nan(nChan, 1);
contrast.mean_diff = nan(nChan, 1);
contrast.n_trials_A = nan(nChan, 1);
contrast.n_trials_B = nan(nChan, 1);
contrast.cond_A = condA;
contrast.cond_B = condB;

avail = unique(sn_obj.condition, 'stable');
condA = resolve_condition_name(avail, {condA}, 'A');
condB = resolve_condition_name(avail, {condB}, 'B');
words = taskCfg.testWords;
words = clamp_words(sn_obj, words, condA);

dataField = contrast_data_field(signalType);
try
    [A_tbl, ~] = sn_obj.get_ave_cond_trial('words', words, 'condition', condA);
    [B_tbl, ~] = sn_obj.get_ave_cond_trial('words', words, 'condition', condB);
catch ME
    warning('compute_task_contrast:Data', '%s: %s', get_subject_id(sn_obj), ME.message);
    return;
end

if ~ismember(dataField, A_tbl.Properties.VariableNames) ...
        || isempty(A_tbl.(dataField){1})
    return;
end

A_data = A_tbl.(dataField){1};
B_data = B_tbl.(dataField){1};
nElec = min(size(A_data, 1), nChan);

for e = 1:nElec
    [rho, md, nA, nB] = channel_contrast_metrics(A_data(e, :, :), B_data(e, :, :));
    contrast.rho(e) = rho;
    contrast.mean_diff(e) = md;
    contrast.n_trials_A(e) = nA;
    contrast.n_trials_B(e) = nB;
end
end


function [rho, meanDiff, nA, nB] = channel_contrast_metrics(A3d, B3d)
A_means = squeeze(mean(A3d, 3));
B_means = squeeze(mean(B3d, 3));
if isrow(A_means); A_means = A_means(:); end
if isrow(B_means); B_means = B_means(:); end

A_means = A_means(isfinite(A_means));
B_means = B_means(isfinite(B_means));
nA = numel(A_means);
nB = numel(B_means);

if nA < 2 || nB < 2
    rho = NaN;
    meanDiff = NaN;
    return;
end

y = [A_means; B_means];
x = [ones(nA, 1); -ones(nB, 1)];
rho = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
if ~isscalar(rho)
    rho = rho(1);
end
meanDiff = mean(A_means, 'omitnan') - mean(B_means, 'omitnan');
end


function fieldName = contrast_data_field(signalType)
if strcmpi(signalType, 'bipolar')
    fieldName = 'bip_elec_data';
else
    fieldName = 'elec_data';
end
end


function words = clamp_words(sn_obj, words, condition)
maxW = 0;
for w = 1:24
    try
        sn_obj.get_ave_cond_trial('words', w, 'condition', condition);
        maxW = w;
    catch
        break;
    end
end
if maxW < 1
    maxW = max(numel(words), 1);
end
words = words(words >= 1 & words <= maxW);
if isempty(words)
    words = 1:maxW;
end
end


function [matchIdx, matchedMask] = match_channels_by_label(srcLabels, tgtLabels, signalType)
nSrc = numel(srcLabels);
matchIdx = nan(nSrc, 1);
matchedMask = false(nSrc, 1);

srcNorm = normalize_labels(srcLabels);
tgtNorm = normalize_labels(tgtLabels);

for i = 1:nSrc
    hit = find(strcmp(srcNorm(i), tgtNorm), 1, 'first');
    if ~isempty(hit)
        matchIdx(i) = hit;
        matchedMask(i) = true;
    end
end

if strcmpi(signalType, 'bipolar')
    return;
end

% Unipolar fallback: allow partial match when exact label missing
for i = 1:nSrc
    if matchedMask(i)
        continue;
    end
    for j = 1:numel(tgtNorm)
        if strlength(srcNorm(i)) >= 3 && contains(tgtNorm(j), srcNorm(i))
            matchIdx(i) = j;
            matchedMask(i) = true;
            break;
        end
    end
end
end


function out = normalize_labels(labels)
out = upper(string(labels));
out = regexprep(out, '[^A-Z0-9]', '');
end


function subSum = summarize_subject_contrast(allRows, subj, tgtTask, sigMaskSrc, matchedMask)
subSum = table();
mask = allRows.subject == string(subj) & allRows.target_task == string(tgtTask);
if ~any(mask)
    return;
end
subTbl = allRows(mask, :);

sigVals = subTbl.target_rho(subTbl.source_sig);
nonsigVals = subTbl.target_rho(~subTbl.source_sig);

subSum = table( ...
    string(subj), string(tgtTask), sum(sigMaskSrc), sum(matchedMask), ...
    sum(subTbl.source_sig), numel(sigVals), numel(nonsigVals), ...
    median(sigVals, 'omitnan'), median(nonsigVals, 'omitnan'), ...
    mean(sigVals, 'omitnan'), mean(nonsigVals, 'omitnan'), ...
    'VariableNames', {'subject', 'target_task', 'n_source_sig', ...
    'n_matched_channels', 'n_matched_source_sig', 'n_sig_rows', 'n_nonsig_rows', ...
    'median_target_rho_sig', 'median_target_rho_nonsig', ...
    'mean_target_rho_sig', 'mean_target_rho_nonsig'});
end


function res = group_compare_channels(allRows, subjFilter, tgtTask, compareGroup)
mask = allRows.target_task == string(tgtTask) & isfinite(allRows.target_rho);
if ~isempty(subjFilter)
    mask = mask & allRows.subject == string(subjFilter);
end
subTbl = allRows(mask, :);

res = struct();
res.target_task = char(tgtTask);
if ~isempty(subjFilter)
    res.subject = char(subjFilter);
else
    res.subject = 'all';
end
res.n_channels = height(subTbl);
res.n_source_sig = sum(subTbl.source_sig);
res.n_source_nonsig = sum(~subTbl.source_sig);

sigVals = subTbl.target_rho(subTbl.source_sig);
nonsigVals = subTbl.target_rho(~subTbl.source_sig);
res.median_rho_sig = median(sigVals, 'omitnan');
res.median_rho_nonsig = median(nonsigVals, 'omitnan');
res.mean_rho_sig = mean(sigVals, 'omitnan');
res.mean_rho_nonsig = mean(nonsigVals, 'omitnan');

mdSig = subTbl.target_mean_diff(subTbl.source_sig);
mdNon = subTbl.target_mean_diff(~subTbl.source_sig);
res.median_mean_diff_sig = median(mdSig, 'omitnan');
res.median_mean_diff_nonsig = median(mdNon, 'omitnan');

if numel(sigVals) >= 2 && numel(nonsigVals) >= 2
    [pTwo, ~, st] = ranksum(sigVals, nonsigVals);
    [pRight, ~, ~] = ranksum(sigVals, nonsigVals, 'Tail', 'right');
    res.ranksum_p_two_sided = pTwo;
    res.ranksum_p_sig_greater = pRight;
    res.ranksum_z = st.zval;
    res.cliffs_delta = cliffs_delta(nonsigVals, sigVals);
else
    res.ranksum_p_two_sided = NaN;
    res.ranksum_p_sig_greater = NaN;
    res.ranksum_z = NaN;
    res.cliffs_delta = NaN;
end

if strcmpi(compareGroup, 'all_valid')
    res.compare_note = 'source_sig vs source_nonsig among matched channels';
else
    res.compare_note = 'source_sig vs source_nonsig among matched channels';
end
end


function print_group_summary(groupStats, targetTasks)
fprintf('\n=== GROUP SUMMARY (target rho: source-sig vs source-nonsig) ===\n');
for t = 1:numel(targetTasks)
    tgtTask = targetTasks{t};
    key = matlab.lang.makeValidName(tgtTask);
    if ~isfield(groupStats, key)
        continue;
    end
    G = groupStats.(key);
    fprintf('\nTarget: %s\n', tgtTask);
    fprintf('  Channels: %d (sig=%d, nonsig=%d)\n', G.n_channels, G.n_source_sig, G.n_source_nonsig);
    fprintf('  Median rho  sig=%.3f  nonsig=%.3f\n', G.median_rho_sig, G.median_rho_nonsig);
    if isfinite(G.ranksum_p_two_sided)
        fprintf('  Mann-Whitney (sig > nonsig): p=%.4g  Cliff''s d=%.3f\n', ...
            G.ranksum_p_sig_greater, G.cliffs_delta);
    end
end
end


function conds = available_task_conditions(sn_obj, taskCfg)
candidates = {taskCfg.S_condition, taskCfg.N_condition, ...
    taskCfg.W_condition, taskCfg.J_condition};
avail = unique(sn_obj.condition, 'stable');
conds = {};
for i = 1:numel(candidates)
    if isempty(candidates{i})
        continue;
    end
    try
        conds{end + 1} = resolve_condition_name(avail, {candidates{i}}, 'cond'); %#ok<AGROW>
    catch
    end
end
if isempty(conds)
    error('No task conditions resolved for %s.', get_subject_id(sn_obj));
end
conds = unique(conds, 'stable');
end


function [groupMean, groupSem] = mean_hg_by_condition(sn_obj, taskCfg, signalType, chanIdx, conditions)
nCond = numel(conditions);
groupMean = nan(1, nCond);
groupSem = nan(1, nCond);
if isempty(chanIdx)
    return;
end

chanIdx = chanIdx(:)';
chanIdx = chanIdx(isfinite(chanIdx) & chanIdx >= 1);
if isempty(chanIdx)
    return;
end

dataField = contrast_data_field(signalType);
words = clamp_words(sn_obj, taskCfg.testWords, conditions{1});
chanScalars = nan(numel(chanIdx), nCond);

for c = 1:nCond
    try
        [tbl, ~] = sn_obj.get_ave_cond_trial('words', words, 'condition', conditions{c});
    catch
        continue;
    end
    if ~ismember(dataField, tbl.Properties.VariableNames) || isempty(tbl.(dataField){1})
        continue;
    end
    data3d = tbl.(dataField){1};
    for k = 1:numel(chanIdx)
        ci = chanIdx(k);
        if ci <= size(data3d, 1)
            chanScalars(k, c) = mean(data3d(ci, :, :), 'all', 'omitnan');
        end
    end
end

for c = 1:nCond
    vals = chanScalars(:, c);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    groupMean(c) = mean(vals, 'omitnan');
    if numel(vals) > 1
        groupSem(c) = std(vals, 0, 'omitnan') / sqrt(numel(vals));
    else
        groupSem(c) = 0;
    end
end
end


function plot_condition_hg_bars(conditions, means, sems, titleStr, outPath)
if isempty(conditions) || all(isnan(means))
    return;
end

labels = cellfun(@short_condition_label, conditions, 'UniformOutput', false);
colors = condition_bar_colors(numel(conditions));

fig = figure('Color', 'w', 'Position', [100 100 560 480], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

x = 1:numel(conditions);
for j = 1:numel(conditions)
    bar(ax, x(j), means(j), 0.65, 'FaceColor', colors(j, :), ...
        'EdgeColor', 'k', 'LineWidth', 1.1);
    if isfinite(sems(j)) && sems(j) > 0
        errorbar(ax, x(j), means(j), sems(j), 'k', 'LineWidth', 1.4, 'CapSize', 10);
    end
end

set(ax, 'XTick', x, 'XTickLabel', labels, 'FontSize', 10);
ylabel(ax, 'Z-scored high-gamma envelope (a.u.)');
title(ax, titleStr, 'Interpreter', 'none');
grid(ax, 'on');
yline(ax, 0, 'k--');
set(ax, 'Box', 'off', 'LineWidth', 1.2);

saveas(fig, outPath);
close(fig);
end


function plot_grouped_condition_bars(conditions, sigMean, sigSem, nonsigMean, nonsigSem, titleStr, outPath)
if isempty(conditions) || all(isnan(sigMean))
    return;
end

labels = cellfun(@short_condition_label, conditions, 'UniformOutput', false);
sigColor = [214 39 40] / 255;
nonsigColor = [187 187 187] / 255;
nCond = numel(conditions);

fig = figure('Color', 'w', 'Position', [100 100 640 480], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

x = 1:nCond;
groupW = 0.36;
offset = groupW / 2;

for j = 1:nCond
    if j == 1
        bar(ax, x(j) - offset, sigMean(j), groupW, 'FaceColor', sigColor, ...
            'EdgeColor', 'k', 'DisplayName', 'Source sig');
    else
        bar(ax, x(j) - offset, sigMean(j), groupW, 'FaceColor', sigColor, ...
            'EdgeColor', 'k', 'HandleVisibility', 'off');
    end
    if isfinite(sigSem(j)) && sigSem(j) > 0
        errorbar(ax, x(j) - offset, sigMean(j), sigSem(j), 'k', 'LineWidth', 1.2, 'CapSize', 8);
    end
    if ~all(isnan(nonsigMean))
        if j == 1
            bar(ax, x(j) + offset, nonsigMean(j), groupW, 'FaceColor', nonsigColor, ...
                'EdgeColor', 'k', 'DisplayName', 'Source nonsig');
        else
            bar(ax, x(j) + offset, nonsigMean(j), groupW, 'FaceColor', nonsigColor, ...
                'EdgeColor', 'k', 'HandleVisibility', 'off');
        end
        if isfinite(nonsigSem(j)) && nonsigSem(j) > 0
            errorbar(ax, x(j) + offset, nonsigMean(j), nonsigSem(j), 'k', 'LineWidth', 1.2, 'CapSize', 8);
        end
    end
end

set(ax, 'XTick', x, 'XTickLabel', labels, 'FontSize', 10);
ylabel(ax, 'Z-scored high-gamma envelope (a.u.)');
title(ax, titleStr, 'Interpreter', 'none');
legend(ax, {'Source sig', 'Source nonsig'}, 'Location', 'best');
grid(ax, 'on');
yline(ax, 0, 'k--');
set(ax, 'Box', 'off', 'LineWidth', 1.2);

saveas(fig, outPath);
close(fig);
end


function plot_group_average_bars(groupBarRows, targetTasks, outputDir)
for t = 1:numel(targetTasks)
    tgtTask = targetTasks{t};
    taskMask = groupBarRows.target_task == string(tgtTask);
    if ~any(taskMask)
        continue;
    end

    for grp = ["source_sig", "source_nonsig"]
        grpMask = taskMask & groupBarRows.channel_group == grp;
        if ~any(grpMask)
            continue;
        end
        subTbl = groupBarRows(grpMask, :);
        conds = unique(subTbl.condition, 'stable');
        means = nan(1, numel(conds));
        sems = nan(1, numel(conds));
        nCh = 0;

        for c = 1:numel(conds)
            cm = subTbl.condition == conds(c);
            vals = subTbl.mean_hg(cm);
            vals = vals(isfinite(vals));
            if isempty(vals)
                continue;
            end
            means(c) = mean(vals, 'omitnan');
            if numel(vals) > 1
                sems(c) = std(vals, 0, 'omitnan') / sqrt(numel(vals));
            else
                sems(c) = 0;
            end
            nCh = max(nCh, max(subTbl.n_channels(cm)));
        end

        titleStr = sprintf('Group average | %s\n%s channels (max n=%d per subject)', ...
            tgtTask, strrep(char(grp), '_', ' '), nCh);
        fname = sprintf('group_%s_%s_hg_bar.png', matlab.lang.makeValidName(tgtTask), grp);
        plot_condition_hg_bars(cellstr(conds), means, sems, titleStr, fullfile(outputDir, fname));
    end

    plot_group_average_sig_vs_nonsig(groupBarRows, tgtTask, outputDir);
end
end


function plot_group_average_sig_vs_nonsig(groupBarRows, tgtTask, outputDir)
taskMask = groupBarRows.target_task == string(tgtTask);
sigTbl = groupBarRows(taskMask & groupBarRows.channel_group == "source_sig", :);
nonTbl = groupBarRows(taskMask & groupBarRows.channel_group == "source_nonsig", :);
if isempty(sigTbl)
    return;
end

conds = unique(sigTbl.condition, 'stable');
sigMean = nan(1, numel(conds));
sigSem = nan(1, numel(conds));
nonMean = nan(1, numel(conds));
nonSem = nan(1, numel(conds));

for c = 1:numel(conds)
    sigVals = sigTbl.mean_hg(sigTbl.condition == conds(c));
    sigVals = sigVals(isfinite(sigVals));
    if ~isempty(sigVals)
        sigMean(c) = mean(sigVals, 'omitnan');
        if numel(sigVals) > 1
            sigSem(c) = std(sigVals, 0, 'omitnan') / sqrt(numel(sigVals));
        else
            sigSem(c) = 0;
        end
    end
    if ~isempty(nonTbl)
        nonVals = nonTbl.mean_hg(nonTbl.condition == conds(c));
        nonVals = nonVals(isfinite(nonVals));
        if ~isempty(nonVals)
            nonMean(c) = mean(nonVals, 'omitnan');
            if numel(nonVals) > 1
                nonSem(c) = std(nonVals, 0, 'omitnan') / sqrt(numel(nonVals));
            else
                nonSem(c) = 0;
            end
        end
    end
end

titleStr = sprintf('Group average | %s\nsource-sig vs source-nonsig', tgtTask);
fname = sprintf('group_%s_sig_vs_nonsig_hg_bar.png', matlab.lang.makeValidName(tgtTask));
plot_grouped_condition_bars(cellstr(conds), sigMean, sigSem, nonMean, nonSem, ...
    titleStr, fullfile(outputDir, fname));
end


function lab = short_condition_label(condName)
lab = char(condName);
if contains(lower(lab), 'sent')
    lab = 'Sentences';
elseif contains(lower(lab), 'nonword')
    lab = 'Nonwords';
elseif contains(lower(lab), 'jabber')
    lab = 'Jabberwocky';
elseif contains(lower(lab), 'word') && ~contains(lower(lab), 'non')
    lab = 'Word lists';
elseif strcmpi(lab, 'Hard')
    lab = 'Hard';
elseif strcmpi(lab, 'Easy')
    lab = 'Easy';
elseif contains(lower(lab), 'intact')
    lab = 'Intact';
elseif contains(lower(lab), 'degrad')
    lab = 'Degraded';
end
end


function colors = condition_bar_colors(nCond)
palette = [
    214 39 40
    31 119 180
    44 160 44
    255 127 14
    148 103 189
    140 86 75] / 255;
colors = palette(mod(0:nCond - 1, size(palette, 1)) + 1, :);
end


function d = cliffs_delta(x, y)
nx = numel(x);
ny = numel(y);
if nx == 0 || ny == 0
    d = NaN;
    return;
end
greater = 0;
less = 0;
for i = 1:nx
    greater = greater + sum(x(i) > y);
    less = less + sum(x(i) < y);
end
d = (greater - less) / (nx * ny);
end


function id = get_subject_id(sn_obj)
id = 'unknown';
if isprop(sn_obj, 'subject') && ~isempty(sn_obj.subject)
    id = sn_obj.subject;
end
end


function out = subjects_to_cellstr(subjects)
if ischar(subjects)
    out = {subjects};
elseif isstring(subjects)
    out = cellstr(subjects(:))';
elseif iscell(subjects)
    out = cellfun(@char, subjects, 'UniformOutput', false);
    out = out(:)';
else
    out = {};
end
end


function subjects = discover_task_subjects(crunchedDir, taskType)
if ~isfolder(crunchedDir)
    subjects = {};
    return;
end

patterns = {
    fullfile(crunchedDir, sprintf('*%s_crunched.mat', taskType))
    fullfile(crunchedDir, '*_MITLangloc_crunched.mat')
};

d = [];
for p = 1:numel(patterns)
    d = [d; dir(patterns{p})]; %#ok<AGROW>
end
if isempty(d)
    subjects = {};
    return;
end

[~, ia] = unique({d.name}, 'stable');
d = d(ia);

subjects = {};
for i = 1:numel(d)
    matFile = fullfile(d(i).folder, d(i).name);
    if ~has_sn_obj_file(matFile)
        continue;
    end
    subj = subject_from_crunched_name(d(i).name, taskType);
    if ~isempty(subj)
        subjects{end + 1} = subj; %#ok<AGROW>
    end
end
subjects = unique(subjects, 'stable');
end


function ok = has_sn_obj_file(fname)
ok = false;
try
    vars = whos('-file', fname);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function subj = subject_from_crunched_name(fname, taskType)
subj = '';
tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched\.mat$'], 'tokens', 'once');
if ~isempty(tok)
    subj = tok{1};
    return;
end
tok = regexp(fname, '^(.+)_MITLangloc_crunched\.mat$', 'tokens', 'once');
if ~isempty(tok)
    subj = tok{1};
end
end
