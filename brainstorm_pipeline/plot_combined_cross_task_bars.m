function plot_combined_cross_task_bars(groupBarRows, targetTasks, outputDir, varargin)
% PLOT_COMBINED_CROSS_TASK_BARS  Single-axis HG bar plot across all tasks.
%
%   All target-task conditions are placed on one x-axis (e.g. vWM Hard,
%   vWM Easy, Math Hard, Math Easy, ...). Figures are saved under the
%   analysis output folder, never next to pipeline scripts.
%
%   plot_combined_cross_task_bars(groupBarRows, targetTasks, outputDir)
%   plot_combined_cross_task_bars   % standalone: reads latest CSV from output/
%
%   Name-Value:
%     sourceTask    - figure title label (default '')
%     channelGroup  - 'source_sig' (default) | 'source_nonsig' | 'both'
%     subjects      - {} group average | {'Subject12'} per subject
%     filePrefix    - output filename prefix (default 'combined')

if nargin < 1 || isempty(groupBarRows)
    cfg = default_combined_plot_settings();
    groupBarRows = load_group_bar_table(cfg.dataFile);
    targetTasks = cfg.targetTasks;
    outputDir = cfg.outputDir;
    varargin = struct_to_name_value(cfg.opts);
end

if nargin < 2 || isempty(targetTasks)
    targetTasks = unique(cellstr(groupBarRows.target_task), 'stable');
end
if nargin < 3 || isempty(outputDir)
    cfg = default_combined_plot_settings();
    outputDir = cfg.outputDir;
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

p = inputParser();
addParameter(p, 'sourceTask', '');
addParameter(p, 'channelGroup', 'source_sig');
addParameter(p, 'subjects', {});
addParameter(p, 'filePrefix', 'combined');
parse(p, varargin{:});
opts = p.Results;

targetTasks = subjects_to_cellstr(targetTasks);
subjects = subjects_to_cellstr(opts.subjects);

if isempty(subjects)
    plot_subject_or_group(groupBarRows, targetTasks, outputDir, opts, '');
else
    for s = 1:numel(subjects)
        subjOpts = opts;
        subjOpts.filePrefix = sprintf('%s_%s', opts.filePrefix, subjects{s});
        plot_subject_or_group(groupBarRows, targetTasks, outputDir, subjOpts, subjects{s});
    end
end
end


function cfg = default_combined_plot_settings()
repoRoot = fileparts(fileparts(mfilename('fullpath')));
cfg = struct();
cfg.workingDir = 'F:\seeg\luohong\analysisEV';
cfg.outputRoot = fullfile(cfg.workingDir, 'output', 'cross_task_sig_test');
cfg.outputDir = fullfile(cfg.outputRoot, ...
    'MITSWJNTask_to_vWM_Math_WM_MSIT_vMSIT_MITSWJNTask');
cfg.dataFile = fullfile(cfg.outputDir, 'condition_hg_bar_data.csv');

if ~isfile(cfg.dataFile)
    [cfg.outputDir, cfg.dataFile] = discover_latest_bar_csv(cfg.outputRoot);
end

cfg.targetTasks = {'vWM', 'Math', 'WM', 'MSIT', 'vMSIT', 'MITSWJNTask'};
cfg.opts = struct( ...
    'sourceTask', 'MITSWJNTask', ...
    'channelGroup', 'source_sig', ...
    'subjects', {{'Subject12'}}, ...
    'filePrefix', 'combined');
end


function [runDir, csvPath] = discover_latest_bar_csv(outputRoot)
runDir = outputRoot;
csvPath = fullfile(outputRoot, 'condition_hg_bar_data.csv');
if ~isfolder(outputRoot)
    error(['Cross-task output folder not found:\n  %s\n' ...
        'Run test_sig_channels_cross_task first.'], outputRoot);
end

d = dir(fullfile(outputRoot, '**', 'condition_hg_bar_data.csv'));
if isempty(d)
    error(['No condition_hg_bar_data.csv under:\n  %s\n' ...
        'Run test_sig_channels_cross_task first.'], outputRoot);
end

[~, ord] = sort([d.datenum], 'descend');
csvPath = fullfile(d(ord(1)).folder, d(ord(1)).name);
runDir = d(ord(1)).folder;
fprintf('Using latest bar data: %s\n', csvPath);
end


function T = load_group_bar_table(csvPath)
if ~isfile(csvPath)
    error('Bar data CSV not found:\n  %s\nRun test_sig_channels_cross_task first.', csvPath);
end
T = readtable(csvPath, 'TextType', 'string');
required = ["subject", "target_task", "condition", "channel_group", "mean_hg"];
missing = required(~ismember(required, string(T.Properties.VariableNames)));
if ~isempty(missing)
    error('CSV missing columns: %s', strjoin(missing, ', '));
end
if ~ismember("sem_hg", string(T.Properties.VariableNames))
    T.sem_hg = nan(height(T), 1);
end
if ~ismember("n_channels", string(T.Properties.VariableNames))
    T.n_channels = nan(height(T), 1);
end
end


function plot_subject_or_group(groupBarRows, targetTasks, outputDir, opts, subjFilter)
if strcmpi(opts.channelGroup, 'both')
    plot_single_axis_sig_vs_nonsig(groupBarRows, targetTasks, outputDir, opts, subjFilter);
else
    plot_single_axis_one_group(groupBarRows, targetTasks, outputDir, opts, subjFilter);
end
end


function plot_single_axis_one_group(groupBarRows, targetTasks, outputDir, opts, subjFilter)
flat = flatten_task_conditions(groupBarRows, targetTasks, subjFilter, opts.channelGroup);
if isempty(flat.xLabels)
    warning('No data for combined plot (%s).', opts.channelGroup);
    return;
end

nBars = numel(flat.xLabels);
figW = max(900, 90 * nBars + 160);
fig = figure('Color', 'w', 'Position', [40 40 figW 520], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

for j = 1:nBars
    bar(ax, j, flat.means(j), 0.72, 'FaceColor', flat.colors(j, :), ...
        'EdgeColor', 'k', 'LineWidth', 1.0);
    if isfinite(flat.sems(j)) && flat.sems(j) > 0
        errorbar(ax, j, flat.means(j), flat.sems(j), 'k', 'LineWidth', 1.2, 'CapSize', 8);
    end
end

decorate_single_axis(ax, flat);
title(ax, build_main_title(opts.sourceTask, opts.channelGroup, subjFilter, flat.meta), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

save_combined_figure(fig, outputDir, opts.filePrefix, opts.channelGroup);
close(fig);
fprintf('Saved combined figure: %s\n', fullfile(outputDir, ...
    combined_filename(opts.filePrefix, opts.channelGroup, '.png')));
end


function plot_single_axis_sig_vs_nonsig(groupBarRows, targetTasks, outputDir, opts, subjFilter)
flatSig = flatten_task_conditions(groupBarRows, targetTasks, subjFilter, 'source_sig');
flatNon = flatten_task_conditions(groupBarRows, targetTasks, subjFilter, 'source_nonsig');
if isempty(flatSig.xLabels)
    warning('No source-sig data for combined sig-vs-nonsig plot.');
    return;
end

nBars = numel(flatSig.xLabels);
figW = max(1000, 100 * nBars + 180);
fig = figure('Color', 'w', 'Position', [40 40 figW 540], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

sigColor = [214 39 40] / 255;
nonsigColor = [187 187 187] / 255;
groupW = 0.38;
offset = groupW / 2;

for j = 1:nBars
    if j == 1
        bar(ax, j - offset, flatSig.means(j), groupW, 'FaceColor', sigColor, ...
            'EdgeColor', 'k', 'DisplayName', 'Source sig');
    else
        bar(ax, j - offset, flatSig.means(j), groupW, 'FaceColor', sigColor, ...
            'EdgeColor', 'k', 'HandleVisibility', 'off');
    end
    if isfinite(flatSig.sems(j)) && flatSig.sems(j) > 0
        errorbar(ax, j - offset, flatSig.means(j), flatSig.sems(j), ...
            'k', 'LineWidth', 1.0, 'CapSize', 7);
    end

    if j <= numel(flatNon.means) && isfinite(flatNon.means(j))
        if j == 1
            bar(ax, j + offset, flatNon.means(j), groupW, 'FaceColor', nonsigColor, ...
                'EdgeColor', 'k', 'DisplayName', 'Source nonsig');
        else
            bar(ax, j + offset, flatNon.means(j), groupW, 'FaceColor', nonsigColor, ...
                'EdgeColor', 'k', 'HandleVisibility', 'off');
        end
        if isfinite(flatNon.sems(j)) && flatNon.sems(j) > 0
            errorbar(ax, j + offset, flatNon.means(j), flatNon.sems(j), ...
                'k', 'LineWidth', 1.0, 'CapSize', 7);
        end
    end
end

decorate_single_axis(ax, flatSig);
legend(ax, 'Location', 'best');
title(ax, build_main_title(opts.sourceTask, 'sig vs nonsig', subjFilter, flatSig.meta), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

prefix = opts.filePrefix;
if strcmpi(prefix, 'combined')
    prefix = 'combined_sig_vs_nonsig';
end
save_combined_figure(fig, outputDir, prefix, 'sig_vs_nonsig');
close(fig);
fprintf('Saved combined figure: %s\n', fullfile(outputDir, ...
    combined_filename(prefix, 'sig_vs_nonsig', '.png')));
end


function flat = flatten_task_conditions(groupBarRows, targetTasks, subjFilter, channelGroup)
flat = struct();
flat.xLabels = {};
flat.means = [];
flat.sems = [];
flat.colors = [];
flat.taskStarts = [];
flat.taskNames = {};
flat.meta = struct('n_subjects', 0, 'n_channels', 0);

x = 0;
for t = 1:numel(targetTasks)
    [conds, means, sems, meta] = aggregate_task_panel( ...
        groupBarRows, targetTasks{t}, subjFilter, channelGroup);
    if isempty(conds) || all(isnan(means))
        continue;
    end

    flat.taskStarts(end + 1) = x + 1; %#ok<AGROW>
    flat.taskNames{end + 1} = targetTasks{t}; %#ok<AGROW>
    flat.meta.n_subjects = max(flat.meta.n_subjects, meta.n_subjects);
    flat.meta.n_channels = max(flat.meta.n_channels, meta.n_channels);

    for c = 1:numel(conds)
        x = x + 1;
        flat.xLabels{end + 1} = short_condition_label(conds{c}); %#ok<AGROW>
        flat.means(end + 1) = means(c); %#ok<AGROW>
        flat.sems(end + 1) = sems(c); %#ok<AGROW>
        flat.colors(end + 1, :) = condition_type_color(conds{c}); %#ok<AGROW>
    end
end
end


function decorate_single_axis(ax, flat)
if isempty(flat.xLabels)
    return;
end

nBars = numel(flat.xLabels);
set(ax, 'XTick', 1:nBars, 'XTickLabel', flat.xLabels, 'FontSize', 9);
ylabel(ax, 'Z-scored high-gamma envelope (a.u.)');
grid(ax, 'on');
yline(ax, 0, 'k--');
set(ax, 'Box', 'off', 'LineWidth', 1.1);
xtickangle(ax, 45);

yl = ylim(ax);
yTop = yl(2);
yPad = 0.06 * max(diff(yl), 0.5);
yAnnot = yTop + yPad;

for k = 1:numel(flat.taskStarts)
    xStart = flat.taskStarts(k);
    if k < numel(flat.taskStarts)
        xEnd = flat.taskStarts(k + 1) - 1;
    else
        xEnd = nBars;
    end
    xMid = (xStart + xEnd) / 2;
    text(ax, xMid, yAnnot, flat.taskNames{k}, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
        'FontSize', 10, 'Clipping', 'off');
    if k > 1
        xline(ax, xStart - 0.5, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    end
end
ylim(ax, [yl(1), yAnnot + yPad]);
xlim(ax, [0.4, nBars + 0.6]);
end


function [conds, means, sems, meta] = aggregate_task_panel(groupBarRows, tgtTask, subjFilter, channelGroup)
conds = {};
means = [];
sems = [];
meta = struct('n_subjects', 0, 'n_channels', 0);

mask = groupBarRows.target_task == string(tgtTask) ...
    & groupBarRows.channel_group == string(channelGroup);
if ~isempty(subjFilter)
    mask = mask & groupBarRows.subject == string(subjFilter);
end
subTbl = groupBarRows(mask, :);
if isempty(subTbl)
    return;
end

conds = cellstr(unique(subTbl.condition, 'stable'));
means = nan(1, numel(conds));
sems = nan(1, numel(conds));

for c = 1:numel(conds)
    cm = subTbl.condition == string(conds{c});
    vals = subTbl.mean_hg(cm);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    means(c) = mean(vals, 'omitnan');
    if numel(vals) > 1
        sems(c) = std(vals, 0, 'omitnan') / sqrt(numel(vals));
    else
        rowSem = subTbl.sem_hg(cm);
        if any(isfinite(rowSem))
            sems(c) = mean(rowSem, 'omitnan');
        else
            sems(c) = 0;
        end
    end
end

meta.n_subjects = numel(unique(subTbl.subject));
meta.n_channels = max(subTbl.n_channels, [], 'omitnan');
if isnan(meta.n_channels)
    meta.n_channels = 0;
end
end


function titleStr = build_main_title(sourceTask, channelGroup, subjFilter, meta)
if ~isempty(subjFilter)
    subjPart = subjFilter;
elseif meta.n_subjects > 0
    subjPart = sprintf('%d subjects', meta.n_subjects);
else
    subjPart = 'group average';
end

if ~isempty(sourceTask)
    srcPart = sprintf('Source-sig from %s', sourceTask);
else
    srcPart = 'Cross-task condition HG';
end

titleStr = sprintf('%s | %s | %s (n=%d ch)', ...
    srcPart, strrep(channelGroup, '_', ' '), subjPart, meta.n_channels);
end


function save_combined_figure(fig, outputDir, filePrefix, channelGroup)
pngPath = fullfile(outputDir, combined_filename(filePrefix, channelGroup, '.png'));
pdfPath = fullfile(outputDir, combined_filename(filePrefix, channelGroup, '.pdf'));
exportgraphics(fig, pngPath, 'Resolution', 200);
try
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
catch
end
end


function fname = combined_filename(prefix, channelGroup, ext)
grpTag = matlab.lang.makeValidName(channelGroup);
fname = sprintf('%s_all_conditions_%s%s', prefix, grpTag, ext);
end


function rgb = condition_type_color(condName)
lab = lower(char(condName));
if contains(lab, 'hard') || contains(lab, 'sent')
    rgb = [214 39 40] / 255;
elseif contains(lab, 'easy') || contains(lab, 'nonword') || contains(lab, 'degrad')
    rgb = [31 119 180] / 255;
elseif contains(lab, 'jabber')
    rgb = [255 127 14] / 255;
elseif contains(lab, 'word')
    rgb = [44 160 44] / 255;
elseif contains(lab, 'intact')
    rgb = [214 39 40] / 255;
else
    rgb = [148 103 189] / 255;
end
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


function args = struct_to_name_value(s)
names = fieldnames(s);
args = cell(1, 2 * numel(names));
for i = 1:numel(names)
    args{2 * i - 1} = names{i};
    args{2 * i} = s.(names{i});
end
end
