function run_ifg_lang_response_analysis(varargin)
% RUN_IFG_LANG_RESPONSE_ANALYSIS  IFG language electrodes + HG plots.
%
%   Extracts language-responsive channels in IFG (HCP-MMP labels or MNI
%   fallback), builds HG word-average profiles, and saves figures.
%
%   Name-Value parameters:
%     analysisTitle   - banner string (default 'IFG LANGUAGE RESPONSES')
%     mniPlotTitle    - cohort MNI figure title prefix
%     taskType        - e.g. 'MITLangloc' or 'MITSWJNTask'
%     cohortName      - e.g. 'MGH' or 'local'
%     crunchedDir     - folder with finished sn_obj .mat files
%     cohortKind      - 'broadband' | 'brainstorm'
%     subjects        - {} auto-detect
%     outputDir       - required output folder
%     repoRoot        - repo root (default: parent of brainstorm_pipeline)
%     signalType      - 'bipolar' (default) | 'unipolar'
%     sigSource       - 's_vs_n' (default) | 'wordwise'
%     words           - default 1:12
%     S_condition     - default 'Sentences'
%     N_condition     - default '' (auto)
%     auxSearchDirs   - MNI enrichment search paths
%     doPerSubjectPlots, doCohortPlots, isPlotVisible, angle

p = inputParser();
addParameter(p, 'analysisTitle', 'IFG LANGUAGE RESPONSES');
addParameter(p, 'mniPlotTitle', 'IFG lang electrodes');
addParameter(p, 'taskType', 'MITLangloc');
addParameter(p, 'cohortName', 'MGH');
addParameter(p, 'crunchedDir', '');
addParameter(p, 'cohortKind', 'broadband');
addParameter(p, 'subjects', {});
addParameter(p, 'outputDir', '');
addParameter(p, 'repoRoot', '');
addParameter(p, 'signalType', 'bipolar');
addParameter(p, 'sigSource', 's_vs_n');
addParameter(p, 'words', 1:12);
addParameter(p, 'S_condition', 'Sentences');
addParameter(p, 'N_condition', '');
addParameter(p, 'auxSearchDirs', {});
addParameter(p, 'doPerSubjectPlots', true);
addParameter(p, 'doCohortPlots', true);
addParameter(p, 'isPlotVisible', false);
addParameter(p, 'angle', 270);
parse(p, varargin{:});
cfg = p.Results;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(cfg.repoRoot)
    cfg.repoRoot = fileparts(scriptDir);
end
if isempty(cfg.crunchedDir)
    error('crunchedDir is required.');
end
if isempty(cfg.outputDir)
    error('outputDir is required.');
end
if isempty(cfg.auxSearchDirs)
    cfg.auxSearchDirs = {fileparts(cfg.crunchedDir), cfg.crunchedDir, 'F:\iEEG_evlab'};
end

addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(cfg.repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

subjects = subjects_to_cellstr(cfg.subjects);
if isempty(subjects)
    subjects = discover_ifg_analysis_subjects(cfg.crunchedDir, cfg.taskType, cfg.cohortKind);
end
if isempty(subjects)
    error('No subjects with sn_obj found in %s (%s, %s).', ...
        cfg.crunchedDir, cfg.taskType, cfg.cohortKind);
end
if ~exist(cfg.outputDir, 'dir'); mkdir(cfg.outputDir); end

fprintf('\n=== %s ===\n', cfg.analysisTitle);
fprintf('Task:      %s (%s)\n', cfg.taskType, cfg.cohortKind);
fprintf('Cohort:    %s\n', cfg.cohortName);
fprintf('Subjects:  %s\n', strjoin(subjects, ', '));
fprintf('Signal:    %s (%s)\n', cfg.signalType, cfg.sigSource);
fprintf('Output:    %s\n', cfg.outputDir);

allProfiles = struct([]);
summaryRows = table();

for i = 1:numel(subjects)
    subj = subjects{i};
    fprintf('\n--- %s ---\n', subj);

    try
        [sn_obj, meta] = load_sn_obj_for_comparison(subj, cfg.crunchedDir, ...
            cfg.repoRoot, cfg.taskType, cfg.auxSearchDirs, cfg.cohortKind);
    catch ME
        warning('Skipping %s: %s', subj, ME.message);
        continue;
    end

    args = {'cohort', cfg.cohortName, 'subject', meta.subject, ...
        'signalType', cfg.signalType, 'sigSource', cfg.sigSource, ...
        'words', cfg.words, 'S_condition', cfg.S_condition, 'sigOnly', true, ...
        'use_odd_for_inference', false};
    if ~isempty(cfg.N_condition)
        args = [args, {'N_condition', cfg.N_condition}]; %#ok<AGROW>
    end

    try
        P = build_lang_channel_profiles(sn_obj, args{:});
    catch ME
        warning('Profile build failed for %s: %s', subj, ME.message);
        continue;
    end

    if isempty(P)
        fprintf('  No language-responsive channels with MNI.\n');
        continue;
    end

    ifgMask = false(numel(P), 1);
    for k = 1:numel(P)
        [ifgMask(k), infoK] = is_ifg_channel(sn_obj, P(k).chanIdx, cfg.signalType);
        P(k).ifg_method = infoK.method;
        if ~isempty(infoK.hcp_labels)
            P(k).hcp_label = strjoin(infoK.hcp_labels, '+');
        else
            P(k).hcp_label = '';
        end
    end

    Pifg = P(ifgMask);
    fprintf('  Lang sig w/ MNI: %d | IFG lang: %d | %s\n', ...
        numel(P), numel(Pifg), meta.file);

    if isempty(Pifg)
        continue;
    end

    allProfiles = append_profiles(allProfiles, Pifg);
    for k = 1:numel(Pifg)
        summaryRows = [summaryRows; profile_to_row(Pifg(k))]; %#ok<AGROW>
    end

    if cfg.doPerSubjectPlots
        subjDir = fullfile(cfg.outputDir, subj);
        if ~exist(subjDir, 'dir'); mkdir(subjDir); end
        plot_subject_ifg_responses(Pifg, subjDir, subj, cfg.signalType, cfg.isPlotVisible);
        plot_lang_sig_mni_distribution(sn_obj, subjDir, ...
            'subjectName', subj, 'signalType', cfg.signalType, ...
            'sigSource', cfg.sigSource, 'plotAllElectrodes', true, ...
            'isPlotVisible', cfg.isPlotVisible, 'angle', cfg.angle);
    end
end

if isempty(allProfiles)
    warning('No IFG language-responsive channels found across subjects.');
    return;
end

writetable(summaryRows, fullfile(cfg.outputDir, 'ifg_lang_profiles.csv'));
save(fullfile(cfg.outputDir, 'ifg_lang_profiles.mat'), ...
    'allProfiles', 'summaryRows', 'subjects', 'cfg', '-v7.3');

fprintf('\nCohort total: %d IFG lang channels across %d subjects\n', ...
    numel(allProfiles), numel(unique({allProfiles.subject})));

if cfg.doCohortPlots
    plot_cohort_ifg_responses(allProfiles, cfg.outputDir, cfg.signalType, cfg.isPlotVisible);
    plot_cohort_ifg_mni(summaryRows, subjects, cfg.crunchedDir, cfg.repoRoot, ...
        cfg.taskType, cfg.auxSearchDirs, cfg.cohortKind, cfg.mniPlotTitle, ...
        cfg.outputDir, cfg.angle, cfg.isPlotVisible);
end

fprintf('\nSaved IFG language response figures to:\n  %s\n', cfg.outputDir);
end


function plot_subject_ifg_responses(profiles, outputDir, subject, signalType, isPlotVisible)
nCh = numel(profiles);
nCol = min(4, nCh);
nRow = ceil(nCh / nCol);

fig = figure('Color', 'w', 'Position', [40 40 1600 900], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
tiled = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:nCh
    ax = nexttile(tiled);
    P = profiles(k);
    if isempty(P.word_ave_S) || isempty(P.word_ave_N)
        text(ax, 0.5, 0.5, 'No word-average data', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        continue;
    end

    t = profile_word_time_axis(P);
    plot(ax, t, P.word_ave_S, 'Color', [0 0.45 0.74], 'LineWidth', 1.4); hold(ax, 'on');
    plot(ax, t, P.word_ave_N, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2, 'LineStyle', '--');
    if ~isempty(P.word_ave_S_sem)
        shaded_error(ax, t, P.word_ave_S, P.word_ave_S_sem, [0 0.45 0.74]);
    end
    if ~isempty(P.word_ave_N_sem)
        shaded_error(ax, t, P.word_ave_N, P.word_ave_N_sem, [0.85 0.33 0.1]);
    end
    grid(ax, 'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'HG (z)');
    title(ax, sprintf('%s\n%s (%s)', P.label, P.hcp_label, P.ifg_method), ...
        'Interpreter', 'none', 'FontSize', 9);
    if k == 1
        legend(ax, {P.S_condition, P.N_condition}, 'Location', 'best', 'FontSize', 8);
    end
end

sgtitle(tiled, sprintf('%s - IFG lang electrodes (%s, n=%d)', subject, signalType, nCh));
pngName = fullfile(outputDir, sprintf('%s_ifg_lang_word_averages_%s.png', subject, signalType));
saveas(fig, pngName);
close(fig);
fprintf('  Saved %s\n', pngName);

fig2 = figure('Color', 'w', 'Position', [40 40 1600 900], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
tiled2 = tiledlayout(fig2, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nCh
    ax = nexttile(tiled2);
    P = profiles(k);
    if isempty(P.word_ave_S)
        continue;
    end
    t = profile_word_time_axis(P);
    diff = profile_s_minus_n(P);
    plot(ax, t, diff, 'k', 'LineWidth', 1.3); hold(ax, 'on');
    yline(ax, 0, ':');
    grid(ax, 'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'S - N (z)');
    title(ax, P.label, 'Interpreter', 'none', 'FontSize', 9);
end
sgtitle(tiled2, sprintf('%s - IFG S minus N (%s)', subject, signalType));
pngName2 = fullfile(outputDir, sprintf('%s_ifg_lang_S_minus_N_%s.png', subject, signalType));
saveas(fig2, pngName2);
close(fig2);
end


function plot_cohort_ifg_responses(profiles, outputDir, signalType, isPlotVisible)
[meanS, meanN, meanDiff, semS, semN, semDiff, tRef] = cohort_average_traces(profiles);

fig = figure('Color', 'w', 'Position', [100 100 1100 420], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
tiled = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled);
plot(ax1, tRef, meanS, 'Color', [0 0.45 0.74], 'LineWidth', 2); hold(ax1, 'on');
plot(ax1, tRef, meanN, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'LineStyle', '--');
shaded_error(ax1, tRef, meanS, semS, [0 0.45 0.74]);
shaded_error(ax1, tRef, meanN, semN, [0.85 0.33 0.1]);
grid(ax1, 'on');
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'HG (z)');
title(ax1, sprintf('IFG lang grand average (%s, n=%d)', signalType, numel(profiles)));
legend(ax1, {profiles(1).S_condition, profiles(1).N_condition}, 'Location', 'best');

ax2 = nexttile(tiled);
plot(ax2, tRef, meanDiff, 'k', 'LineWidth', 2); hold(ax2, 'on');
shaded_error(ax2, tRef, meanDiff, semDiff, [0.2 0.2 0.2]);
yline(ax2, 0, ':');
grid(ax2, 'on');
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'S - N (z)');
title(ax2, 'Grand average S - N');

pngName = fullfile(outputDir, sprintf('cohort_ifg_lang_grand_average_%s.png', signalType));
saveas(fig, pngName);
close(fig);
fprintf('Saved cohort average: %s\n', pngName);

fig3 = figure('Color', 'w', 'Position', [80 80 1400 900], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
nPlot = min(numel(profiles), 30);
nCol = 5;
nRow = ceil(nPlot / nCol);
tiled3 = tiledlayout(fig3, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nPlot
    ax = nexttile(tiled3);
    P = profiles(k);
    t = profile_word_time_axis(P);
    plot(ax, t, profile_s_minus_n(P), 'LineWidth', 1.1);
    grid(ax, 'on');
    title(ax, sprintf('%s %s', P.subject, P.label), 'Interpreter', 'none', 'FontSize', 8);
    if mod(k, nCol) == 1
        ylabel(ax, 'S-N');
    end
    if k > nPlot - nCol
        xlabel(ax, 'Time (s)');
    end
end
sgtitle(tiled3, sprintf('Individual IFG S-N traces (%s)', signalType));
pngName3 = fullfile(outputDir, sprintf('cohort_ifg_lang_individual_S_minus_N_%s.png', signalType));
saveas(fig3, pngName3);
close(fig3);
end


function plot_cohort_ifg_mni(T, subjects, crunchedDir, repoRoot, taskType, ...
    auxSearchDirs, cohortKind, plotTitle, outputDir, angle, isPlotVisible)

templateSn = [];
for i = 1:numel(subjects)
    try
        [sn, ~] = load_sn_obj_for_comparison(subjects{i}, crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, cohortKind);
        if isprop(sn, 'anatomy') && isstruct(sn.anatomy) ...
                && isfield(sn.anatomy, 'template_brain') ...
                && ~isempty(sn.anatomy.template_brain)
            templateSn = sn;
            break;
        end
    catch
    end
end

subjList = cellstr(unique(T.subject, 'stable'));
cmap = lines(numel(subjList));

fig = figure('Color', 'w', 'Position', [80 80 1300 900], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

if ~isempty(templateSn)
    try_plot_cortex_from_sn(ax, templateSn, 0.3);
    axis(ax, 'off');
else
    xlabel(ax, 'MNI X'); ylabel(ax, 'MNI Y'); zlabel(ax, 'MNI Z');
    grid(ax, 'on');
end

hLeg = gobjects(numel(subjList), 1);
for s = 1:numel(subjList)
    subj = subjList{s};
    mask = strcmp(string(T.subject), string(subj));
    coords = [T.mni_x(mask), T.mni_y(mask), T.mni_z(mask)];
    if isempty(coords)
        continue;
    end
    hLeg(s) = scatter3(ax, coords(:,1), coords(:,2), coords(:,3), ...
        110, cmap(s,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
end

axis(ax, 'equal');
view(ax, angle, 0);
if ~isempty(templateSn)
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
end
title(ax, sprintf('%s (n=%d)', plotTitle, height(T)));
legend(hLeg, subjList, 'Location', 'eastoutside', 'Interpreter', 'none');

pngName = fullfile(outputDir, 'cohort_ifg_lang_mni.png');
saveas(fig, pngName);
close(fig);
fprintf('Saved cohort MNI: %s\n', pngName);
end


function [meanS, meanN, meanDiff, semS, semN, semDiff, tRef] = cohort_average_traces(profiles)
nRef = 0;
for k = 1:numel(profiles)
    nRef = max(nRef, numel(profiles(k).word_ave_S));
end
if nRef < 1
    meanS = []; meanN = []; meanDiff = [];
    semS = []; semN = []; semDiff = [];
    tRef = [];
    return;
end

fs = profiles(1).sample_freq;
if isempty(fs) || fs <= 0
    fs = 1;
end
tRef = (0:nRef - 1) / fs;

matS = nan(numel(profiles), nRef);
matN = nan(numel(profiles), nRef);
for k = 1:numel(profiles)
    s = profiles(k).word_ave_S(:)';
    n = profiles(k).word_ave_N(:)';
    matS(k, 1:numel(s)) = s;
    matN(k, 1:numel(n)) = n;
end

meanS = mean(matS, 1, 'omitnan');
meanN = mean(matN, 1, 'omitnan');
meanDiff = meanS - meanN;
semS = std(matS, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(matS), 1));
semN = std(matN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(matN), 1));
semDiff = std(matS - matN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(matS - matN), 1));
end


function row = profile_to_row(P)
row = table(string(P.subject), string(P.cohort), string(P.label), ...
    P.chanIdx, P.mni(1), P.mni(2), P.mni(3), ...
    string(P.hcp_label), string(P.ifg_method), P.p_ratio, ...
    'VariableNames', {'subject', 'cohort', 'label', 'chan_idx', ...
    'mni_x', 'mni_y', 'mni_z', 'hcp_label', 'ifg_method', 'p_ratio'});
end


function profiles = append_profiles(profiles, P)
if isempty(P)
    return;
end
if isempty(profiles)
    profiles = P(:);
else
    profiles = [profiles(:); P(:)]; %#ok<AGROW>
end
end


function t = profile_word_time_axis(profile)
nSamp = numel(profile.word_ave_S);
fs = profile.sample_freq;
if isempty(fs) || fs <= 0
    fs = 1;
end
t = (0:nSamp - 1) / fs;
end


function diff = profile_s_minus_n(profile)
n = min(numel(profile.word_ave_S), numel(profile.word_ave_N));
diff = profile.word_ave_S(1:n) - profile.word_ave_N(1:n);
end


function shaded_error(ax, t, y, ysem, color)
if nargin < 5 || isempty(ysem) || all(isnan(ysem))
    return;
end
lo = y - ysem;
hi = y + ysem;
fill(ax, [t, fliplr(t)], [lo, fliplr(hi)], color, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none');
end


function tf = try_plot_cortex_from_sn(ax, sn_obj, alphaVal)
tf = false;
try
    tmpl = sn_obj.anatomy.template_brain;
    if exist('plot3DModel', 'file') == 2 && isfield(tmpl, 'cortex')
        plot3DModel(ax, tmpl.cortex, [], alphaVal);
        colormap(ax, gray);
        tf = true;
        return;
    end
    if isfield(tmpl, 'cortex')
        cortex = tmpl.cortex;
        if isfield(cortex, 'faces') && isfield(cortex, 'vertices')
            F = cortex.faces; V = cortex.vertices;
        elseif isfield(cortex, 'tri') && isfield(cortex, 'vert')
            F = cortex.tri; V = cortex.vert;
        else
            return;
        end
        patch(ax, 'Faces', F, 'Vertices', V, ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', alphaVal);
        tf = true;
    end
catch
end
end


function subjects = discover_ifg_analysis_subjects(crunchedDir, taskType, cohortKind)
switch lower(cohortKind)
    case 'brainstorm'
        patterns = {
            fullfile(crunchedDir, sprintf('*%s_crunched.mat', taskType))
            fullfile(crunchedDir, '*_MITLangloc_crunched.mat')
        };
    case 'broadband'
        patterns = {
            fullfile(crunchedDir, sprintf('*%s_crunched_HG_ZScore.mat', taskType))
            fullfile(crunchedDir, sprintf('*%s_crunched_HG.mat', taskType))
            fullfile(crunchedDir, '*_MITLangloc_crunched_HG_ZScore.mat')
            fullfile(crunchedDir, '*_MITLangloc_crunched_HG.mat')
        };
    otherwise
        error('Unknown cohortKind: %s', cohortKind);
end

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
    if ~has_sn_obj_in_file(matFile)
        continue;
    end
    subj = subject_from_analysis_filename(d(i).name, taskType, cohortKind);
    if ~isempty(subj)
        subjects{end+1} = subj; %#ok<AGROW>
    end
end
subjects = unique(subjects, 'stable');
end


function ok = has_sn_obj_in_file(matFile)
ok = false;
try
    vars = whos('-file', matFile);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function subj = subject_from_analysis_filename(fname, taskType, cohortKind)
subj = '';
switch lower(cohortKind)
    case 'brainstorm'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched\.mat$'], 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
            return;
        end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
        end
    case 'broadband'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched_HG(_ZScore)?\.mat$'], 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
            return;
        end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched_HG(_ZScore)?\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
        end
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
    error('Subject list must be char, string, or cell array.');
end
end


function v = ternary(cond, a, b)
if cond
    v = a;
else
    v = b;
end
end
