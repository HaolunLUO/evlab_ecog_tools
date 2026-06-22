function compare_top_lang_electrodes()
% COMPARE_TOP_LANG_ELECTRODES  Top-10% lang electrodes: MITLangloc vs MITSWJNTask.
%
%   Selects the best language-responsive bipolar channels (top 10% among
%   s_vs_n significant electrodes) ranked by S-vs-N correlation and
%   split-half reliability, then compares cohort-level HG word-average
%   profiles between MITSWJNTask (local) and MITLangloc (MGH).
%
%   Local:  <workingDir>/<Subject>_MITSWJNTask_crunched.mat
%   Other:  <workingDir>/crunched/MITLangloc/<Subject>_MITLangloc_crunched_HG_ZScore.mat
%
%   Run:
%     compare_top_lang_electrodes

%% USER SETTINGS
repoRoot   = fileparts(fileparts(mfilename('fullpath')));
workingDir = 'F:\seeg\luohong\analysisEV';

localTaskType  = 'MITSWJNTask';
otherTaskType  = 'MITLangloc';

localCohortName  = 'local';
localCrunchedDir = workingDir;
otherCohortName  = 'MGH';
otherCrunchedDir = fullfile(workingDir, 'crunched', 'MITLangloc');

localSubjects = {
    'Subject04', ...
    'Subject06', ...
    'Subject07', ...
};
otherSubjects = {
    'AMC091', ...
    'AMC092', ...
    'AMC099', ...
};

topPct        = 0.10;
signalType    = 'bipolar';
sigSource     = 's_vs_n';
S_condition   = 'Sentences';
N_condition   = '';

outputDir = fullfile(workingDir, 'output', 'cross_task_top_lang', ...
    sprintf('top%d_pct_%s', round(topPct * 100), lower(signalType)));
doPlots   = true;
auxSearchDirs = {workingDir, otherCrunchedDir, 'F:\iEEG_evlab'};

%% SETUP
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

if isempty(localSubjects)
    localSubjects = discover_cohort_subjects(localCrunchedDir, localTaskType, 'brainstorm');
end
if isempty(otherSubjects)
    otherSubjects = discover_cohort_subjects(otherCrunchedDir, otherTaskType, 'broadband');
end
localSubjects = subjects_to_cellstr(localSubjects);
otherSubjects = subjects_to_cellstr(otherSubjects);

if ~exist(outputDir, 'dir'); mkdir(outputDir); end

fprintf('\n=== TOP-%d%% LANG ELECTRODE COMPARISON ===\n', round(topPct * 100));
fprintf('Ranking: S-vs-N correlation + split-half reliability\n');
fprintf('Local task:  %s (n=%d subjects)\n', localTaskType, numel(localSubjects));
fprintf('Other task:  %s (n=%d subjects)\n', otherTaskType, numel(otherSubjects));
fprintf('Top pct:     %.0f%%\n', topPct * 100);
fprintf('Output:      %s\n', outputDir);

[profilesLocal, rankLocal, summaryLocal] = load_top_profiles( ...
    localSubjects, localCrunchedDir, localCohortName, repoRoot, localTaskType, ...
    signalType, sigSource, S_condition, N_condition, topPct, auxSearchDirs, 'brainstorm');

[profilesOther, rankOther, summaryOther] = load_top_profiles( ...
    otherSubjects, otherCrunchedDir, otherCohortName, repoRoot, otherTaskType, ...
    signalType, sigSource, S_condition, N_condition, topPct, auxSearchDirs, 'broadband');

fprintf('\nTop-%d%% channels:  local=%d  other=%d\n', ...
    round(topPct * 100), numel(profilesLocal), numel(profilesOther));

rankTable = [summaryLocal; summaryOther];
writetable(rankTable, fullfile(outputDir, 'top_lang_electrode_ranking.csv'));

if isempty(profilesLocal) || isempty(profilesOther)
    warning('Not enough top channels for comparison.');
    save(fullfile(outputDir, 'top_lang_comparison.mat'), ...
        'profilesLocal', 'profilesOther', 'rankLocal', 'rankOther', ...
        'rankTable', 'topPct', 'localTaskType', 'otherTaskType', '-v7.3');
    return;
end

save(fullfile(outputDir, 'top_lang_comparison.mat'), ...
    'profilesLocal', 'profilesOther', 'rankTable', ...
    'rankLocal', 'rankOther', 'topPct', 'localTaskType', 'otherTaskType', ...
    'signalType', '-v7.3');

if doPlots
    plot_cohort_top_averages(profilesLocal, profilesOther, outputDir, signalType, topPct);
    plot_individual_top_traces(profilesLocal, outputDir, localTaskType, signalType, topPct);
    plot_individual_top_traces(profilesOther, outputDir, otherTaskType, signalType, topPct);
end

fprintf('\nSaved results to: %s\n', outputDir);
end


function [profiles, allRankings, summaryRows] = load_top_profiles( ...
    subjects, crunchedDir, cohortName, repoRoot, taskType, ...
    signalType, sigSource, S_condition, N_condition, topPct, auxSearchDirs, cohortKind)

profiles = struct([]);
allRankings = struct([]);
summaryRows = table();

for s = 1:numel(subjects)
    subj = subjects{s};
    try
        [sn_obj, meta] = load_sn_obj_for_comparison(subj, crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, cohortKind);
    catch ME
        warning('Skipping %s: %s', subj, ME.message);
        continue;
    end

    rankArgs = {'signalType', signalType, 'sigSource', sigSource, ...
        'S_condition', S_condition, 'topPct', topPct};
    if ~isempty(N_condition)
        rankArgs = [rankArgs, {'N_condition', N_condition}]; %#ok<AGROW>
    end
    try
        R = compute_lang_elec_ranking(sn_obj, rankArgs{:});
    catch ME
        warning('Ranking failed for %s: %s', subj, ME.message);
        continue;
    end

    if isempty(allRankings)
        allRankings = R;
    else
        allRankings(end + 1) = R; %#ok<AGROW>
    end

    summaryRows = [summaryRows; ranking_to_rows(R, cohortName, taskType)]; %#ok<AGROW>

    if isempty(R.top_chan_idx)
        fprintf('  %s (%s): 0 sig -> 0 top\n', meta.subject, cohortName);
        continue;
    end

    profArgs = {'cohort', cohortName, 'subject', meta.subject, ...
        'signalType', signalType, 'sigSource', sigSource, ...
        'words', R.words, 'S_condition', S_condition, 'sigOnly', true, ...
        'use_odd_for_inference', false, 'chanIdx', R.top_chan_idx, ...
        'requireMni', false};
    if ~isempty(N_condition)
        profArgs = [profArgs, {'N_condition', N_condition}]; %#ok<AGROW>
    end

    try
        P = build_lang_channel_profiles(sn_obj, profArgs{:});
    catch ME
        warning('Profile build failed for %s: %s', subj, ME.message);
        continue;
    end

    for k = 1:numel(P)
        ci = P(k).chanIdx;
        P(k).taskType = taskType;
        P(k).s_vs_n_corr = R.s_vs_n_corr(ci);
        P(k).heldout_corr = R.heldout_corr(ci);
        P(k).split_half_reliab = R.split_half_reliab(ci);
        P(k).composite_score = R.composite_score(ci);
        P(k).rank_among_sig = R.rank_among_sig(ci);
    end

    fprintf('  %s (%s): %d sig -> %d top | %s\n', ...
        meta.subject, cohortName, R.n_sig, numel(P), meta.file);
    profiles = append_profiles(profiles, P);
end
end


function T = ranking_to_rows(R, cohortName, taskType)
sigIdx = find(R.sig_mask);
nSig = numel(sigIdx);

subjCol = strings(nSig, 1);
cohortCol = strings(nSig, 1);
taskCol = strings(nSig, 1);
labelCol = strings(nSig, 1);
chanCol = zeros(nSig, 1);
isTopCol = false(nSig, 1);
corrCol = zeros(nSig, 1);
heldCol = zeros(nSig, 1);
relCol = zeros(nSig, 1);
compCol = zeros(nSig, 1);
rankCol = zeros(nSig, 1);
pCol = zeros(nSig, 1);

for i = 1:nSig
    ci = sigIdx(i);
    subjCol(i) = string(R.subject);
    cohortCol(i) = string(cohortName);
    taskCol(i) = string(taskType);
    labelCol(i) = string(R.labels{ci});
    chanCol(i) = ci;
    isTopCol(i) = R.top_mask(ci);
    corrCol(i) = R.s_vs_n_corr(ci);
    heldCol(i) = R.heldout_corr(ci);
    relCol(i) = R.split_half_reliab(ci);
    compCol(i) = R.composite_score(ci);
    rankCol(i) = R.rank_among_sig(ci);
    pCol(i) = R.p_ratio(ci);
end

T = table(subjCol, cohortCol, taskCol, labelCol, chanCol, isTopCol, ...
    corrCol, heldCol, relCol, compCol, rankCol, pCol, ...
    'VariableNames', {'subject', 'cohort', 'task', 'label', 'chan_idx', 'is_top', ...
    's_vs_n_corr', 'heldout_corr', 'split_half_reliab', 'composite_score', ...
    'rank_among_sig', 'p_ratio'});
end


function plot_cohort_top_averages(profilesLocal, profilesOther, outputDir, signalType, topPct)
[meanS_L, meanN_L, meanDiff_L, semS_L, semN_L, semDiff_L, tL] = cohort_average_traces(profilesLocal);
[meanS_O, meanN_O, meanDiff_O, semS_O, semN_O, semDiff_O, tO] = cohort_average_traces(profilesOther);

fig = figure('Color', 'w', 'Position', [80 80 1300 480], 'Visible', 'off');
tiled = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled);
plot(ax1, tL, meanS_L, 'Color', [0 0.45 0.74], 'LineWidth', 2); hold(ax1, 'on');
plot(ax1, tL, meanN_L, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'LineStyle', '--');
shaded_error(ax1, tL, meanS_L, semS_L, [0 0.45 0.74]);
shaded_error(ax1, tL, meanN_L, semN_L, [0.85 0.33 0.1]);
grid(ax1, 'on');
title(ax1, sprintf('MITSWJNTask top-%d%% (n=%d)', round(topPct*100), numel(profilesLocal)));
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'HG (z)');
legend(ax1, {profilesLocal(1).S_condition, profilesLocal(1).N_condition}, 'Location', 'best');

ax2 = nexttile(tiled);
plot(ax2, tO, meanS_O, 'Color', [0 0.45 0.74], 'LineWidth', 2); hold(ax2, 'on');
plot(ax2, tO, meanN_O, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'LineStyle', '--');
shaded_error(ax2, tO, meanS_O, semS_O, [0 0.45 0.74]);
shaded_error(ax2, tO, meanN_O, semN_O, [0.85 0.33 0.1]);
grid(ax2, 'on');
title(ax2, sprintf('MITLangloc top-%d%% (n=%d)', round(topPct*100), numel(profilesOther)));
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'HG (z)');

ax3 = nexttile(tiled);
plot(ax3, tL, meanDiff_L, 'b', 'LineWidth', 2); hold(ax3, 'on');
plot(ax3, tO, meanDiff_O, 'r', 'LineWidth', 2);
shaded_error(ax3, tL, meanDiff_L, semDiff_L, [0 0.45 0.74]);
shaded_error(ax3, tO, meanDiff_O, semDiff_O, [0.85 0.33 0.1]);
yline(ax3, 0, 'k:');
grid(ax3, 'on');
title(ax3, 'S - N grand average');
xlabel(ax3, 'Time (s)'); ylabel(ax3, 'S - N (z)');
legend(ax3, {'MITSWJNTask', 'MITLangloc'}, 'Location', 'best');

pngName = fullfile(outputDir, sprintf('cohort_top%d_grand_average_%s.png', ...
    round(topPct*100), signalType));
saveas(fig, pngName);
close(fig);
fprintf('Saved cohort grand average: %s\n', pngName);
end


function plot_individual_top_traces(profiles, outputDir, taskType, signalType, topPct)
nPlot = min(numel(profiles), 30);
if nPlot < 1
    return;
end
nCol = 5;
nRow = ceil(nPlot / nCol);

fig = figure('Color', 'w', 'Position', [80 80 1400 900], 'Visible', 'off');
tiled = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:nPlot
    ax = nexttile(tiled);
    P = profiles(k);
    if isempty(P.word_ave_S)
        continue;
    end
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

sgtitle(tiled, sprintf('%s top-%d%% S-N traces (%s, n=%d)', ...
    taskType, round(topPct*100), signalType, numel(profiles)));
pngName = fullfile(outputDir, sprintf('%s_top%d_individual_S_minus_N_%s.png', ...
    taskType, round(topPct*100), signalType));
saveas(fig, pngName);
close(fig);
fprintf('Saved individual traces: %s\n', pngName);
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


function subjects = discover_cohort_subjects(crunchedDir, taskType, cohortKind)
if ~isfolder(crunchedDir)
    subjects = {};
    return;
end
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
    if ~has_sn_obj(fullfile(d(i).folder, d(i).name))
        continue;
    end
    subj = subject_from_filename(d(i).name, taskType, cohortKind);
    if ~isempty(subj)
        subjects{end+1} = subj; %#ok<AGROW>
    end
end
subjects = unique(subjects, 'stable');
end


function ok = has_sn_obj(matFile)
ok = false;
try
    vars = whos('-file', matFile);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function subj = subject_from_filename(fname, taskType, cohortKind)
subj = '';
switch lower(cohortKind)
    case 'brainstorm'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched\.mat$'], 'tokens', 'once');
        if ~isempty(tok), subj = tok{1}; return; end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched\.mat$', 'tokens', 'once');
        if ~isempty(tok), subj = tok{1}; end
    case 'broadband'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched_HG(_ZScore)?\.mat$'], 'tokens', 'once');
        if ~isempty(tok), subj = tok{1}; return; end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched_HG(_ZScore)?\.mat$', 'tokens', 'once');
        if ~isempty(tok), subj = tok{1}; end
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
