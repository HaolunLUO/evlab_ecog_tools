function test_encoding_group_heldout_corr()
% TEST_ENCODING_GROUP_HELDOUT_CORR
%   Test whether encoding_group categories align with language-localizer
%   metrics from MITSWJNTask crunched data:
%     heldout_corr (even-trial S vs N rho)
%     rho_odd (odd-trial S vs N rho)
%     oddeven_trial_reliab (fixed odd/even split-half, S+N trials pooled)
%     oddeven_rho_abs_gap = |rho_odd - rho_even|
%
%   Primary hypothesis (loc_only vs both):
%     loc_only -> lower heldout_corr / lower odd-even reliability
%     both     -> higher heldout_corr / higher odd-even reliability
%
%   Reads:  <workingDir>/group_encoding_results/functional_taxonomy_electrode_table.csv
%   Loads:  <workingDir>/<Subject>_MITSWJNTask_crunched.mat
%   Saves:  <workingDir>/group_encoding_results/encoding_group_heldout_corr_test/

%% USER SETTINGS
repoRoot    = fileparts(fileparts(mfilename('fullpath')));
workingDir  = 'F:\seeg\luohong\analysisEV';
taskType    = 'MITSWJNTask';
signalType  = 'bipolar';
sigSource   = 's_vs_n';
S_condition = 'Sentences';
N_condition = '';
encFeature  = 'gpt2cn_l24';  % must match encoding_group sig column

electrodeCsv = fullfile(workingDir, 'group_encoding_results', ...
    'functional_taxonomy_electrode_table.csv');
outputDir = fullfile(workingDir, 'group_encoding_results', ...
    'encoding_group_heldout_corr_test');
crunchedDir = workingDir;
auxSearchDirs = {workingDir, fullfile(workingDir, 'output', taskType)};

%% SETUP
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

if ~isfile(electrodeCsv)
    error('Electrode table not found:\n  %s', electrodeCsv);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

encTbl = readtable(electrodeCsv, 'TextType', 'string');
requiredCols = ["subject", "channel", "ch_idx", "category", "is_lang"];
missingCols = requiredCols(~ismember(requiredCols, string(encTbl.Properties.VariableNames)));
if ~isempty(missingCols)
    error('Missing columns in electrode table: %s', strjoin(missingCols, ', '));
end

sigCol = "sig_" + string(encFeature);
if ~ismember(sigCol, string(encTbl.Properties.VariableNames))
    error('Encoding significance column not found: %s', sigCol);
end

subjects = unique(encTbl.subject);
fprintf('\n=== ENCODING GROUP vs LANGLOC METRICS ===\n');
fprintf('Electrode table: %s\n', electrodeCsv);
fprintf('Subjects: %s\n', strjoin(subjects, ', '));

%% Attach langloc metrics per channel (heldout + odd/even reliability)
nRows = height(encTbl);
heldoutCol = nan(nRows, 1);
sVsNCol = nan(nRows, 1);
rhoOddCol = nan(nRows, 1);
rhoEvenCol = nan(nRows, 1);
oddevenRelCol = nan(nRows, 1);
oddevenGapCol = nan(nRows, 1);
isLangSigCol = false(nRows, 1);

for s = 1:numel(subjects)
    subj = subjects(s);
    try
        [sn_obj, meta] = load_sn_obj_for_comparison(char(subj), crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, 'brainstorm');
    catch ME
        warning('Skipping %s: %s', subj, ME.message);
        continue;
    end

    held = extract_heldout_corr(sn_obj, signalType);
    sigMask = get_lang_sig_mask(sn_obj, signalType, sigSource);
    corrVec = extract_s_vs_n_corr(sn_obj, signalType);
    relVec = nan(channel_count(sn_obj, signalType), 1);
    rhoOddVec = nan(channel_count(sn_obj, signalType), 1);
    try
        [relVec, rhoOddVec] = extract_oddeven_lang_metrics(sn_obj, signalType, S_condition, N_condition);
    catch ME
        warning('oddeven metrics skipped for %s: %s', subj, ME.message);
    end

    subMask = encTbl.subject == subj;
    subRows = find(subMask);
    for i = 1:numel(subRows)
        ri = subRows(i);
        ci = encTbl.ch_idx(ri);
        if isnan(ci) || ci < 1
            continue;
        end
        if ci <= numel(held)
            heldoutCol(ri) = pick_scalar(held, ci);
            rhoEvenCol(ri) = pick_scalar(held, ci);
        end
        if ci <= numel(corrVec) && isfinite(pick_scalar(corrVec, ci))
            rhoOddCol(ri) = pick_scalar(corrVec, ci);
            sVsNCol(ri) = pick_scalar(corrVec, ci);
        elseif ci <= numel(rhoOddVec)
            rhoOddCol(ri) = pick_scalar(rhoOddVec, ci);
            sVsNCol(ri) = pick_scalar(rhoOddVec, ci);
        end
        if ci <= numel(relVec)
            oddevenRelCol(ri) = pick_scalar(relVec, ci);
        end
        if isfinite(rhoOddCol(ri)) && isfinite(rhoEvenCol(ri))
            oddevenGapCol(ri) = abs(rhoOddCol(ri) - rhoEvenCol(ri));
        end
        if ci <= numel(sigMask)
            isLangSigCol(ri) = logical(pick_scalar(sigMask, ci));
        end
    end
    fprintf('  %s: heldout=%d  oddeven_reliab=%d / %d | %s\n', ...
        subj, sum(isfinite(heldoutCol(subMask))), ...
        sum(isfinite(oddevenRelCol(subMask))), sum(subMask), meta.file);
end

encTbl.heldout_corr = heldoutCol;
encTbl.rho_odd = rhoOddCol;
encTbl.rho_even = rhoEvenCol;
encTbl.oddeven_trial_reliab = oddevenRelCol;
encTbl.oddeven_rho_abs_gap = oddevenGapCol;
encTbl.s_vs_n_corr = sVsNCol;
encTbl.is_lang_sig_ranking = isLangSigCol;

writetable(encTbl, fullfile(outputDir, 'electrodes_with_heldout_corr.csv'));

%% Focus on loc_only vs both (user hypothesis)
testMask = ismember(encTbl.category, ["loc_only", "both"]);
testTbl = encTbl(testMask, :);
testTbl = testTbl(isfinite(testTbl.heldout_corr), :);

fprintf('\nChannels with finite metrics:\n');
for cat = ["loc_only", "both", "gpt2cn_only", "neither"]
    nHeld = sum(encTbl.category == cat & isfinite(encTbl.heldout_corr));
    nRel = sum(encTbl.category == cat & isfinite(encTbl.oddeven_trial_reliab));
    fprintf('  %-12s  heldout=%d  oddeven_reliab=%d\n', cat, nHeld, nRel);
end

statsOut = struct();
metricDefs = {
    'heldout_corr',         'rho_even (held-out S vs N)',           'left',  'both higher';
    'rho_odd',              'rho_odd (inference S vs N)',           'left',  'both higher';
    'oddeven_trial_reliab', 'fixed odd/even trial reliability (S+N pooled)', 'left', 'both higher';
    'oddeven_rho_abs_gap',  '|rho_odd - rho_even| (lower=stable)',  'right', 'loc_only higher (less stable)'
    };

for m = 1:size(metricDefs, 1)
    fieldName = metricDefs{m, 1};
    label = metricDefs{m, 2};
    tail = metricDefs{m, 3};
    subTbl = encTbl(testMask, :);
    subTbl = subTbl(isfinite(subTbl.(fieldName)), :);
    locVals = subTbl.(fieldName)(subTbl.category == "loc_only");
    bothVals = subTbl.(fieldName)(subTbl.category == "both");
    res = compare_loc_vs_both(locVals, bothVals, label, tail);
    statsOut.(fieldName) = res;
end

% Kruskal-Wallis across all four categories (finite heldout only)
allMask = isfinite(encTbl.heldout_corr);
allCats = encTbl.category(allMask);
allHeld = encTbl.heldout_corr(allMask);
catList = ["both", "gpt2cn_only", "loc_only", "neither"];
groups = cell(1, numel(catList));
for k = 1:numel(catList)
    groups{k} = allHeld(allCats == catList(k));
end
validGroups = groups(~cellfun(@(g) numel(g) < 2, groups));
if numel(validGroups) >= 2
  groupVec = [];
  labelVec = [];
  for g = 1:numel(validGroups)
      groupVec = [groupVec; validGroups{g}(:)]; %#ok<AGROW>
      labelVec = [labelVec; g * ones(numel(validGroups{g}), 1)]; %#ok<AGROW>
  end
  [pKw, ~, kwStats] = kruskalwallis(groupVec, labelVec, 'off');
  statsOut.kruskal_p = pKw;
  if isfield(kwStats, 'chi2stat')
      statsOut.kruskal_h = kwStats.chi2stat;
  elseif isfield(kwStats, 'chistat')
      statsOut.kruskal_h = kwStats.chistat;
  else
      statsOut.kruskal_h = NaN;
  end
  fprintf('\nKruskal-Wallis (all categories): p=%.4g  chi2=%.3f\n', ...
      pKw, statsOut.kruskal_h);
end

save(fullfile(outputDir, 'heldout_corr_test_stats.mat'), 'statsOut', 'testTbl', 'encTbl');

%% Plots
plot_category_violin(testTbl, outputDir, 'heldout_corr', ...
    'heldout\_corr (even-trial S vs N \rho)', 'heldout_corr_loc_only_vs_both.png');
plot_category_violin(encTbl(encTbl.category == "loc_only" | encTbl.category == "both", :), ...
    outputDir, 'oddeven_trial_reliab', ...
    'Odd/even trial reliability', 'oddeven_reliab_loc_only_vs_both.png');
plot_category_violin(encTbl(encTbl.category == "loc_only" | encTbl.category == "both", :), ...
    outputDir, 'oddeven_rho_abs_gap', ...
    '|rho_{odd} - rho_{even}|', 'oddeven_rho_gap_loc_only_vs_both.png');
plot_oddeven_scatter(encTbl, outputDir);
plot_oddeven_trial_scatter_category_average(encTbl, subjects, crunchedDir, repoRoot, ...
    taskType, auxSearchDirs, signalType, S_condition, N_condition, outputDir);
plot_oddeven_trial_scatter_typical_channel(encTbl, subjects, crunchedDir, repoRoot, ...
    taskType, auxSearchDirs, signalType, S_condition, N_condition, outputDir);
plot_all_categories_violin(encTbl, outputDir);

fprintf('\nSaved results to: %s\n', outputDir);
end


function [reliabVec, rhoOddVec] = extract_oddeven_lang_metrics(sn_obj, signalType, S_condition, N_condition)
% Fixed odd/even metrics from per-trial HG:
%   oddeven_trial_reliab - corr(odd, even) on all trials pooled (S + N)
%   rho_odd             - S vs N Spearman rho on odd trials only
nChan = channel_count(sn_obj, signalType);
reliabVec = nan(nChan, 1);
rhoOddVec = nan(nChan, 1);

[S_data, N_data, ~, ~] = load_channel_trial_tensors( ...
    sn_obj, signalType, S_condition, N_condition, []);
if isempty(S_data)
    return;
end
nElec = size(S_data, 1);

for e = 1:min(nElec, nChan)
    reliabVec(e) = fixed_oddeven_trial_reliab_pooled(S_data(e, :, :), N_data(e, :, :));
    rhoOddVec(e) = sn_rho_fixed_half_3d(S_data(e, :, :), N_data(e, :, :), 'odd');
end
end


function [S_data, N_data, S_cond, N_cond] = load_channel_trial_tensors( ...
    sn_obj, signalType, S_condition, N_condition, words)
avail = unique(sn_obj.condition, 'stable');
S_cond = resolve_condition_name(avail, ...
    {S_condition, 'Sentences', 'SENTENCES', 'sentence'}, 'S');
if isempty(N_condition)
    try
        N_cond = resolve_condition_name(avail, ...
            {'Jabberwocky', 'JABBERWOCKY', 'jabberwocky'}, 'N');
    catch
        N_cond = resolve_condition_name(avail, ...
            {'Nonword-lists', 'Nonwords', 'NONWORDS', 'nonwords'}, 'N');
    end
else
    N_cond = resolve_condition_name(avail, ...
        {N_condition, 'Jabberwocky', 'Nonword-lists', 'Nonwords', 'NONWORDS'}, 'N');
end
if nargin < 5 || isempty(words)
    words = default_words_for_sn(sn_obj, S_cond);
end

[S_tbl, ~] = sn_obj.get_ave_cond_trial('words', words, 'condition', S_cond);
[N_tbl, ~] = sn_obj.get_ave_cond_trial('words', words, 'condition', N_cond);

if strcmpi(signalType, 'bipolar')
    dataField = 'bip_elec_data';
else
    dataField = 'elec_data';
end
if ~ismember(dataField, S_tbl.Properties.VariableNames) ...
        || isempty(S_tbl.(dataField){1})
    S_data = [];
    N_data = [];
    return;
end
S_data = S_tbl.(dataField){1};
N_data = N_tbl.(dataField){1};
end


function r = fixed_oddeven_trial_reliab_pooled(S3d, N3d)
% Pool S and N per-trial means, then fixed odd/even split-half reliability.
S_means = squeeze(mean(S3d, 3));
N_means = squeeze(mean(N3d, 3));
if isrow(S_means); S_means = S_means(:); end
if isrow(N_means); N_means = N_means(:); end
trialMeans = [S_means; N_means];
r = fixed_oddeven_trial_reliab_vector(trialMeans);
end


function r = fixed_oddeven_trial_reliab_vector(trialMeans)
trialMeans = trialMeans(:);
n = numel(trialMeans);
if n < 4
    r = NaN;
    return;
end
odd = trialMeans(1:2:end);
even = trialMeans(2:2:end);
nHalf = min(numel(odd), numel(even));
if nHalf < 2
    r = NaN;
    return;
end
r = corr(odd(1:nHalf), even(1:nHalf), 'Type', 'Spearman', 'Rows', 'pairwise');
if ~isscalar(r)
    r = r(1);
end
end


function rho = sn_rho_fixed_half_3d(S3d, N3d, whichHalf)
S_means = squeeze(mean(S3d, 3));
N_means = squeeze(mean(N3d, 3));
if isrow(S_means); S_means = S_means(:); end
if isrow(N_means); N_means = N_means(:); end
if strcmpi(whichHalf, 'odd')
    S_h = S_means(1:2:end);
    N_h = N_means(1:2:end);
else
    S_h = S_means(2:2:end);
    N_h = N_means(2:2:end);
end
if numel(S_h) < 2 || numel(N_h) < 2
    rho = NaN;
    return;
end
y = [S_h; N_h];
x = [ones(numel(S_h), 1); -ones(numel(N_h), 1)];
rho = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
if ~isscalar(rho); rho = rho(1); end
end


function r = fixed_oddeven_trial_reliab_3d(data3d)
% Single-condition fixed odd/even reliability (kept for reference).
trialMeans = squeeze(mean(data3d, 3));
if isrow(trialMeans)
    trialMeans = trialMeans(:);
end
r = fixed_oddeven_trial_reliab_vector(trialMeans);
end


function words = default_words_for_sn(sn_obj, S_condition)
words = 1:12;
if isprop(sn_obj, 's_vs_n_ops') && ~isempty(sn_obj.s_vs_n_ops) ...
        && isstruct(sn_obj.s_vs_n_ops) && isfield(sn_obj.s_vs_n_ops, 'words')
    words = sn_obj.s_vs_n_ops.words;
end
maxW = 0;
for w = 1:24
    try
        sn_obj.get_ave_cond_trial('words', w, 'condition', S_condition);
        maxW = w;
    catch
        break;
    end
end
if maxW >= 1
    words = words(words >= 1 & words <= maxW);
    if isempty(words)
        words = 1:maxW;
    end
end
end


function n = channel_count(sn_obj, signalType)
if strcmpi(signalType, 'bipolar')
    n = numel(sn_obj.bip_ch_label);
else
    n = numel(sn_obj.elec_ch_label);
end
end


function res = compare_loc_vs_both(locVals, bothVals, label, tail)
fprintf('\n--- loc_only vs both: %s ---\n', label);
fprintf('  loc_only: n=%d  median=%.4f  mean=%.4f\n', ...
    numel(locVals), median(locVals, 'omitnan'), mean(locVals, 'omitnan'));
fprintf('  both:     n=%d  median=%.4f  mean=%.4f\n', ...
    numel(bothVals), median(bothVals, 'omitnan'), mean(bothVals, 'omitnan'));

res = struct();
res.label = label;
res.n_loc_only = numel(locVals);
res.n_both = numel(bothVals);
res.median_loc_only = median(locVals, 'omitnan');
res.median_both = median(bothVals, 'omitnan');
res.mean_loc_only = mean(locVals, 'omitnan');
res.mean_both = mean(bothVals, 'omitnan');

if numel(locVals) >= 2 && numel(bothVals) >= 2
    [pTwo, ~, statsTwo] = ranksum(locVals, bothVals);
    res.ranksum_p_two_sided = pTwo;
    res.ranksum_z = statsTwo.zval;
    res.cliffs_delta = cliffs_delta(locVals, bothVals);
    if strcmpi(tail, 'left')
        [pOne, ~, ~] = ranksum(locVals, bothVals, 'Tail', 'left');
        res.ranksum_p_directional = pOne;
        dirNote = 'loc < both';
    elseif strcmpi(tail, 'right')
        [pOne, ~, ~] = ranksum(locVals, bothVals, 'Tail', 'right');
        res.ranksum_p_directional = pOne;
        dirNote = 'loc > both';
    else
        pOne = NaN;
        res.ranksum_p_directional = NaN;
        dirNote = '';
    end
    fprintf('  Mann-Whitney U (two-sided): p=%.4g  z=%.3f\n', pTwo, statsTwo.zval);
    if ~isnan(pOne)
        fprintf('  Mann-Whitney U (%s):      p=%.4g\n', dirNote, pOne);
    end
    fprintf('  Cliff''s delta:             %.3f\n', res.cliffs_delta);
else
    warning('Not enough channels for rank-sum test: %s', label);
end
end


function heldoutVec = extract_heldout_corr(sn_obj, signalType)
if strcmpi(signalType, 'bipolar')
    nChan = numel(sn_obj.bip_ch_label);
else
    nChan = numel(sn_obj.elec_ch_label);
end
heldoutVec = nan(nChan, 1);
if ~isprop(sn_obj, 'hg_sn_corr') || isempty(sn_obj.hg_sn_corr) ...
        || ~isfield(sn_obj.hg_sn_corr, 'results')
    return;
end
if strcmpi(signalType, 'bipolar') ...
        && isfield(sn_obj.hg_sn_corr.results, 'bipolar')
    raw = sn_obj.hg_sn_corr.results.bipolar.corr;
elseif isfield(sn_obj.hg_sn_corr.results, 'unipolar')
    raw = sn_obj.hg_sn_corr.results.unipolar.corr;
else
    return;
end
raw = raw(:);
heldoutVec(1:min(nChan, numel(raw))) = raw(1:min(nChan, numel(raw)));
end


function corrVec = extract_s_vs_n_corr(sn_obj, signalType)
if strcmpi(signalType, 'bipolar')
    nChan = numel(sn_obj.bip_ch_label);
    dataField = 'bip_elec_data';
else
    nChan = numel(sn_obj.elec_ch_label);
    dataField = 'elec_data';
end
corrVec = nan(nChan, 1);
if istable(sn_obj.s_vs_n_corr) && ismember(dataField, sn_obj.s_vs_n_corr.Properties.VariableNames)
    raw = sn_obj.s_vs_n_corr.(dataField){1, 1};
    raw = raw(:);
    corrVec(1:min(nChan, numel(raw))) = raw(1:min(nChan, numel(raw)));
end
end


function v = pick_scalar(vec, idx)
v = NaN;
if idx < 1 || idx > numel(vec)
    return;
end
v = vec(idx);
if ~isscalar(v)
    v = v(1);
end
end


function d = cliffs_delta(x, y)
% Cliff's delta: P(x>y) - P(x<y); positive => y tends larger than x
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


function plot_category_violin(testTbl, outputDir, fieldName, yLabel, fileName)
cats = ["loc_only", "both"];
displayNames = {'Lan-Loc only', 'Both'};
colors = [31 119 180; 214 39 40] / 255;
if ~ismember(fieldName, testTbl.Properties.VariableNames)
    return;
end
testTbl = testTbl(isfinite(testTbl.(fieldName)), :);

fig = figure('Color', 'w', 'Position', [100 100 520 420], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');
positions = 1:2;
for k = 1:2
    vals = testTbl.(fieldName)(testTbl.category == cats(k));
    if numel(vals) < 1
        continue;
    end
    if numel(vals) >= 2
        [f, xi] = ksdensity(vals);
        f = f / max(f) * 0.35;
        fill(ax, [positions(k) - f, positions(k) + fliplr(f)], ...
            [xi, fliplr(xi)], colors(k, :), 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    end
    boxchart(ax, positions(k) * ones(size(vals)), vals, ...
        'BoxFaceColor', colors(k, :), 'BoxFaceAlpha', 0.2, ...
        'LineWidth', 1.2, 'MarkerStyle', '.', 'JitterOutliers', 'on', ...
        'BoxWidth', 0.18);
end
set(ax, 'XTick', positions, 'XTickLabel', displayNames);
ylabel(ax, yLabel);
title(ax, strrep(yLabel, '\', ''));
grid(ax, 'on');
saveas(fig, fullfile(outputDir, fileName));
close(fig);
end


function plot_oddeven_scatter(encTbl, outputDir)
mask = isfinite(encTbl.rho_odd) & isfinite(encTbl.rho_even);
if ~any(mask)
    return;
end
cats = ["loc_only", "both", "gpt2cn_only", "neither"];
colors = [31 119 180; 214 39 40; 255 127 14; 187 187 187] / 255;

fig = figure('Color', 'w', 'Position', [100 100 560 480], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');
for k = 1:numel(cats)
    cm = mask & encTbl.category == cats(k);
    if ~any(cm)
        continue;
    end
    scatter(ax, encTbl.rho_odd(cm), encTbl.rho_even(cm), 36, ...
        colors(k, :), 'filled', 'DisplayName', char(cats(k)));
end
lims = [min([ax.XLim ax.YLim]) max([ax.XLim ax.YLim])];
plot(ax, lims, lims, 'k--');
axis(ax, 'square');
xlabel(ax, '\rho_{odd} (inference trials)');
ylabel(ax, '\rho_{even} (held-out trials)');
title(ax, 'Odd vs even S-vs-N correlation by encoding group');
legend(ax, 'Location', 'best');
grid(ax, 'on');
saveas(fig, fullfile(outputDir, 'rho_odd_vs_rho_even_scatter.png'));
close(fig);
end


function plot_all_categories_violin(encTbl, outputDir)
cats = ["both", "gpt2cn_only", "loc_only", "neither"];
displayNames = {'Both', 'GPT2-CN only', 'Lan-Loc only', 'Neither'};
colors = [214 39 40; 255 127 14; 31 119 180; 187 187 187] / 255;
mask = isfinite(encTbl.heldout_corr);

fig = figure('Color', 'w', 'Position', [100 100 720 420], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');
positions = 1:numel(cats);
for k = 1:numel(cats)
    vals = encTbl.heldout_corr(mask & encTbl.category == cats(k));
    if isempty(vals)
        continue;
    end
    if numel(vals) >= 2
        [f, xi] = ksdensity(vals);
        f = f / max(f) * 0.3;
        fill(ax, [positions(k) - f, positions(k) + fliplr(f)], ...
            [xi, fliplr(xi)], colors(k, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    boxchart(ax, positions(k) * ones(size(vals)), vals, ...
        'BoxFaceColor', colors(k, :), 'BoxFaceAlpha', 0.15, ...
        'LineWidth', 1.1, 'MarkerStyle', '.', 'JitterOutliers', 'on', ...
        'BoxWidth', 0.16);
end
set(ax, 'XTick', positions, 'XTickLabel', displayNames);
ylabel(ax, 'heldout\_corr');
title(ax, 'heldout\_corr by encoding\_group category');
grid(ax, 'on');
saveas(fig, fullfile(outputDir, 'heldout_corr_all_categories.png'));
close(fig);
end


function plot_oddeven_trial_scatter_category_average(encTbl, subjects, crunchedDir, ...
    repoRoot, taskType, auxSearchDirs, signalType, S_condition, N_condition, outputDir)
% Grand-average odd/even trial pattern per category: z-score each channel's
% trial profile, average across channels, then plot S vs N colored points.

catOrder = ["both", "gpt2cn_only", "loc_only", "neither"];
catTitles = {'Both', 'GPT2-CN sig only', 'Lan-Loc only', 'Neither'};
condColors = struct('S', [0 0.45 0.74], 'N', [0.85 0.33 0.1]);
markerSize = 52;
alphaPt = 0.9;

subjCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
summaryRows = table();
fig = figure('Color', 'w', 'Position', [60 60 1100 900], 'Visible', 'off');
tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

S_label = 'Sentences (S)';
N_label = 'Nonwords (N)';

for k = 1:numel(catOrder)
    ax = nexttile(tiled);
    hold(ax, 'on');
    cat = catOrder(k);

    catMask = encTbl.category == cat & isfinite(encTbl.oddeven_trial_reliab);
    if ~any(catMask)
        title(ax, sprintf('%s\n(no channels)', catTitles{k}), 'FontSize', 10);
        grid(ax, 'on');
        continue;
    end

    catTbl = encTbl(catMask, :);
    medR = median(catTbl.oddeven_trial_reliab);
    S_vecs = {};
    N_vecs = {};

    for i = 1:height(catTbl)
        subj = char(catTbl.subject(i));
        ci = catTbl.ch_idx(i);
        try
            trialData = get_cached_trial_tensors(subj, crunchedDir, repoRoot, ...
                taskType, auxSearchDirs, signalType, S_condition, N_condition, subjCache);
        catch
            continue;
        end
        if i == 1
            S_label = short_cond_label(trialData.S_label);
            N_label = short_cond_label(trialData.N_label);
        end
        if isnan(ci) || ci < 1 || ci > size(trialData.S_data, 1)
            continue;
        end
        [sVec, nVec] = channel_trial_mean_vectors(trialData.S_data, trialData.N_data, ci);
        S_vecs{end + 1} = zscore_trial_vector(sVec); %#ok<AGROW>
        N_vecs{end + 1} = zscore_trial_vector(nVec); %#ok<AGROW>
    end

    nCh = numel(S_vecs);
    if nCh < 1
        title(ax, sprintf('%s\n(no trial data)', catTitles{k}), 'FontSize', 10);
        grid(ax, 'on');
        continue;
    end

    nTrialsS = min(cellfun(@numel, S_vecs));
    nTrialsN = min(cellfun(@numel, N_vecs));
    S_stack = zeros(nCh, nTrialsS);
    N_stack = zeros(nCh, nTrialsN);
    for j = 1:nCh
        S_stack(j, :) = S_vecs{j}(1:nTrialsS);
        N_stack(j, :) = N_vecs{j}(1:nTrialsN);
    end

    S_avg = mean(S_stack, 1, 'omitnan');
    N_avg = mean(N_stack, 1, 'omitnan');
    S3d = reshape(S_avg, 1, nTrialsS, 1);
    N3d = reshape(N_avg, 1, nTrialsN, 1);

    [oS, eS, oN, eN] = oddeven_trial_pairs(S3d, N3d);
    rAvg = fixed_oddeven_trial_reliab_pooled(S3d, N3d);
    rS = corr(oS, eS, 'Type', 'Spearman', 'Rows', 'pairwise');
    rN = corr(oN, eN, 'Type', 'Spearman', 'Rows', 'pairwise');
    if ~isscalar(rS); rS = rS(1); end
    if ~isscalar(rN); rN = rN(1); end

    if ~isempty(oS)
        scatter(ax, oS, eS, markerSize, condColors.S, 'filled', ...
            'MarkerFaceAlpha', alphaPt, 'MarkerEdgeColor', condColors.S * 0.7, ...
            'DisplayName', S_label);
    end
    if ~isempty(oN)
        scatter(ax, oN, eN, markerSize, condColors.N, 'filled', ...
            'MarkerFaceAlpha', alphaPt, 'MarkerEdgeColor', condColors.N * 0.7, ...
            'DisplayName', N_label);
    end

    allX = [oS(:); oN(:)];
    allY = [eS(:); eN(:)];
    if ~isempty(allX)
        lims = [min([allX; allY]), max([allX; allY])];
        pad = 0.12 * max(diff(lims), eps);
        lims = [lims(1) - pad, lims(2) + pad];
        plot(ax, lims, lims, 'k--', 'HandleVisibility', 'off');
        xlim(ax, lims);
        ylim(ax, lims);
        axis(ax, 'square');
    end

    title(ax, sprintf(['%s\nmean pattern (r=%.2f, median r=%.2f, n=%d ch)\n' ...
        'r_S=%.2f  r_N=%.2f'], ...
        catTitles{k}, rAvg, medR, nCh, rS, rN), 'FontSize', 9);
    xlabel(ax, 'Odd-trial HG (z-scored, channel-mean)');
    ylabel(ax, 'Even-trial HG (z-scored, channel-mean)');
    grid(ax, 'on');
    if k == 1
        legend(ax, 'Location', 'best');
    end

    summaryRows = [summaryRows; table( ...
        cat, nCh, medR, rAvg, rS, rN, nTrialsS, nTrialsN, ...
        'VariableNames', {'category', 'n_channels', 'median_r_channels', ...
        'r_mean_pattern', 'r_S_mean_pattern', 'r_N_mean_pattern', ...
        'n_S_trials', 'n_N_trials'})]; %#ok<AGROW>
end

sgtitle(fig, ['Category-mean odd/even trial pattern (z-scored per channel, then averaged)' newline ...
    'Color = langloc condition; r = reliability of the mean profile'], ...
    'FontWeight', 'bold');

pngPath = fullfile(outputDir, 'oddeven_trial_scatter_category_average.png');
pdfPath = fullfile(outputDir, 'oddeven_trial_scatter_category_average.pdf');
exportgraphics(fig, pngPath, 'Resolution', 200);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
close(fig);

if ~isempty(summaryRows)
    writetable(summaryRows, fullfile(outputDir, 'oddeven_category_average_summary.csv'));
end
fprintf('Saved category-average oddeven scatter: %s\n', pngPath);
end


function trialData = get_cached_trial_tensors(subj, crunchedDir, repoRoot, ...
    taskType, auxSearchDirs, signalType, S_condition, N_condition, subjCache)
if isKey(subjCache, subj)
    trialData = subjCache(subj);
    return;
end
[sn_obj, ~] = load_sn_obj_for_comparison(subj, crunchedDir, repoRoot, ...
    taskType, auxSearchDirs, 'brainstorm');
[S_data, N_data, S_label, N_label] = load_channel_trial_tensors( ...
    sn_obj, signalType, S_condition, N_condition, []);
if isempty(S_data)
    error('No trial tensors for %s.', subj);
end
trialData = struct('S_data', S_data, 'N_data', N_data, ...
    'S_label', S_label, 'N_label', N_label);
subjCache(subj) = trialData;
end


function [sVec, nVec] = channel_trial_mean_vectors(S_data, N_data, ci)
sVec = squeeze(mean(S_data(ci, :, :), 3));
nVec = squeeze(mean(N_data(ci, :, :), 3));
sVec = sVec(:)';
nVec = nVec(:)';
end


function v = zscore_trial_vector(v)
v = v(:)';
mu = mean(v, 'omitnan');
sd = std(v, 0, 'omitnan');
if sd < eps || ~isfinite(sd)
    v = zeros(size(v));
else
    v = (v - mu) / sd;
end
end


function plot_oddeven_trial_scatter_typical_channel(encTbl, subjects, crunchedDir, ...
    repoRoot, taskType, auxSearchDirs, signalType, S_condition, N_condition, outputDir)
% One exemplar channel per encoding category (closest to category-median
% oddeven_trial_reliab). Shows odd vs even trial HG; color = S vs N condition.

catOrder = ["both", "gpt2cn_only", "loc_only", "neither"];
catTitles = {'Both', 'GPT2-CN sig only', 'Lan-Loc only', 'Neither'};
condColors = struct('S', [0 0.45 0.74], 'N', [0.85 0.33 0.1]);
markerSize = 48;
alphaPt = 0.85;

exemplarRows = table();
fig = figure('Color', 'w', 'Position', [60 60 1100 900], 'Visible', 'off');
tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

S_label = 'Sentences (S)';
N_label = 'Nonwords (N)';

for k = 1:numel(catOrder)
    ax = nexttile(tiled);
    hold(ax, 'on');
    cat = catOrder(k);

    catMask = encTbl.category == cat & isfinite(encTbl.oddeven_trial_reliab);
    if ~any(catMask)
        title(ax, sprintf('%s\n(no channels with reliability)', catTitles{k}), 'FontSize', 10);
        xlabel(ax, 'Odd-trial mean HG');
        ylabel(ax, 'Even-trial mean HG');
        grid(ax, 'on');
        continue;
    end

    catTbl = encTbl(catMask, :);
    medR = median(catTbl.oddeven_trial_reliab);
    [~, pickLocal] = min(abs(catTbl.oddeven_trial_reliab - medR));
    ex = catTbl(pickLocal, :);

    try
        [sn_obj, ~] = load_sn_obj_for_comparison(char(ex.subject), crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, 'brainstorm');
        [S_data, N_data, S_res, N_res] = load_channel_trial_tensors( ...
            sn_obj, signalType, S_condition, N_condition, []);
    catch ME
        title(ax, sprintf('%s\n(load failed: %s)', catTitles{k}, ME.message), 'FontSize', 9);
        continue;
    end

    S_label = short_cond_label(S_res);
    N_label = short_cond_label(N_res);

    ci = ex.ch_idx;
    if isnan(ci) || ci < 1 || ci > size(S_data, 1)
        title(ax, sprintf('%s\n(invalid ch_idx)', catTitles{k}), 'FontSize', 10);
        continue;
    end

    [oS, eS, oN, eN] = oddeven_trial_pairs(S_data(ci, :, :), N_data(ci, :, :));
    rCh = fixed_oddeven_trial_reliab_pooled(S_data(ci, :, :), N_data(ci, :, :));
    rS = corr(oS, eS, 'Type', 'Spearman', 'Rows', 'pairwise');
    rN = corr(oN, eN, 'Type', 'Spearman', 'Rows', 'pairwise');
    if ~isscalar(rS); rS = rS(1); end
    if ~isscalar(rN); rN = rN(1); end

    if ~isempty(oS)
        scatter(ax, oS, eS, markerSize, condColors.S, 'filled', ...
            'MarkerFaceAlpha', alphaPt, 'MarkerEdgeColor', condColors.S * 0.7, ...
            'DisplayName', S_label);
    end
    if ~isempty(oN)
        scatter(ax, oN, eN, markerSize, condColors.N, 'filled', ...
            'MarkerFaceAlpha', alphaPt, 'MarkerEdgeColor', condColors.N * 0.7, ...
            'DisplayName', N_label);
    end

    allX = [oS(:); oN(:)];
    allY = [eS(:); eN(:)];
    if ~isempty(allX)
        lims = [min([allX; allY]), max([allX; allY])];
        pad = 0.08 * max(diff(lims), eps);
        lims = [lims(1) - pad, lims(2) + pad];
        plot(ax, lims, lims, 'k--', 'HandleVisibility', 'off');
        xlim(ax, lims);
        ylim(ax, lims);
        axis(ax, 'square');
    end

    title(ax, sprintf(['%s  |  %s %s\n' ...
        'typical ch (r=%.2f, cat median=%.2f, n=%d ch)\n' ...
        'r_S=%.2f  r_N=%.2f'], ...
        catTitles{k}, ex.subject, ex.channel, rCh, medR, height(catTbl), rS, rN), ...
        'FontSize', 9, 'Interpreter', 'none');
    xlabel(ax, 'Odd-trial mean HG');
    ylabel(ax, 'Even-trial mean HG');
    grid(ax, 'on');
    if k == 1
        legend(ax, 'Location', 'best');
    end

    exemplarRows = [exemplarRows; table( ...
        cat, ex.subject, ex.channel, ex.ch_idx, ex.oddeven_trial_reliab, ...
        medR, rCh, rS, rN, ...
        'VariableNames', {'category', 'subject', 'channel', 'ch_idx', ...
        'oddeven_trial_reliab', 'category_median_r', 'exemplar_r_pooled', ...
        'exemplar_r_S', 'exemplar_r_N'})]; %#ok<AGROW>
end

sgtitle(fig, ['Typical channel per category: odd vs even trial HG' newline ...
    '(exemplar = closest to category-median oddeven\_trial\_reliab; color = condition)'], ...
    'FontWeight', 'bold');

pngPath = fullfile(outputDir, 'oddeven_trial_scatter_typical_channel.png');
pdfPath = fullfile(outputDir, 'oddeven_trial_scatter_typical_channel.pdf');
exportgraphics(fig, pngPath, 'Resolution', 200);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');
close(fig);

if ~isempty(exemplarRows)
    writetable(exemplarRows, fullfile(outputDir, 'oddeven_typical_channel_exemplars.csv'));
end
fprintf('Saved typical-channel oddeven scatter: %s\n', pngPath);
end


function [oddS, evenS, oddN, evenN] = oddeven_trial_pairs(S3d, N3d)
S_means = squeeze(mean(S3d, 3));
N_means = squeeze(mean(N3d, 3));
if isrow(S_means); S_means = S_means(:); end
if isrow(N_means); N_means = N_means(:); end
[oddS, evenS] = split_oddeven_pairs(S_means);
[oddN, evenN] = split_oddeven_pairs(N_means);
end


function [oddV, evenV] = split_oddeven_pairs(trialMeans)
trialMeans = trialMeans(:)';
oddV = trialMeans(1:2:end);
evenV = trialMeans(2:2:end);
n = min(numel(oddV), numel(evenV));
oddV = oddV(1:n);
evenV = evenV(1:n);
end


function lab = short_cond_label(condName)
lab = char(condName);
if contains(lower(lab), 'sent')
    lab = 'Sentences (S)';
elseif contains(lower(lab), 'jabber')
    lab = 'Jabberwocky (N)';
elseif contains(lower(lab), 'nonword')
    lab = 'Nonwords (N)';
end
end
