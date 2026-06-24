function ranking = compute_lang_elec_ranking(sn_obj, varargin)
% COMPUTE_LANG_ELEC_RANKING  Rank language electrodes by correlation + reliability.
%
%   ranking = compute_lang_elec_ranking(sn_obj)
%   ranking = compute_lang_elec_ranking(sn_obj, Name, Value, ...)
%
%   Among language-responsive channels (s_vs_n_sig), computes:
%     - s_vs_n_corr: inference Spearman rho (odd trials from test_s_vs_n)
%     - heldout_corr: held-out S-vs-N rho on even trials (hg_sn_corr if present)
%     - split_half_reliab: split-half reliability of per-trial mean HG
%
%   Composite score = mean of within-subject z-scores of heldout_corr and
%   split_half_reliab (higher = stronger held-out effect + stable response).
%   Falls back to s_vs_n_corr when heldout_corr is unavailable.
%
%   Name-Value:
%     signalType   - 'bipolar' (default) | 'unipolar'
%     sigSource    - 's_vs_n' (default) | 'wordwise'
%     words        - word indices (default [] -> auto from trial timing)
%     S_condition  - default 'Sentences'
%     N_condition  - default '' (auto)
%     topPct       - fraction to mark as top (default 0.10)
%     splitSeed    - RNG seed for split-half (default 42)

p = inputParser();
addParameter(p, 'signalType', 'bipolar');
addParameter(p, 'sigSource', 's_vs_n');
addParameter(p, 'words', []);
addParameter(p, 'S_condition', 'Sentences');
addParameter(p, 'N_condition', '');
addParameter(p, 'topPct', 0.10);
addParameter(p, 'splitSeed', 42);
parse(p, varargin{:});
ops = p.Results;

avail = unique(sn_obj.condition, 'stable');
ops.S_condition = resolve_condition_name(avail, ...
    {ops.S_condition, 'Sentences', 'SENTENCES', 'sentence'}, 'S');
if isempty(ops.N_condition)
    try
        ops.N_condition = resolve_condition_name(avail, ...
            {'Jabberwocky', 'JABBERWOCKY', 'jabberwocky'}, 'N');
    catch
        ops.N_condition = resolve_condition_name(avail, ...
            {'Nonword-lists', 'Nonwords', 'NONWORDS', 'non-word lists', 'nonwords'}, 'N');
    end
else
    ops.N_condition = resolve_condition_name(avail, ...
        {ops.N_condition, 'Jabberwocky', 'JABBERWOCKY', ...
        'Nonword-lists', 'Nonwords', 'NONWORDS', 'non-word lists'}, 'N');
end

if isempty(ops.words)
    ops.words = [];
end
ops.words = resolve_words_for_sn_obj(sn_obj, ops.words, ops.S_condition);

sigMask = get_lang_sig_mask(sn_obj, ops.signalType, ops.sigSource);
sigMask = sigMask(:);

if strcmpi(ops.signalType, 'bipolar')
    if ~isprop(sn_obj, 'bip_ch_label') || isempty(sn_obj.bip_ch_label)
        error('bip_ch_label required for bipolar ranking.');
    end
    labels = sn_obj.bip_ch_label(:);
    dataField = 'bip_elec_data';
else
    labels = sn_obj.elec_ch_label(:);
    dataField = 'elec_data';
end

nChan = numel(labels);
sigIdx = find(sigMask);
if isempty(sigIdx)
    ranking = empty_ranking_struct(sn_obj, ops, labels);
    return;
end

corrVec = nan(nChan, 1);
if istable(sn_obj.s_vs_n_corr) && ismember(dataField, sn_obj.s_vs_n_corr.Properties.VariableNames)
    raw = sn_obj.s_vs_n_corr.(dataField){1, 1};
    corrVec(1:min(nChan, numel(raw))) = raw(1:min(nChan, numel(raw)));
end

heldoutVec = nan(nChan, 1);
if isprop(sn_obj, 'hg_sn_corr') && ~isempty(sn_obj.hg_sn_corr) ...
        && isfield(sn_obj.hg_sn_corr, 'results')
    if strcmpi(ops.signalType, 'bipolar') ...
            && isfield(sn_obj.hg_sn_corr.results, 'bipolar')
        heldoutVec = expand_to_n(sn_obj.hg_sn_corr.results.bipolar.corr, nChan);
    elseif isfield(sn_obj.hg_sn_corr.results, 'unipolar')
        heldoutVec = expand_to_n(sn_obj.hg_sn_corr.results.unipolar.corr, nChan);
    end
end

pRatio = nan(nChan, 1);
if istable(sn_obj.s_vs_n_p_ratio) && ismember(dataField, sn_obj.s_vs_n_p_ratio.Properties.VariableNames)
    raw = sn_obj.s_vs_n_p_ratio.(dataField){1, 1};
    pRatio(1:min(nChan, numel(raw))) = raw(1:min(nChan, numel(raw)));
end

reliabVec = compute_split_half_reliability(sn_obj, ops, dataField, nChan);

effectVec = heldoutVec;
if sum(isfinite(effectVec(sigIdx))) < 2
    warning('compute_lang_elec_ranking:NoHeldout', ...
        'heldout_corr unavailable for %s; using s_vs_n_corr for composite ranking.', ...
        get_subject_id(sn_obj));
    effectVec = corrVec;
end

zEffect = zscore_within(sigIdx, effectVec);
zRel = zscore_within(sigIdx, reliabVec);
composite = nan(nChan, 1);
composite(sigIdx) = (zEffect(sigIdx) + zRel(sigIdx)) / 2;

[~, sortOrder] = sort(composite(sigIdx), 'descend', 'MissingPlacement', 'last');
rankAmongSig = nan(nChan, 1);
rankAmongSig(sigIdx(sortOrder)) = 1:numel(sigIdx);

nTop = max(1, ceil(numel(sigIdx) * ops.topPct));
topMask = false(nChan, 1);
topMask(sigIdx(sortOrder(1:nTop))) = true;

ranking = struct();
ranking.subject = get_subject_id(sn_obj);
ranking.signalType = ops.signalType;
ranking.sigSource = ops.sigSource;
ranking.S_condition = ops.S_condition;
ranking.N_condition = ops.N_condition;
ranking.words = ops.words;
ranking.topPct = ops.topPct;
ranking.n_sig = numel(sigIdx);
ranking.n_top = nTop;
ranking.labels = labels;
ranking.sig_mask = sigMask;
ranking.top_mask = topMask;
ranking.chan_idx = (1:nChan)';
ranking.s_vs_n_corr = corrVec;
ranking.heldout_corr = heldoutVec;
ranking.split_half_reliab = reliabVec;
ranking.p_ratio = pRatio;
ranking.composite_score = composite;
ranking.rank_among_sig = rankAmongSig;
ranking.top_chan_idx = sigIdx(sortOrder(1:nTop));
end


function ranking = empty_ranking_struct(sn_obj, ops, labels)
ranking = struct();
ranking.subject = get_subject_id(sn_obj);
ranking.signalType = ops.signalType;
ranking.sigSource = ops.sigSource;
ranking.S_condition = ops.S_condition;
ranking.N_condition = ops.N_condition;
ranking.words = ops.words;
ranking.topPct = ops.topPct;
ranking.n_sig = 0;
ranking.n_top = 0;
ranking.labels = labels;
ranking.sig_mask = false(numel(labels), 1);
ranking.top_mask = false(numel(labels), 1);
ranking.chan_idx = (1:numel(labels))';
ranking.s_vs_n_corr = nan(numel(labels), 1);
ranking.heldout_corr = nan(numel(labels), 1);
ranking.split_half_reliab = nan(numel(labels), 1);
ranking.p_ratio = nan(numel(labels), 1);
ranking.composite_score = nan(numel(labels), 1);
ranking.rank_among_sig = nan(numel(labels), 1);
ranking.top_chan_idx = [];
end


function reliabVec = compute_split_half_reliability(sn_obj, ops, dataField, nChan)
reliabVec = nan(nChan, 1);

try
    [~, S_tbl] = sn_obj.get_ave_cond_trial('words', ops.words, ...
        'condition', ops.S_condition);
    [~, N_tbl] = sn_obj.get_ave_cond_trial('words', ops.words, ...
        'condition', ops.N_condition);
catch ME
    warning('compute_lang_elec_ranking:SplitHalfData', ...
        'Split-half reliability skipped for %s: %s', get_subject_id(sn_obj), ME.message);
    return;
end

if ~ismember(dataField, S_tbl.Properties.VariableNames) ...
        || isempty(S_tbl.(dataField){1})
    return;
end

S_data = S_tbl.(dataField){1};  % [nElec x nTrials x nWords]
N_data = N_tbl.(dataField){1};
nElec = size(S_data, 1);

rng(ops.splitSeed);
for e = 1:min(nElec, nChan)
    rS = split_half_corr_trials(S_data(e, :, :));
    rN = split_half_corr_trials(N_data(e, :, :));
    reliabVec(e) = mean([rS, rN], 'omitnan');
end
end


function r = split_half_corr_trials(data3d)
% data3d: [1 x nTrials x nWords]
trialMeans = squeeze(mean(data3d, 3));
if isrow(trialMeans)
    trialMeans = trialMeans(:);
end
n = numel(trialMeans);
if n < 4
    r = NaN;
    return;
end
perm = randperm(n);
nHalf = floor(n / 2);
h1 = perm(1:nHalf);
h2 = perm(nHalf + 1:2 * nHalf);
if numel(h1) < 2 || numel(h2) < 2
    r = NaN;
    return;
end
r = corr(trialMeans(h1), trialMeans(h2), 'Type', 'Spearman', 'Rows', 'pairwise');
end


function z = zscore_within(idx, x)
z = nan(size(x));
vals = x(idx);
if sum(isfinite(vals)) < 2
    z(idx) = 0;
    return;
end
mu = mean(vals, 'omitnan');
sd = std(vals, 0, 'omitnan');
if sd < eps
    z(idx) = 0;
else
    z(idx) = (vals - mu) / sd;
end
end


function out = expand_to_n(vec, n)
out = nan(n, 1);
out(1:min(n, numel(vec))) = vec(1:min(n, numel(vec)));
end


function words = resolve_words_for_sn_obj(sn_obj, words, S_condition)
% Prefer word indices from test_s_vs_n; clamp to available data positions.
if isempty(words) && isprop(sn_obj, 's_vs_n_ops') && ~isempty(sn_obj.s_vs_n_ops) ...
        && isstruct(sn_obj.s_vs_n_ops) && isfield(sn_obj.s_vs_n_ops, 'words') ...
        && ~isempty(sn_obj.s_vs_n_ops.words)
    words = sn_obj.s_vs_n_ops.words;
end

maxW = probe_max_word_positions(sn_obj, S_condition);
if maxW < 1
    wordsFallback = default_words_from_timing(sn_obj);
    maxW = numel(wordsFallback);
end

if isempty(words)
    words = 1:maxW;
else
    words = words(words >= 1 & words <= maxW);
    if isempty(words)
        words = 1:maxW;
    end
end
end


function maxW = probe_max_word_positions(sn_obj, condition)
% Find the largest word index that get_ave_cond_trial can index without error.
maxW = 0;
for w = 1:24
    try
        sn_obj.get_ave_cond_trial('words', w, 'condition', condition);
        maxW = w;
    catch
        break;
    end
end
end


function words = default_words_from_timing(sn_obj)
maxW = 0;
if isprop(sn_obj, 'trial_timing') && ~isempty(sn_obj.trial_timing)
    for i = 1:numel(sn_obj.trial_timing)
        tt = sn_obj.trial_timing{i};
        if istable(tt) && ismember('end', tt.Properties.VariableNames)
            maxW = max(maxW, trial_timing_max_word_index(tt));
        elseif isstruct(tt) && isfield(tt, 'end')
            maxW = max(maxW, numel(tt.end));
        end
    end
end
if maxW < 1
    maxW = 12;
end
words = 1:maxW;
end


function id = get_subject_id(sn_obj)
id = 'unknown';
if isprop(sn_obj, 'subject') && ~isempty(sn_obj.subject)
    id = sn_obj.subject;
end
end
