function results = artifact_svN_crosscheck(sn_obj, varargin)
% ARTIFACT_SVN_CROSSCHECK  Compare sharp-artifact stats to S-vs-N significance.
%
%   results = artifact_svN_crosscheck(sn_obj)
%   results = artifact_svN_crosscheck(sn_obj, 'montage', 'both')
%
%   Inputs
%     sn_obj   - ecog_sn_data_seeg after complete_mit_pipeline_brainstorm (Step 7.5+)
%
%   Name-value options
%     montage               - 'unipolar' (default), 'bipolar', or 'both'
%     runDetectionIfMissing - true (default): call detect_sharp_artifacts if stats absent
%     min_amplitude         - passed to detect_sharp_artifacts (default 15)
%     min_slope             - passed to detect_sharp_artifacts (default 10)
%
%   Output struct fields
%     .unipolar / .bipolar  - per-montage summary + table
%     .subject, .taskType   - copied from sn_obj when available

p = inputParser();
addParameter(p, 'montage', 'unipolar');
addParameter(p, 'runDetectionIfMissing', true);
addParameter(p, 'min_amplitude', 15);
addParameter(p, 'min_slope', 10);
parse(p, varargin{:});
ops = p.Results;

montage = lower(ops.montage);
if ~ismember(montage, {'unipolar', 'bipolar', 'both'})
    error('montage must be ''unipolar'', ''bipolar'', or ''both''.');
end

sn_obj = ensure_artifact_stats(sn_obj, ops);

results = struct();
if isprop(sn_obj, 'subject')
    results.subject = sn_obj.subject;
else
    results.subject = '';
end
if isprop(sn_obj, 'experiment')
    results.taskType = sn_obj.experiment;
else
    results.taskType = '';
end

if ismember(montage, {'unipolar', 'both'})
    results.unipolar = crosscheck_one_montage(sn_obj, 'unipolar');
end

if ismember(montage, {'bipolar', 'both'})
    if has_bipolar_montage(sn_obj)
        results.bipolar = crosscheck_one_montage(sn_obj, 'bipolar');
    else
        results.bipolar = empty_montage_result('No bipolar data on sn_obj.');
    end
end

print_crosscheck_summary(results, montage);
end


function sn_obj = ensure_artifact_stats(sn_obj, ops)
if isprop(sn_obj, 'stats') && isstruct(sn_obj.stats) ...
        && isfield(sn_obj.stats, 'artifact_stats_unipolar') ...
        && ~isempty(sn_obj.stats.artifact_stats_unipolar)
    return;
end

if ~ops.runDetectionIfMissing
    error(['sn_obj.stats.artifact_stats_unipolar is missing. ' ...
        'Re-run the pipeline with detectSharpArtifacts = true, or set ' ...
        'runDetectionIfMissing = true.']);
end

fprintf('artifact_stats missing — running detect_sharp_artifacts ...\n');
sn_obj.detect_sharp_artifacts('min_amplitude', ops.min_amplitude, ...
    'min_slope', ops.min_slope);
end


function out = crosscheck_one_montage(sn_obj, montage)
montage = lower(montage);

if strcmp(montage, 'unipolar')
    sigMask = get_sig_mask(sn_obj, 'elec_data');
    artStats = sn_obj.stats.artifact_stats_unipolar;
    chNums = sn_obj.elec_ch(:);
    chLabels = sn_obj.elec_ch_label(:);
    validMask = logical(sn_obj.elec_ch_valid(:));
    pRatio = get_p_ratio(sn_obj, 'elec_data', numel(chNums));
else
    sigMask = get_sig_mask(sn_obj, 'bip_elec_data');
    artStats = sn_obj.stats.artifact_stats_bipolar;
    chNums = sn_obj.bip_ch(:);
    chLabels = sn_obj.bip_ch_label(:);
    validMask = logical(sn_obj.bip_ch_valid(:));
    pRatio = get_p_ratio(sn_obj, 'bip_elec_data', numel(chNums));
end

nCh = numel(chNums);
artCount = artifact_counts_from_stats(artStats, nCh);
flaggedMask = artCount > 0;

overlapIdx = find(sigMask & flaggedMask);
sigOnlyIdx = find(sigMask & ~flaggedMask);
artifactOnlyIdx = find(~sigMask & flaggedMask);

T = table( ...
    chNums, ...
    chLabels, ...
    sigMask, ...
    flaggedMask, ...
    artCount, ...
    pRatio, ...
    validMask, ...
    'VariableNames', {'chan', 'label', 's_vs_n_sig', 'has_artifacts', ...
    'n_artifacts', 'p_ratio', 'valid'});

out = struct();
out.montage = montage;
out.table = T;
out.overlap = T(sigMask & flaggedMask, :);
out.sig_only = T(sigMask & ~flaggedMask, :);
out.artifact_only = T(~sigMask & flaggedMask, :);
out.overlap_idx = overlapIdx(:);
out.sig_only_idx = sigOnlyIdx(:);
out.artifact_only_idx = artifactOnlyIdx(:);
out.n_channels = nCh;
out.n_sig = sum(sigMask);
out.n_flagged = sum(flaggedMask);
out.n_overlap = numel(overlapIdx);
out.n_sig_only = numel(sigOnlyIdx);
out.n_artifact_only = numel(artifactOnlyIdx);
end


function sigMask = get_sig_mask(sn_obj, fieldName)
if ~istable(sn_obj.s_vs_n_sig) ...
        || ~ismember(fieldName, sn_obj.s_vs_n_sig.Properties.VariableNames)
    error('sn_obj.s_vs_n_sig.%s is missing. Run test_s_vs_n or apply_roi first.', fieldName);
end

raw = sn_obj.s_vs_n_sig.(fieldName){1};
sigMask = logical(raw(:));
end


function pRatio = get_p_ratio(sn_obj, fieldName, nExpected)
pRatio = nan(nExpected, 1);
if ~istable(sn_obj.s_vs_n_p_ratio) ...
        || ~ismember(fieldName, sn_obj.s_vs_n_p_ratio.Properties.VariableNames)
    return;
end

raw = sn_obj.s_vs_n_p_ratio.(fieldName){1};
raw = raw(:);
n = min(numel(raw), nExpected);
pRatio(1:n) = raw(1:n);
end


function artCount = artifact_counts_from_stats(artStats, nChannels)
artCount = zeros(nChannels, 1);
if isempty(artStats)
    return;
end

for k = 1:numel(artStats)
    idx = artStats(k).channel_idx;
    if isempty(idx) || idx < 1 || idx > nChannels
        continue;
    end
    artCount(idx) = artStats(k).artifact_count;
end
end


function tf = has_bipolar_montage(sn_obj)
tf = ~isempty(sn_obj.bip_elec_data) ...
    && isprop(sn_obj, 'stats') && isstruct(sn_obj.stats) ...
    && isfield(sn_obj.stats, 'artifact_stats_bipolar') ...
    && ~isempty(sn_obj.stats.artifact_stats_bipolar);
end


function out = empty_montage_result(msg)
out = struct('montage', 'bipolar', 'message', msg, ...
    'table', table(), 'overlap', table(), 'sig_only', table(), ...
    'artifact_only', table(), 'overlap_idx', [], 'sig_only_idx', [], ...
    'artifact_only_idx', [], 'n_channels', 0, 'n_sig', 0, ...
    'n_flagged', 0, 'n_overlap', 0, 'n_sig_only', 0, 'n_artifact_only', 0);
end


function print_crosscheck_summary(results, montage)
fprintf('\n=== artifact vs S-vs-N cross-check ===\n');
if ~isempty(results.subject)
    fprintf('Subject: %s\n', results.subject);
end
if ~isempty(results.taskType)
    fprintf('Task:    %s\n', results.taskType);
end

if ismember(montage, {'unipolar', 'both'}) && isfield(results, 'unipolar')
    print_one_summary('Unipolar', results.unipolar);
end
if ismember(montage, {'bipolar', 'both'}) && isfield(results, 'bipolar')
    if isfield(results.bipolar, 'message')
        fprintf('\nBipolar: %s\n', results.bipolar.message);
    else
        print_one_summary('Bipolar', results.bipolar);
    end
end
fprintf('======================================\n\n');
end


function print_one_summary(label, block)
fprintf('\n%s\n', label);
fprintf('  Channels:              %d\n', block.n_channels);
fprintf('  S-vs-N significant:    %d\n', block.n_sig);
fprintf('  Sharp-artifact flagged:%d\n', block.n_flagged);
fprintf('  Overlap (sig + art):   %d\n', block.n_overlap);
fprintf('  Sig only:              %d\n', block.n_sig_only);
fprintf('  Artifacts only:        %d\n', block.n_artifact_only);

if block.n_overlap > 0
    fprintf('  Overlap labels:\n');
    for i = 1:height(block.overlap)
        fprintf('    %s  (nArtifacts=%d, pRatio=%.4f)\n', ...
            block.overlap.label{i}, block.overlap.n_artifacts(i), ...
            block.overlap.p_ratio(i));
    end
end
end
