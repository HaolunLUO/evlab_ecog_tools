function groupResult = build_group_result_from_sn_obj(sn_obj, subjectName, taskType, taskConfig, opts)
% BUILD_GROUP_RESULT_FROM_SN_OBJ  Pack sn_obj localization outputs for group analysis.
%
%   Mirrors the groupResult .mat written by the legacy MIT_multi_single workflow
%   (sig masks, bipolar labels, MNI coords, held-out effect sizes, timecourses).

if nargin < 5 || isempty(opts)
    opts = struct();
end
opts = fill_default_opts(opts);

if ~istable(sn_obj.s_vs_n_sig)
    error('sn_obj.s_vs_n_sig is missing or invalid.');
end

groupResult = struct();
groupResult.subject  = subjectName;
groupResult.taskType = taskType;
if ~isempty(opts.sourceFile)
    groupResult.sourceFile = opts.sourceFile;
end

groupResult.sig_uni     = logical(sn_obj.s_vs_n_sig.elec_data{1});
groupResult.p_ratio_uni = sn_obj.s_vs_n_p_ratio.elec_data{1};
groupResult.elec_labels = sn_obj.elec_ch_label;
groupResult.elec_valid  = sn_obj.elec_ch_valid;
groupResult.elec_clean  = sn_obj.elec_ch_clean;
groupResult.nClean      = numel(sn_obj.elec_ch_clean);
groupResult.nSig        = sum(groupResult.sig_uni);

if ~isempty(sn_obj.bip_elec_data) && ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
    groupResult.sig_bip     = logical(sn_obj.s_vs_n_sig.bip_elec_data{1});
    groupResult.p_ratio_bip = sn_obj.s_vs_n_p_ratio.bip_elec_data{1};
    groupResult.bip_labels  = sn_obj.bip_ch_label;
    groupResult.nSigBip     = sum(groupResult.sig_bip);
else
    groupResult.sig_bip  = [];
    groupResult.bip_labels = {};
    groupResult.nSigBip  = 0;
end

if isfield(sn_obj, 'anatomy') && isfield(sn_obj.anatomy, 'mni_space') ...
        && isfield(sn_obj.anatomy.mni_space, 'tala')
    allMNI = sn_obj.anatomy.mni_space.tala.electrodes;
    groupResult.mni_coords = allMNI;
    sigIdx = find(groupResult.sig_uni);
    mappedCells = sn_obj.anatomy.mapping(sigIdx);
    ok = ~cellfun(@isempty, mappedCells);
    sigAnatomyIdx = cell2mat(mappedCells(ok));
    sigAnatomyIdx = sigAnatomyIdx(sigAnatomyIdx >= 1 & sigAnatomyIdx <= size(allMNI,1));
    sigMask = false(size(allMNI,1), 1);
    sigMask(sigAnatomyIdx) = true;
    groupResult.mni_sig_mask = sigMask;
end

groupResult.sample_freq = sn_obj.sample_freq;
groupResult.nWords      = numel(taskConfig.testWords);

if opts.includeEffectSizes && isprop(sn_obj, 'hg_power_diff') && ~isempty(sn_obj.hg_power_diff)
    groupResult.hg_power_diff = sn_obj.hg_power_diff.results;
end
if opts.includeEffectSizes && isprop(sn_obj, 'hg_sn_corr') && ~isempty(sn_obj.hg_sn_corr)
    groupResult.hg_sn_corr = sn_obj.hg_sn_corr.results;
end

if opts.includeWordwise && isprop(sn_obj, 'langloc_wordwise') && ~isempty(sn_obj.langloc_wordwise)
    groupResult.langloc_wordwise = sn_obj.langloc_wordwise.results;
    if isfield(sn_obj.langloc_wordwise.results, 'unipolar')
        groupResult.sig_wordwise_uni = sn_obj.langloc_wordwise.results.unipolar.is_sig;
        groupResult.nSigWordwise     = sum(groupResult.sig_wordwise_uni);
    end
end

if opts.includeWordBoundaries && isprop(sn_obj, 's_vs_n_wordboundaries_sigUnipolarChannels') ...
        && ~isempty(sn_obj.s_vs_n_wordboundaries_sigUnipolarChannels)
    groupResult.sig_wordboundaries_uni = sn_obj.s_vs_n_wordboundaries_sigUnipolarChannels;
end

try
    [S_tc, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, ...
        'condition', taskConfig.S_condition, 'signalType', 'unipolar');
    [N_tc, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, ...
        'condition', taskConfig.N_condition, 'signalType', 'unipolar');
    groupResult.S_timecourse_mean = mean(S_tc, 1);
    groupResult.S_timecourse_sem  = std(S_tc, [], 1) / sqrt(size(S_tc,1));
    groupResult.N_timecourse_mean = mean(N_tc, 1);
    groupResult.N_timecourse_sem  = std(N_tc, [], 1) / sqrt(size(N_tc,1));
    groupResult.nSigElecs_tc      = size(S_tc, 1);
catch ME
    warning('Monopolar timecourses failed: %s', ME.message);
end

if ~isempty(sn_obj.bip_elec_data) && ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames) ...
        && sum(sn_obj.s_vs_n_sig.bip_elec_data{1}) > 0
    try
        [S_tc_bip, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, ...
            'condition', taskConfig.S_condition, 'signalType', 'bipolar');
        [N_tc_bip, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, ...
            'condition', taskConfig.N_condition, 'signalType', 'bipolar');
        groupResult.S_timecourse_mean_bip = mean(S_tc_bip, 1);
        groupResult.S_timecourse_sem_bip  = std(S_tc_bip, [], 1) / sqrt(size(S_tc_bip,1));
        groupResult.N_timecourse_mean_bip = mean(N_tc_bip, 1);
        groupResult.N_timecourse_sem_bip  = std(N_tc_bip, [], 1) / sqrt(size(N_tc_bip,1));
        groupResult.nSigElecs_tc_bip      = size(S_tc_bip, 1);
    catch ME
        warning('Bipolar timecourses failed: %s', ME.message);
    end
end
end


function opts = fill_default_opts(opts)
defaults = struct( ...
    'sourceFile', '', ...
    'includeEffectSizes', true, ...
    'includeWordwise', true, ...
    'includeWordBoundaries', false, ...
    'useOddForInference', false);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end
end
