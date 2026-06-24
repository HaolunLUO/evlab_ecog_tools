function profiles = build_lang_channel_profiles(sn_obj, varargin)
% BUILD_LANG_CHANNEL_PROFILES  Extract lang-responsive channel profiles + MNI.
%
%   profiles = build_lang_channel_profiles(sn_obj)
%   profiles = build_lang_channel_profiles(sn_obj, Name, Value, ...)
%
%   Name-Value:
%     cohort          - label string (default 'local')
%     subject         - override subject id
%     signalType      - 'unipolar' (default) or 'bipolar'
%     words           - word indices (default 1:12)
%     S_condition     - default 'Sentences'
%     N_condition     - default '' -> auto (Jabberwocky or Nonword-lists)
%     sigOnly         - logical (default true)
%     sigSource       - 's_vs_n' (trial-averaged S vs N) or 'wordwise'
%     chanIdx         - optional channel index subset (e.g. top-10% selection)
%     requireMni      - keep only channels with MNI coords (default true)
%     use_odd_for_inference - passed to get_word_averages (default false for plots)

p = inputParser();
addParameter(p, 'cohort', 'local');
addParameter(p, 'subject', '');
addParameter(p, 'signalType', 'unipolar');
addParameter(p, 'words', 1:12);
addParameter(p, 'S_condition', 'Sentences');
addParameter(p, 'N_condition', '');
addParameter(p, 'sigOnly', true);
addParameter(p, 'sigSource', 's_vs_n');
addParameter(p, 'use_odd_for_inference', false);
addParameter(p, 'chanIdx', []);  % optional subset of channel indices (overrides sigOnly mask)
addParameter(p, 'requireMni', true);
parse(p, varargin{:});
ops = p.Results;

if isempty(ops.subject)
    if isprop(sn_obj, 'subject') && ~isempty(sn_obj.subject)
        ops.subject = sn_obj.subject;
    else
        ops.subject = 'unknown';
    end
end

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

ops.words = clamp_words_to_trial_timing(sn_obj, ops.words);

if strcmpi(ops.signalType, 'bipolar')
    if ~isprop(sn_obj, 'bip_ch_label') || isempty(sn_obj.bip_ch_label)
        error('bip_ch_label required for bipolar profiles.');
    end
    labels = sn_obj.bip_ch_label(:);
else
    labels = sn_obj.elec_ch_label(:);
end

if ops.sigOnly
    sigMask = get_lang_sig_mask(sn_obj, ops.signalType, ops.sigSource);
else
    sigMask = true(numel(labels), 1);
end

if strcmpi(ops.signalType, 'bipolar') ...
        && ismember('bip_elec_data', sn_obj.s_vs_n_p_ratio.Properties.VariableNames)
    pRatio = sn_obj.s_vs_n_p_ratio.bip_elec_data{1, 1};
elseif ismember('elec_data', sn_obj.s_vs_n_p_ratio.Properties.VariableNames)
    pRatio = sn_obj.s_vs_n_p_ratio.elec_data{1, 1};
else
    pRatio = nan(numel(labels), 1);
end
pRatio = pRatio(:);

if ~isempty(ops.chanIdx)
    chanIdx = ops.chanIdx(:);
    chanIdx = chanIdx(chanIdx >= 1 & chanIdx <= numel(labels));
    if ops.sigOnly
        chanIdx = chanIdx(sigMask(chanIdx));
    end
else
    chanIdx = find(sigMask);
    if ~ops.sigOnly
        if strcmpi(ops.signalType, 'bipolar')
            chanIdx = find(sn_obj.bip_ch_valid);
        else
            chanIdx = find(sn_obj.elec_ch_valid);
        end
    end
end

if isempty(chanIdx)
    profiles = struct([]);
    warning('No channels to profile for %s (%s).', ops.subject, ops.cohort);
    return;
end

tcArgs = {'words', ops.words, ...
    'useLangElecs', false, ...
    'allElecs', true};
waArgs = [tcArgs, {'use_odd_for_inference', ops.use_odd_for_inference}];

[S_tc, S_sem] = sn_obj.get_timecourses(tcArgs{:}, ...
    'condition', ops.S_condition, 'signalType', ops.signalType);
[N_tc, N_sem] = sn_obj.get_timecourses(tcArgs{:}, ...
    'condition', ops.N_condition, 'signalType', ops.signalType);

[S_w, S_w_sem] = sn_obj.get_word_averages(waArgs{:}, ...
    'condition', ops.S_condition);
[N_w, N_w_sem] = sn_obj.get_word_averages(waArgs{:}, ...
    'condition', ops.N_condition);

if strcmpi(ops.signalType, 'bipolar') && numel(S_w) >= 2
    sigCol = 2;
else
    sigCol = 1;
end

profiles = repmat(struct(), numel(chanIdx), 1);
for k = 1:numel(chanIdx)
    ci = chanIdx(k);
    [mni, ok] = get_lang_channel_mni_coords(sn_obj, ci, ops.signalType);

    profiles(k).cohort = ops.cohort;
    profiles(k).subject = ops.subject;
    profiles(k).signalType = ops.signalType;
    profiles(k).chanIdx = ci;
    profiles(k).label = labels{ci};
    profiles(k).mni = mni;
    profiles(k).hasMni = ok;
    profiles(k).sig = sigMask(ci);
    profiles(k).sigSource = ops.sigSource;
    if ci <= numel(pRatio)
        profiles(k).p_ratio = pRatio(ci);
    else
        profiles(k).p_ratio = NaN;
    end
    profiles(k).S_condition = ops.S_condition;
    profiles(k).N_condition = ops.N_condition;
    profiles(k).sample_freq = sn_obj.sample_freq;
    profiles(k).words = ops.words;

    if ci <= size(S_tc, 1)
        profiles(k).trial_timecourse_S = S_tc(ci, :);
        profiles(k).trial_timecourse_N = N_tc(ci, :);
        profiles(k).trial_timecourse_S_sem = S_sem(ci, :);
        profiles(k).trial_timecourse_N_sem = N_sem(ci, :);
    else
        profiles(k).trial_timecourse_S = [];
        profiles(k).trial_timecourse_N = [];
    end

    if ci <= size(S_w{sigCol}, 1)
        profiles(k).word_ave_S = S_w{sigCol}(ci, :);
        profiles(k).word_ave_N = N_w{sigCol}(ci, :);
        profiles(k).word_ave_S_sem = S_w_sem{sigCol}(ci, :);
        profiles(k).word_ave_N_sem = N_w_sem{sigCol}(ci, :);
    else
        profiles(k).word_ave_S = [];
        profiles(k).word_ave_N = [];
    end
end

if ops.requireMni
    profiles = profiles([profiles.hasMni]);
    if isempty(profiles)
        warning('No channels with valid MNI for %s (%s).', ops.subject, ops.cohort);
    end
end

end


function words = clamp_words_to_trial_timing(sn_obj, words)
maxW = max_word_index_in_timing(sn_obj);
if isempty(words)
    words = 1:maxW;
    return;
end
words = words(words >= 1 & words <= maxW);
if isempty(words)
    words = 1:maxW;
end
end


function maxW = max_word_index_in_timing(sn_obj)
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
end
