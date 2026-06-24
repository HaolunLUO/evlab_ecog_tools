function [obj, changed] = ensure_fixation_in_obj_trial_timing(obj, fixationDurationSec)
% ENSURE_FIXATION_IN_OBJ_TRIAL_TIMING  Add fix row to each trial if missing.
changed = false;
if nargin < 2 || isempty(fixationDurationSec)
    fixationDurationSec = 0.5;
end
if isempty(obj.trial_timing)
    return;
end
fs = obj.sample_freq;
if isempty(fs) || fs <= 0
    fs = 500;
end
for i = 1:numel(obj.trial_timing)
    tbl = obj.trial_timing{i};
    if ~istable(tbl)
        continue;
    end
    if trial_timing_has_fixation(tbl)
        continue;
    end
    obj.trial_timing{i} = prepend_fixation_row(tbl, fixationDurationSec, fs);
    changed = true;
end
if ~changed
    return;
end
fprintf('Prepended fixation row to trial_timing (%d trials, %.2f s pre-word_1).\n', ...
    numel(obj.trial_timing), fixationDurationSec);
obj.trial_data = [];
if isfield(obj, 'for_preproc') && isstruct(obj.for_preproc) ...
        && isfield(obj.for_preproc, 'trial_timing_raw')
    rawFs = obj.for_preproc.sample_freq_raw;
    if isempty(rawFs); rawFs = fs; end
    for i = 1:min(numel(obj.trial_timing), numel(obj.for_preproc.trial_timing_raw))
        if istable(obj.for_preproc.trial_timing_raw{i})
            obj.for_preproc.trial_timing_raw{i} = prepend_fixation_row( ...
                obj.for_preproc.trial_timing_raw{i}, fixationDurationSec, rawFs);
        end
    end
end
if isfield(obj, 'for_preproc') && isstruct(obj.for_preproc) ...
        && isfield(obj.for_preproc, 'trial_timing_dec')
    decFs = obj.for_preproc.decimation_freq;
    if isempty(decFs); decFs = fs; end
    for i = 1:min(numel(obj.trial_timing), numel(obj.for_preproc.trial_timing_dec))
        if istable(obj.for_preproc.trial_timing_dec{i})
            obj.for_preproc.trial_timing_dec{i} = prepend_fixation_row( ...
                obj.for_preproc.trial_timing_dec{i}, fixationDurationSec, decFs);
        end
    end
end
end
