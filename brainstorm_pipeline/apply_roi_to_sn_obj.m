function sn_obj = apply_roi_to_sn_obj(sn_obj, roi)
% APPLY_ROI_TO_SN_OBJ  Apply a saved ROI onto a new task's sn_obj.
%
%   sn_obj = apply_roi_to_sn_obj(sn_obj, roi)
%
%   Matches channel labels (case/punctuation insensitive) between the ROI
%   and the current task.  Populates sn_obj.s_vs_n_sig and
%   sn_obj.s_vs_n_p_ratio in the format expected by lang_resp_plots().

roiUni = normalize_labels(roi.elec_labels);
curUni = normalize_labels(sn_obj.elec_ch_label);

uniMask = ismember(curUni, roiUni);
uniMask = uniMask(:) & logical(sn_obj.elec_ch_valid(:));

bipMask = [];
if ~isempty(sn_obj.bip_elec_data)
    if isfield(roi, 'bip_labels') && ~isempty(roi.bip_labels)
        roiBip = normalize_labels(roi.bip_labels);
        curBip = normalize_labels(sn_obj.bip_ch_label);
        bipMask = ismember(curBip, roiBip);
    else
        % Fallback: mark bipolar pairs that include any ROI unipolar contact
        bipMask = false(numel(sn_obj.bip_ch_label), 1);
        for i = 1:numel(sn_obj.bip_ch_label)
            parts = split(string(sn_obj.bip_ch_label{i}), "-");
            partsNorm = normalize_labels(cellstr(parts));
            bipMask(i) = any(ismember(partsNorm, roiUni));
        end
    end
    bipMask = bipMask(:) & logical(sn_obj.bip_ch_valid(:));
end

if ~isempty(bipMask)
    sn_obj.s_vs_n_sig = table("roi_sig", {uniMask}, {bipMask}, ...
        'VariableNames', {'key','elec_data','bip_elec_data'});
    sn_obj.s_vs_n_p_ratio = table("roi_p_ratio", {nan(size(uniMask))}, {nan(size(bipMask))}, ...
        'VariableNames', {'key','elec_data','bip_elec_data'});
else
    sn_obj.s_vs_n_sig = table("roi_sig", {uniMask}, ...
        'VariableNames', {'key','elec_data'});
    sn_obj.s_vs_n_p_ratio = table("roi_p_ratio", {nan(size(uniMask))}, ...
        'VariableNames', {'key','elec_data'});
end
end


function out = normalize_labels(labels)
% NORMALIZE_LABELS  Uppercase and strip non-alphanumeric characters.
out = upper(string(labels));
out = regexprep(out, '[^A-Z0-9]', '');
end
