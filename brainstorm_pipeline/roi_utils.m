% ROI_UTILS  Helper functions for cross-task ROI transfer.
%
% These functions are called by complete_mit_pipeline_brainstorm.m.
% They are defined as standalone scripts rather than local functions so
% that they can be reused across pipeline runs without re-running the
% full pipeline script.
%
% Functions in this file:
%   save_roi_from_sn_obj   - save significant channels as an ROI struct
%   apply_roi_to_sn_obj    - apply a saved ROI to a new task's sn_obj
%   normalize_labels       - case-insensitive label matching helper

% ------------------------------------------------------------------
% Nothing runs at the script level; all logic is in the functions.
% ------------------------------------------------------------------


function save_roi_from_sn_obj(sn_obj, roiFile, roiSourceTask)
% SAVE_ROI_FROM_SN_OBJ  Save significant-channel labels as an ROI file.
%
%   save_roi_from_sn_obj(sn_obj, roiFile, roiSourceTask)
%
%   sn_obj        - ecog_sn_data object after test_s_vs_n()
%   roiFile       - full path to save the .mat file
%   roiSourceTask - string label of the source task (for provenance)

roi = struct();
roi.subject    = sn_obj.subject;
roi.sourceTask = roiSourceTask;
roi.created    = datestr(now, 'yyyy-mm-dd HH:MM:SS');

sigUni = find(sn_obj.s_vs_n_sig.elec_data{1});
roi.elec_labels = sn_obj.elec_ch_label(sigUni(:));

roi.bip_labels = {};
if ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames) && ~isempty(sn_obj.bip_elec_data)
    sigBip = find(sn_obj.s_vs_n_sig.bip_elec_data{1});
    roi.bip_labels = sn_obj.bip_ch_label(sigBip(:));
end

roiDir = fileparts(roiFile);
if ~exist(roiDir,'dir'); mkdir(roiDir); end
save(roiFile, 'roi', '-v7.3');
end


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
%
%   Enables robust label matching across systems that may differ in
%   separator style (underscore vs hyphen) or capitalisation.
out = upper(string(labels));
out = regexprep(out, '[^A-Z0-9]', '');
end
