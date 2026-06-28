function save_roi_from_sn_obj(sn_obj, roiFile, roiSourceTask)
% SAVE_ROI_FROM_SN_OBJ  Save significant-channel labels as an ROI file.
%
%   save_roi_from_sn_obj(sn_obj, roiFile, roiSourceTask)
%
%   sn_obj        - ecog_sn_data_seeg object after test_s_vs_n()
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
