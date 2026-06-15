function obj = sync_obj_from_sn_obj(obj, sn_obj)
% SYNC_OBJ_FROM_SN_OBJ  Copy processed data and stats from sn_obj onto obj.
%
%   obj = sync_obj_from_sn_obj(obj, sn_obj)
%
%   After Step 7.5+, all preprocessing and engine stats live on sn_obj.
%   This keeps the saved `obj` variable in the crunched .mat aligned with
%   the analysis object.

obj.elec_data     = sn_obj.elec_data;
obj.bip_elec_data = sn_obj.bip_elec_data;
obj.stitch_index  = sn_obj.stitch_index;
obj.sample_freq   = sn_obj.sample_freq;
obj.trial_data    = sn_obj.trial_data;

obj.trial_timing = sn_obj.trial_timing;
obj.condition    = sn_obj.condition;
obj.session      = sn_obj.session;

if isprop(sn_obj, 'events_table')
    obj.events_table = sn_obj.events_table;
end
if isprop(sn_obj, 'for_preproc')
    obj.for_preproc = sn_obj.for_preproc;
end

obj.elec_ch_with_IED      = sn_obj.elec_ch_with_IED;
obj.elec_ch_with_noise    = sn_obj.elec_ch_with_noise;
obj.elec_ch_user_deselect = sn_obj.elec_ch_user_deselect;
obj.elec_ch_clean         = sn_obj.elec_ch_clean;
obj.elec_ch_valid         = sn_obj.elec_ch_valid;

obj.bip_ch           = sn_obj.bip_ch;
obj.bip_ch_label     = sn_obj.bip_ch_label;
obj.bip_ch_valid     = sn_obj.bip_ch_valid;
obj.bip_ch_grp       = sn_obj.bip_ch_grp;
obj.bip_ch_label_grp = sn_obj.bip_ch_label_grp;

if isprop(sn_obj, 'anatomy') && ~isempty(sn_obj.anatomy)
    obj.anatomy = sn_obj.anatomy;
end

obj.stats = sn_obj.stats;
end
