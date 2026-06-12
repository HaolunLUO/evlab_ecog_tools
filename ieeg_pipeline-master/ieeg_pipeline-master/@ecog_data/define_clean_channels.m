function define_clean_channels(obj)
    % Defines the clean channels to be used for analysis.
    %
    % Returns : obj
    
    bad_elecs = union(obj.elec_ch_prelim_deselect,obj.elec_ch_with_IED);
    bad_elecs = union(bad_elecs,obj.elec_ch_with_noise);
    bad_elecs = union(bad_elecs,obj.elec_ch_user_deselect);

    obj.elec_ch_clean = setdiff(obj.elec_ch,bad_elecs,'stable');
    obj.elec_ch_valid = ismember(obj.elec_ch,obj.elec_ch_clean);
end