function output_data_structures(obj)
    % OUTPUT_DATA_STRUCTURES Save ecog_data object to HDF5 file
    %
    
   
    filename = [obj.subject '_' obj.experiment '_ecog_data.h5']
    
    
    % Ensure crunched_file_path ends with file separator
    if ~strcmp(obj.crunched_file_path(end), filesep)
        save_path = [obj.crunched_file_path filesep filename]
    else
        save_path = [obj.crunched_file_path filename]
    end

    % % Delete existing file if present
    % if exist(save_path, 'file')
    %     delete(save_path);
    % end
    
    % Create file using low-level HDF5 functions
    file_id = H5F.create(save_path, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    
    
    % Define group structure
    groups = {
        '/DATA'
        '/INFO'
        '/LABELS'
        '/ANATOMY'
    };
    
    % Create groups
    for i = 1:numel(groups)
        group_id = H5G.create(file_id, groups{i}, 'H5P_DEFAULT');
        H5G.close(group_id);
    end
    
    % Define property mappings
    property_map.DATA = {
        'elec_data', 'bip_elec_data', 'stitch_index', ...
        'sample_freq', 'for_preproc', 'trial_data', 'stats'
    };
    
    property_map.INFO = {
        'subject', 'experiment', 'trial_timing', ...
        'events_table', 'condition', 'session', ...
        'crunched_file_name', 'crunched_file_path', ...
        'raw_file_name', 'raw_file_path'
    };
    
    property_map.LABELS = {
        'elec_ch', 'elec_ch_label', 'elec_ch_prelim_deselect', ...
        'elec_ch_with_IED', 'elec_ch_with_noise', 'elec_ch_user_deselect', ...
        'elec_ch_clean', 'elec_ch_valid', 'elec_ch_type', ...
        'bip_ch', 'bip_ch_label', 'bip_ch_valid', ...
        'bip_ch_grp', 'bip_ch_label_grp'
    };
    
    property_map.ANATOMY = {'anatomy'};
    
    % Save properties to HDF5
    group_names = fieldnames(property_map);
    for g = 1:numel(group_names)
        group_name = group_names{g};
        group_path = ['/' group_name];
        properties = property_map.(group_name);
        
        for p = 1:numel(properties)
            prop_name = properties{p};
            value = obj.(prop_name);
            dataset_path = [group_path '/' prop_name];
            
            % Handle different data types
            if isempty(value)
                % Skip empty properties
                continue;
            elseif isstruct(value)
                save_struct_h5(file_id, dataset_path, value);
            elseif ischar(value)
                save_string_h5(file_id, dataset_path, value);
            elseif iscell(value) && all(cellfun(@ischar, value))
                save_cellstring_h5(file_id, dataset_path, value);
            elseif isnumeric(value) || islogical(value)
                save_numeric_h5(file_id, dataset_path, value);
            else
                warning('Skipping unsupported property: %s', prop_name);
            end
        end
    end
    
    % Add class description attribute
    desc = 'ecog_data class instance';
    space_id = H5S.create('H5S_SCALAR');
    type_id = H5T.copy('H5T_C_S1');
    H5T.set_size(type_id, numel(desc));
    attr_id = H5A.create(file_id, 'class_description', type_id, space_id, 'H5P_DEFAULT');
    H5A.write(attr_id, type_id, desc);
    
    % Close resources
    H5A.close(attr_id);
    H5T.close(type_id);
    H5S.close(space_id);
    H5F.close(file_id);
    
   fprintf('Saved ecog_data to: %s\n', save_path);
    
   
end

 % Nested helper functions
    


