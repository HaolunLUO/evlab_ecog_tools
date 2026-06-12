function output_xarray_minimal(obj, filename)
    % OUTPUT_XARRAY_MINIMAL Save selected ecog_data properties as HDF5 for Python xarray

    if nargin < 2 || isempty(filename)
        filename = [obj.subject '_' obj.experiment '_ecog_data.h5'];
    end
    if ~strcmp(obj.crunched_file_path(end), filesep)
        save_path = [obj.crunched_file_path filesep filename];
    else
        save_path = [obj.crunched_file_path filename];
    end
    if exist(save_path, 'file')
        delete(save_path);
    end

    % List of properties to save
    keep_props = { ...
        'elec_data', 'bip_elec_data', 'stitch_index', 'sample_freq', 'trial_data', 'stats', ...
        'subject', 'experiment', 'trial_timing', 'events_table', 'condition', 'session', ...
        'elec_ch', 'elec_ch_label', 'elec_ch_prelim_deselect', 'elec_ch_with_IED', ...
        'elec_ch_with_noise', 'elec_ch_user_deselect', 'elec_ch_clean', 'elec_ch_valid', ...
        'elec_ch_type', 'bip_ch', 'bip_ch_label', 'bip_ch_valid', 'bip_ch_grp', 'bip_ch_label_grp' ...
    };

    for i = 1:numel(keep_props)
        pname = keep_props{i};
        if ~isprop(obj, pname), continue; end
        value = obj.(pname);
        if isempty(value), continue; end
        try
            save_value(save_path, ['/' pname], value);
        catch ME
            warning('Failed to save %s: %s', pname, ME.message);
        end
    end

    h5writeatt(save_path, '/', 'class_description', 'ecog_data minimal export for Python xarray');
    fprintf('Saved selected ecog_data fields as xarray-compatible HDF5: %s\n', save_path);

    function save_value(file_path, h5_path, value)
        if isnumeric(value)
            h5create(file_path, h5_path, size(value), 'Datatype', 'double');
            h5write(file_path, h5_path, double(value));
        elseif islogical(value)
            h5create(file_path, h5_path, size(value), 'Datatype', 'uint8');
            h5write(file_path, h5_path, uint8(value));
        elseif ischar(value)
            if isrow(value)
                h5writeatt(file_path, '/', h5_path(2:end), value);
            else
                h5create(file_path, h5_path, size(value), 'Datatype', 'char');
                h5write(file_path, h5_path, value);
            end
        elseif isstring(value) || iscellstr(value) || (iscell(value) && all(cellfun(@ischar,value)))
            char_data = char(value);
            h5create(file_path, h5_path, size(char_data), 'Datatype', 'char');
            h5write(file_path, h5_path, char_data);
        elseif isstruct(value)
            save_struct(file_path, h5_path, value);
        elseif istable(value)
            save_struct(file_path, h5_path, table2struct(value));
        else
            warning('Skipping unsupported property: %s', h5_path);
        end
    end

    function save_struct(file_path, h5_path, S)
        for idx = 1:numel(S)
            item_path = h5_path;
            if numel(S) > 1
                item_path = [h5_path '/item_' num2str(idx)];
            end
            fields = fieldnames(S(idx));
            for f = 1:numel(fields)
                field_name = fields{f};
                field_value = S(idx).(field_name);
                field_path = [item_path '/' field_name];
                save_value(file_path, field_path, field_value);
            end
        end
    end
end
