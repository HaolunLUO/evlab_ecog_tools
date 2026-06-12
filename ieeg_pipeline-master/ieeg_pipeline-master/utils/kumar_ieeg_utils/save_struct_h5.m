function save_struct_h5(loc_id, path, S)
        group_id = H5G.create(loc_id, path, 'H5P_DEFAULT');
        fields = fieldnames(S);
        for f = 1:numel(fields)
            field_name = fields{f};
            field_value = S.(field_name);
            field_path = [path '/' field_name];
            
            % Recursively save struct fields
            if isstruct(field_value)
                save_struct_h5(loc_id, field_path, field_value);
            elseif isnumeric(field_value)
                save_numeric_h5(loc_id, field_path, field_value);
            elseif ischar(field_value)
                save_string_h5(loc_id, field_path, field_value);
            end
        end
        H5G.close(group_id);
    end