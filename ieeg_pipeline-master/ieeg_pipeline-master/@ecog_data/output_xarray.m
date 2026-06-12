function output_xarray(obj, filename)
    % OUTPUT_XARRAY Save ecog_data object as HDF5 for Python xarray
    % Converts all data to HDF5-compatible formats for Python
    
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
    
    % Save object properties directly to HDF5
    props = properties(obj);
    for i = 1:numel(props)
        pname = props{i};
        value = obj.(pname);
        if isempty(value)
            continue;
        end
        try
            save_hdf5_compatible(save_path, ['/' pname], value);
        catch ME
            warning('Failed to save %s: %s', pname, ME.message);
        end
    end
    
    h5writeatt(save_path, '/', 'class_description', 'ecog_data MATLAB export for Python');
    h5writeatt(save_path, '/', 'matlab_class', 'ecog_data');
    fprintf('Saved ecog_data as Python-compatible HDF5: %s\n', save_path);
    
    function save_hdf5_compatible(file_path, h5_path, value)
        % Save values in HDF5-compatible format for Python
        
        if isnumeric(value) || islogical(value)
            % Save numeric/logical as datasets
            if islogical(value)
                value = uint8(value);  % Convert logical to uint8
            end
            create_and_write_dataset(file_path, h5_path, value);
            
        elseif ischar(value)
            % Save char as fixed-length string dataset
            if isrow(value)
                % Single string - save as scalar string dataset
                h5create(file_path, h5_path, [1], 'Datatype', 'string');
                h5write(file_path, h5_path, string(value));
            else
                % Multi-row char matrix - save each row as string array
                strings = cellstr(value);
                h5create(file_path, h5_path, [length(strings)], 'Datatype', 'string');
                h5write(file_path, h5_path, string(strings));
            end
            
        elseif isstring(value)
            % Save string array as string dataset
            h5create(file_path, h5_path, size(value), 'Datatype', 'string');
            h5write(file_path, h5_path, value);
            
        elseif iscellstr(value) || (iscell(value) && all(cellfun(@ischar, value)))
            % Save cellstr as string array
            h5create(file_path, h5_path, [length(value)], 'Datatype', 'string');
            h5write(file_path, h5_path, string(value));
            
        elseif iscell(value)
            % Handle mixed cell arrays by flattening to datasets
            save_cell_array(file_path, h5_path, value);
            
        elseif isstruct(value)
            % Save struct as HDF5 group with subfields
            save_struct_as_group(file_path, h5_path, value);
            
        elseif istable(value)
            % Save table as group with columns as datasets
            save_table_as_group(file_path, h5_path, value);
            
        elseif isa(value, 'dfilt.df2sos')
            % Save filter as coefficients with metadata
            create_and_write_dataset(file_path, h5_path, value.sosMatrix);
            h5writeatt(file_path, h5_path, 'filter_type', 'dfilt_df2sos');
            
        elseif isa(value, 'fdesign')
            % Save fdesign as string specification
            h5create(file_path, h5_path, [1], 'Datatype', 'string');
            h5write(file_path, h5_path, string(value.Specification));
            h5writeatt(file_path, h5_path, 'filter_type', 'fdesign');
            
        else
            % Convert unknown types to string representation
            try
                str_repr = evalc('disp(value)');
                h5create(file_path, h5_path, [1], 'Datatype', 'string');
                h5write(file_path, h5_path, string(str_repr));
                h5writeatt(file_path, h5_path, 'original_type', class(value));
            catch
                warning('Could not save %s of type %s', h5_path, class(value));
            end
        end
    end
    
    function create_and_write_dataset(file_path, h5_path, data)
        % Create and write numeric dataset with chunking for large arrays
        sz = size(data);
        if numel(data) > 1000
            chunk_size = min(sz, 100);
            h5create(file_path, h5_path, sz, 'Datatype', class(data), 'ChunkSize', chunk_size);
        else
            h5create(file_path, h5_path, sz, 'Datatype', class(data));
        end
        h5write(file_path, h5_path, data);
    end
    
    function save_cell_array(file_path, h5_path, cell_value)
        % Save cell array by creating indexed datasets
        h5writeatt(file_path, '/', [h5_path(2:end) '_type'], 'cell_array');
        h5writeatt(file_path, '/', [h5_path(2:end) '_length'], length(cell_value));
        
        % Check if homogeneous numeric data
        if all(cellfun(@isnumeric, cell_value)) && all(cellfun(@isscalar, cell_value))
            % Save as numeric array
            numeric_data = cell2mat(cell_value);
            create_and_write_dataset(file_path, h5_path, numeric_data);
            h5writeatt(file_path, h5_path, 'original_type', 'cell_numeric');
        else
            % Save each cell as separate indexed dataset
            for i = 1:length(cell_value)
                item_path = sprintf('%s_item_%d', h5_path, i);
                save_hdf5_compatible(file_path, item_path, cell_value{i});
            end
        end
    end
    
    function save_struct_as_group(file_path, h5_path, struct_value)
        % Save struct as HDF5 group with fields as datasets
        h5writeatt(file_path, '/', [h5_path(2:end) '_type'], 'struct');
        
        for idx = 1:numel(struct_value)
            if numel(struct_value) > 1
                item_path = sprintf('%s/item_%d', h5_path, idx);
            else
                item_path = h5_path;
            end
            
            fields = fieldnames(struct_value(idx));
            for f = 1:numel(fields)
                field_name = fields{f};
                field_value = struct_value(idx).(field_name);
                field_path = sprintf('%s/%s', item_path, field_name);
                save_hdf5_compatible(file_path, field_path, field_value);
            end
        end
    end
    
    function save_table_as_group(file_path, h5_path, table_value)
        % Save table as group with each column as dataset
        h5writeatt(file_path, '/', [h5_path(2:end) '_type'], 'table');
        
        % Save variable names
        var_names = table_value.Properties.VariableNames;
        h5create(file_path, [h5_path '/variable_names'], [length(var_names)], 'Datatype', 'string');
        h5write(file_path, [h5_path '/variable_names'], string(var_names));
        
        % Save each column
        for i = 1:width(table_value)
            var_name = var_names{i};
            var_data = table_value.(var_name);
            col_path = sprintf('%s/%s', h5_path, var_name);
            save_hdf5_compatible(file_path, col_path, var_data);
        end
        
        % Save table metadata
        h5writeatt(file_path, h5_path, 'n_rows', height(table_value));
        h5writeatt(file_path, h5_path, 'n_cols', width(table_value));
    end
end
