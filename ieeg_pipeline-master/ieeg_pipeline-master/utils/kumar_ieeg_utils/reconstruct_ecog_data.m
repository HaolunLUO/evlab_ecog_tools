function obj = reconstruct_ecog_data(filename)
% RECONSTRUCT_ECOG_DATA Reconstruct ecog_data object from HDF5 file
%   obj = reconstruct_ecog_data(filename) recreates the ecog_data object
%   from an HDF5 file created with output_xarray()

% Read entire file structure
data = h5load_recursive(filename);

% Create object with required properties
obj = ecog_data(...
    data.for_preproc,...
    data.subject,...
    data.experiment,...
    data.crunched_file_name,...
    data.crunched_file_path,...
    data.raw_file_name,...
    data.raw_file_path,...
    data.elec_ch_label,...
    data.elec_ch,...
    data.elec_ch_prelim_deselect,...
    data.elec_ch_type...
);

% Set all other properties
props = fieldnames(data);
for i = 1:numel(props)
    pname = props{i};
    if isprop(obj, pname) && ~any(strcmp(pname, {'for_preproc', 'subject', 'experiment', ...
            'crunched_file_name', 'crunched_file_path', 'raw_file_name', 'raw_file_path', ...
            'elec_ch_label', 'elec_ch', 'elec_ch_prelim_deselect', 'elec_ch_type'}))
        obj.(pname) = data.(pname);
    end
end

% Nested helper functions
    function s = h5load_recursive(filename, path)
        % Recursively load HDF5 file into struct
        if nargin < 2
            path = '/';
        end
        s = struct();
        info = h5info(filename, path);
        
        % Process datasets
        for i = 1:numel(info.Datasets)
            dset_name = info.Datasets(i).Name;
            dset_path = [path '/' dset_name];
            value = h5read(filename, dset_path);
            
            % Handle special types via attributes
            try
                attrs = h5readatt(filename, dset_path, 'MATLABClass');
                if strcmp(attrs, 'dfilt.df2sos')
                    value = dfilt.df2sos(value);
                elseif strcmp(attrs, 'fdesign')
                    value = fdesign(value);
                end
            catch
            end
            
            s.(dset_name) = value;
        end
        
        % Process groups (structs)
        for i = 1:numel(info.Groups)
            group = info.Groups(i);
            group_path = group.Name;
            [~, name] = fileparts(group_path);
            s.(name) = h5load_recursive(filename, group_path);
        end
        
        % Process attributes (special cases)
        for i = 1:numel(info.Attributes)
            attr = info.Attributes(i);
            if strcmp(attr.Name, 'VariableNames') && contains(path, 'events_table')
                % Reconstruct table from struct and variable names
               % Fix for VariableNames mismatch in reconstruct_ecog_data/h5load_recursive
               varnames = strsplit(attr.Value, '|');
                fields = fieldnames(s);
                common_vars = intersect(varnames, fields, 'stable');
                missing_vars = setdiff(fields, common_vars, 'stable');
                all_vars = [common_vars; missing_vars];
                s = struct2table(s, 'AsArray', true);
                s.Properties.VariableNames = all_vars;


            end
        end
    end
end
