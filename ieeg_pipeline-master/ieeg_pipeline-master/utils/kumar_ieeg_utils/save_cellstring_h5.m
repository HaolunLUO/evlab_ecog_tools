
    function save_cellstring_h5(loc_id, path, cellstr)
        % Convert to 2D char array
        char_data = char(cellstr);
        dims = size(char_data);
        space_id = H5S.create_simple(2, fliplr(dims), []);
        type_id = H5T.copy('H5T_C_S1');
        H5T.set_size(type_id, size(char_data, 2));
        dset_id = H5D.create(loc_id, path, type_id, space_id, 'H5P_DEFAULT');
        H5D.write(dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', char_data);
        H5D.close(dset_id);
        H5S.close(space_id);
        H5T.close(type_id);
    end