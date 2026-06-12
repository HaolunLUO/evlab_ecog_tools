
    function save_string_h5(loc_id, path, str)
        type_id = H5T.copy('H5T_C_S1');
        H5T.set_size(type_id, numel(str));
        space_id = H5S.create('H5S_SCALAR');
        dset_id = H5D.create(loc_id, path, type_id, space_id, 'H5P_DEFAULT');
        H5D.write(dset_id, type_id, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);
        H5D.close(dset_id);
        H5S.close(space_id);
        H5T.close(type_id);
    end
