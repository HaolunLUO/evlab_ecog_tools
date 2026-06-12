
    function save_numeric_h5(loc_id, path, data)
        dims = size(data);
        rank = numel(dims);
        dims_flip = fliplr(dims); % HDF5 uses C-order (row-major)
        space_id = H5S.create_simple(rank, dims_flip, []);
        type_id = get_h5_datatype(data);
        dset_id = H5D.create(loc_id, path, type_id, space_id, 'H5P_DEFAULT');
        H5D.write(dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
        H5D.close(dset_id);
        H5S.close(space_id);
        H5T.close(type_id);
    end