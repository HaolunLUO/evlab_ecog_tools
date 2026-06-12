
    function type_id = get_h5_datatype(data)
        if isa(data, 'single')
            type_id = H5T.copy('H5T_NATIVE_FLOAT');
        elseif isa(data, 'double')
            type_id = H5T.copy('H5T_NATIVE_DOUBLE');
        elseif isa(data, 'int8')
            type_id = H5T.copy('H5T_NATIVE_INT8');
        elseif isa(data, 'uint8')
            type_id = H5T.copy('H5T_NATIVE_UINT8');
        elseif isa(data, 'int16')
            type_id = H5T.copy('H5T_NATIVE_INT16');
        elseif isa(data, 'uint16')
            type_id = H5T.copy('H5T_NATIVE_UINT16');
        elseif isa(data, 'int32')
            type_id = H5T.copy('H5T_NATIVE_INT32');
        elseif isa(data, 'uint32')
            type_id = H5T.copy('H5T_NATIVE_UINT32');
        elseif isa(data, 'int64')
            type_id = H5T.copy('H5T_NATIVE_INT64');
        elseif isa(data, 'uint64')
            type_id = H5T.copy('H5T_NATIVE_UINT64');
        elseif islogical(data)
            type_id = H5T.copy('H5T_NATIVE_UINT8');
        else
            error('Unsupported numeric type: %s', class(data));
        end
    end
