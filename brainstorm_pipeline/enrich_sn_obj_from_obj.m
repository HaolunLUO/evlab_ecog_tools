function sn_obj = enrich_sn_obj_from_obj(sn_obj, obj)
% ENRICH_SN_OBJ_FROM_OBJ  Copy MNI / anatomy fields from obj onto sn_obj if missing.

copyNames = {'elec_ch_pos_mni', 'bip_ch_pos_mni', 'filt_ops', 'anatomy'};
for i = 1:numel(copyNames)
    name = copyNames{i};
    if ~has_usable_prop(sn_obj, name) && has_usable_prop(obj, name)
        sn_obj = set_object_prop(sn_obj, name, obj.(name));
    end
end
end


function tf = has_usable_prop(obj, name)
tf = false;
if ~isprop(obj, name)
    return;
end
val = obj.(name);
if iscell(val)
    tf = ~isempty(val);
elseif isnumeric(val)
    tf = ~isempty(val);
elseif isstruct(val)
    tf = ~isempty(fieldnames(val));
else
    tf = ~isempty(val);
end
end


function obj = set_object_prop(obj, name, value)
if isprop(obj, name)
    obj.(name) = value;
else
    obj.addprop(name);
    obj.(name) = value;
end
end
