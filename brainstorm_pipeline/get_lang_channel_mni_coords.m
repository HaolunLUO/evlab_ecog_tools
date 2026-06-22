function [mni, ok] = get_lang_channel_mni_coords(sn_obj, chanIdx, signalType)
% GET_LANG_CHANNEL_MNI_COORDS  MNI XYZ for a channel on sn_obj.
%
%   Unipolar: elec_ch_pos_mni, or VERA anatomy.mapping -> mni_space electrodes.
%   Bipolar:  bip_ch_pos_mni, or midpoint of the two contacts in bip_ch_label
%             (e.g. 'LF12-LF13') using unipolar MNI coords.

if nargin < 3 || isempty(signalType)
    signalType = 'unipolar';
end

mni = [NaN NaN NaN];
ok = false;

if isempty(chanIdx) || chanIdx < 1
    return;
end

if strcmpi(signalType, 'bipolar')
    [mni, ok] = bipolar_mni_coords(sn_obj, chanIdx);
else
    [mni, ok] = unipolar_mni_coords(sn_obj, chanIdx);
end

end


function [mni, ok] = unipolar_mni_coords(sn_obj, uniIdx)
mni = [NaN NaN NaN];
ok = false;

posMni = get_object_mni_field(sn_obj, 'elec_ch_pos_mni');
if ~isempty(posMni)
    [mni, ok] = mni_at_index(posMni, uniIdx);
    if ok; return; end
end

[mni, ok] = mni_from_anatomy_mapping(sn_obj, uniIdx);
if ok; return; end

[mni, ok] = mni_from_clinical_info(sn_obj, uniIdx);
end


function [mni, ok] = bipolar_mni_coords(sn_obj, bipIdx)
mni = [NaN NaN NaN];
ok = false;

posMni = get_object_mni_field(sn_obj, 'bip_ch_pos_mni');
if ~isempty(posMni)
    [mni, ok] = mni_at_index(posMni, bipIdx);
    if ok; return; end
end

labels = get_object_cellstr(sn_obj, 'bip_ch_label');
pairLabels = get_bipolar_pair_labels(sn_obj, bipIdx, labels);
if numel(pairLabels) >= 2
    elecLabels = get_object_cellstr(sn_obj, 'elec_ch_label');
    idx1 = label_to_index(pairLabels{1}, elecLabels);
    idx2 = label_to_index(pairLabels{2}, elecLabels);
    pairIdx = [idx1, idx2];
    pairIdx = pairIdx(pairIdx > 0);
else
    if isempty(labels) || bipIdx > numel(labels)
        [mni, ok] = mni_from_clinical_info_bipolar(sn_obj, pairLabels);
        return;
    end
    pairIdx = bipolar_contact_indices(sn_obj, labels{bipIdx});
end
if numel(pairIdx) < 2
    [mni, ok] = mni_from_clinical_info_bipolar(sn_obj, pairLabels);
    return;
end

[mni1, ok1] = unipolar_mni_coords(sn_obj, pairIdx(1));
[mni2, ok2] = unipolar_mni_coords(sn_obj, pairIdx(2));
if ok1 && ok2
    mni = mean([mni1; mni2], 1);
    ok = true;
elseif ok1
    mni = mni1;
    ok = true;
elseif ok2
    mni = mni2;
    ok = true;
end
end


function pairLabels = get_bipolar_pair_labels(sn_obj, bipIdx, labels)
pairLabels = {};

if ~isempty(labels) && bipIdx <= numel(labels)
    parts = split_bipolar_label(labels{bipIdx});
    if numel(parts) >= 2
        pairLabels = parts(1:2);
        return;
    end
end

if isprop(sn_obj, 'bip_ch_label_valid') && ~isempty(sn_obj.bip_ch_label_valid)
    bvl = sn_obj.bip_ch_label_valid;
    if size(bvl, 2) >= 2 && bipIdx <= size(bvl, 1)
        pairLabels = {bvl{bipIdx, 1}, bvl{bipIdx, 2}};
        return;
    end
end

if ~isempty(labels) && bipIdx <= numel(labels)
    pairLabels = {labels{bipIdx}};
end
end


function pairIdx = bipolar_contact_indices(sn_obj, bipLabel)
pairIdx = [];

parts = split_bipolar_label(bipLabel);
if numel(parts) < 2
    return;
end

elecLabels = get_object_cellstr(sn_obj, 'elec_ch_label');
if isempty(elecLabels)
    return;
end

pairIdx = [label_to_index(parts{1}, elecLabels), ...
           label_to_index(parts{2}, elecLabels)];
pairIdx = pairIdx(pairIdx > 0);
if numel(pairIdx) >= 2
    return;
end

% Fallback: match full bipolar label against unipolar names
hit = find(strcmpi(elecLabels, bipLabel), 1);
if ~isempty(hit)
    pairIdx = hit;
end
end


function idx = label_to_index(label, elecLabels)
idx = 0;
if isempty(label) || isempty(elecLabels)
    return;
end

hit = find(strcmpi(strtrim(elecLabels), strtrim(label)), 1, 'first');
if ~isempty(hit)
    idx = hit;
    return;
end

norm = @(s) regexprep(lower(strtrim(s)), '[^a-z0-9]', '');
target = norm(label);
for i = 1:numel(elecLabels)
    if strcmp(norm(elecLabels{i}), target)
        idx = i;
        return;
    end
end
end


function [mni, ok] = mni_from_anatomy_mapping(sn_obj, chanIdx)
mni = [NaN NaN NaN];
ok = false;

if ~has_anatomy(sn_obj)
    return;
end

mapping = sn_obj.anatomy.mapping;
if numel(mapping) < chanIdx
    return;
end

mapped = mapping{chanIdx};
if isempty(mapped)
    return;
end

try
    coords = anatomy_electrode_coords(sn_obj, 'mni');
catch
    return;
end

idx = mapped(1);
if idx >= 1 && idx <= size(coords, 1)
    mni = coords(idx, :);
    ok = all(isfinite(mni));
end
end


function tf = has_anatomy(sn_obj)
tf = isprop(sn_obj, 'anatomy') && ~isempty(sn_obj.anatomy) && isstruct(sn_obj.anatomy);
end


function arr = get_object_mni_field(obj, name)
arr = [];
if isprop(obj, name)
    val = obj.(name);
    if (isnumeric(val) || iscell(val)) && ~isempty(val)
        arr = val;
    end
end
end


function [row, ok] = mni_at_index(posMni, idx)
row = [NaN NaN NaN];
ok = false;
if isempty(posMni) || idx < 1
    return;
end
if isnumeric(posMni) && idx <= size(posMni, 1)
    row = posMni(idx, 1:min(3, size(posMni, 2)));
    ok = numel(row) == 3 && all(isfinite(row));
elseif iscell(posMni) && idx <= numel(posMni)
    v = posMni{idx};
    if isnumeric(v)
        row = v(1:min(3, numel(v)));
        if numel(row) < 3
            row = [row, NaN(1, 3 - numel(row))];
        end
        ok = all(isfinite(row));
    end
end
end


function [mni, ok] = mni_from_clinical_info(sn_obj, uniIdx)
mni = [NaN NaN NaN];
ok = false;

labels = get_object_cellstr(sn_obj, 'elec_ch_label');
if uniIdx > numel(labels)
    return;
end

[mni, ok] = mni_from_clinical_label(sn_obj, labels{uniIdx});
end


function [mni, ok] = mni_from_clinical_info_bipolar(sn_obj, pairLabels)
mni = [NaN NaN NaN];
ok = false;

if isempty(pairLabels)
    return;
end
if numel(pairLabels) == 1
    parts = split_bipolar_label(pairLabels{1});
    if numel(parts) >= 2
        pairLabels = parts(1:2);
    else
        return;
    end
end
if numel(pairLabels) < 2
    return;
end

[mni1, ok1] = mni_from_clinical_label(sn_obj, pairLabels{1});
[mni2, ok2] = mni_from_clinical_label(sn_obj, pairLabels{2});
if ok1 && ok2
    mni = mean([mni1; mni2], 1);
    ok = true;
elseif ok1
    mni = mni1;
    ok = true;
elseif ok2
    mni = mni2;
    ok = true;
end
end


function [mni, ok] = mni_from_clinical_label(sn_obj, label)
mni = [NaN NaN NaN];
ok = false;

if ~isprop(sn_obj, 'filt_ops') || isempty(sn_obj.filt_ops)
    return;
end
ops = sn_obj.filt_ops;
if ~isstruct(ops) || ~isfield(ops, 'elec_clinic_info') || isempty(ops.elec_clinic_info)
    return;
end

ci = ops.elec_clinic_info;
if ~istable(ci) && ~isstruct(ci)
    return;
end

if istable(ci)
    hit = find_clinical_label_row(ci, label);
    if isempty(hit)
        return;
    end
    row = ci(hit, :);
    if ismember('mni_linear_x', ci.Properties.VariableNames)
        mni = [row.mni_linear_x, row.mni_linear_y, row.mni_linear_z];
    elseif all(ismember({'x', 'y', 'z'}, ci.Properties.VariableNames))
        mni = [row.x, row.y, row.z];
    end
else
    labels = {ci.label};
    hit = find_clinical_label_row_struct(labels, label);
    if isempty(hit)
        return;
    end
    if isfield(ci, 'mni_linear_x')
        mni = [ci(hit).mni_linear_x, ci(hit).mni_linear_y, ci(hit).mni_linear_z];
    end
end

ok = numel(mni) == 3 && all(isfinite(mni));
end


function hit = find_clinical_label_row(ci, label)
hit = find(strcmpi(strtrim(ci.label), strtrim(label)), 1, 'first');
if ~isempty(hit)
    return;
end
norm = @(s) regexprep(lower(strtrim(string(s))), '[^a-z0-9]', '');
target = norm(label);
for i = 1:height(ci)
    if norm(ci.label(i)) == target
        hit = i;
        return;
    end
end
hit = [];
end


function hit = find_clinical_label_row_struct(labels, label)
hit = find(strcmpi(strtrim(labels), strtrim(label)), 1, 'first');
if ~isempty(hit)
    return;
end
norm = @(s) regexprep(lower(strtrim(s)), '[^a-z0-9]', '');
target = norm(label);
for i = 1:numel(labels)
    if norm(labels{i}) == target
        hit = i;
        return;
    end
end
hit = [];
end


function arr = get_object_array(obj, name)
arr = [];
if isprop(obj, name)
    val = obj.(name);
    if isnumeric(val) && ~isempty(val)
        arr = val;
    end
end
end


function labels = get_object_cellstr(obj, name)
labels = {};
if isprop(obj, name)
    val = obj.(name);
    if iscell(val) && ~isempty(val)
        labels = val;
    end
end
end


function coords = anatomy_electrode_coords(sn_obj, space)
space = lower(space);
switch space
    case 'mni'
        S = sn_obj.anatomy.mni_space;
    case 'subject'
        S = sn_obj.anatomy.subject_space;
    otherwise
        error('space must be ''mni'' or ''subject''.');
end
if isfield(S, 'tala') && isfield(S.tala, 'electrodes')
    coords = S.tala.electrodes;
elseif isfield(S, 'vera_mat_minimal') && isfield(S.vera_mat_minimal, 'tala')
    coords = S.vera_mat_minimal.tala.electrodes;
else
    error('Could not find electrode coordinates in anatomy.%s_space', space);
end
end
