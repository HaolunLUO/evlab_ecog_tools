function labels = get_channel_hcp_labels(sn_obj, chanIdx, signalType)
% GET_CHANNEL_HCP_LABELS  HCP-MMP labels for one channel (0..N labels).
%
%   Unipolar: one label from elec_ch_HCP_label or filt_ops.elec_clinic_info.
%   Bipolar:  bip_ch_HCP_label when valid, otherwise labels of both contacts.

labels = {};

if nargin < 3 || isempty(signalType)
    signalType = 'bipolar';
end
if isempty(chanIdx) || chanIdx < 1
    return;
end

if strcmpi(signalType, 'bipolar')
    labels = bipolar_hcp_labels(sn_obj, chanIdx);
else
    labels = unipolar_hcp_labels(sn_obj, chanIdx);
end

labels = labels(~cellfun(@(s) isempty(s) || strcmpi(strtrim(s), 'N/A'), labels));
end


function labels = unipolar_hcp_labels(sn_obj, uniIdx)
labels = {};

if isprop(sn_obj, 'elec_ch_HCP_label') && ~isempty(sn_obj.elec_ch_HCP_label) ...
        && uniIdx <= numel(sn_obj.elec_ch_HCP_label)
    val = sn_obj.elec_ch_HCP_label{uniIdx};
    if ~isempty(val) && ~strcmpi(strtrim(char(string(val))), 'N/A')
        labels = {char(string(val))};
        return;
    end
end

elecLabels = get_object_cellstr(sn_obj, 'elec_ch_label');
if uniIdx <= numel(elecLabels)
    label = hcp_from_clinical(sn_obj, elecLabels{uniIdx});
    if ~isempty(label)
        labels = {label};
    end
end
end


function labels = bipolar_hcp_labels(sn_obj, bipIdx)
labels = {};

if isprop(sn_obj, 'bip_ch_HCP_label') && ~isempty(sn_obj.bip_ch_HCP_label) ...
        && bipIdx <= numel(sn_obj.bip_ch_HCP_label)
    val = sn_obj.bip_ch_HCP_label{bipIdx};
    if ~isempty(val) && ~strcmpi(strtrim(char(string(val))), 'N/A')
        labels = {char(string(val))};
        return;
    end
end

bipLabels = get_object_cellstr(sn_obj, 'bip_ch_label');
pairLabels = get_bipolar_pair_labels(sn_obj, bipIdx, bipLabels);
for i = 1:numel(pairLabels)
    label = hcp_from_clinical(sn_obj, pairLabels{i});
    if ~isempty(label)
        labels{end+1} = label; %#ok<AGROW>
    end
end
labels = unique(labels, 'stable');
end


function label = hcp_from_clinical(sn_obj, contactLabel)
label = '';
if ~isprop(sn_obj, 'filt_ops') || isempty(sn_obj.filt_ops)
    return;
end
ops = sn_obj.filt_ops;
if ~isstruct(ops) || ~isfield(ops, 'elec_clinic_info') || isempty(ops.elec_clinic_info)
    return;
end

ci = ops.elec_clinic_info;
if istable(ci)
    if ~ismember('label', ci.Properties.VariableNames)
        return;
    end
    hit = find(strcmpi(strtrim(ci.label), strtrim(contactLabel)), 1, 'first');
    if isempty(hit)
        hit = find_clinical_label_row(ci, contactLabel);
    end
    if isempty(hit)
        return;
    end
    row = ci(hit, :);
    if ismember('HCPMMP1_label_1', ci.Properties.VariableNames)
        label = char(string(row.HCPMMP1_label_1));
    end
elseif isstruct(ci) && isfield(ci, 'label')
    labels = {ci.label};
    hit = find(strcmpi(strtrim(labels), strtrim(contactLabel)), 1, 'first');
    if isempty(hit)
        return;
    end
    if isfield(ci, 'HCPMMP1_label_1')
        label = char(string(ci(hit).HCPMMP1_label_1));
    end
end

if strcmpi(strtrim(label), 'N/A')
    label = '';
end
end


function hit = find_clinical_label_row(ci, label)
norm = @(s) regexprep(lower(strtrim(string(s))), '[^a-z0-9]', '');
target = norm(label);
hit = [];
for i = 1:height(ci)
    if norm(ci.label(i)) == target
        hit = i;
        return;
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
    end
end
end
