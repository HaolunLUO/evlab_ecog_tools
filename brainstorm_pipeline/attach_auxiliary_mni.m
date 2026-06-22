function sn_obj = attach_auxiliary_mni(sn_obj, subject, taskType, searchDirs)
% ATTACH_AUXILIARY_MNI  Pull MNI coords from broadband / source crunched files.

if nargin < 4 || isempty(searchDirs)
    searchDirs = {};
end

if sn_obj_has_direct_mni(sn_obj)
    return;
end

files = auxiliary_mni_candidate_files(subject, taskType, searchDirs);
for i = 1:numel(files)
    if ~isfile(files{i})
        continue;
    end
    try
        S = load(files{i}, 'obj');
        if isfield(S, 'obj') && ~isempty(S.obj)
            sn_obj = enrich_sn_obj_from_obj(sn_obj, S.obj);
        end
    catch
    end
    if sn_obj_has_direct_mni(sn_obj)
        return;
    end
end
end


function tf = sn_obj_has_direct_mni(sn_obj)
tf = has_usable_prop(sn_obj, 'elec_ch_pos_mni') ...
    || has_usable_prop(sn_obj, 'bip_ch_pos_mni') ...
    || has_usable_clinical_info(sn_obj) ...
    || has_usable_anatomy_mni(sn_obj);
end


function tf = has_usable_anatomy_mni(sn_obj)
tf = false;
if ~isprop(sn_obj, 'anatomy') || isempty(sn_obj.anatomy) || ~isstruct(sn_obj.anatomy)
    return;
end
if ~isfield(sn_obj.anatomy, 'mapping') || isempty(sn_obj.anatomy.mapping)
    return;
end
if ~isfield(sn_obj.anatomy, 'mni_space') || isempty(sn_obj.anatomy.mni_space)
    return;
end
ms = sn_obj.anatomy.mni_space;
if isfield(ms, 'tala') && isfield(ms.tala, 'electrodes') && ~isempty(ms.tala.electrodes)
    tf = true;
elseif isfield(ms, 'vera_mat_minimal') && isfield(ms.vera_mat_minimal, 'tala') ...
        && isfield(ms.vera_mat_minimal.tala, 'electrodes') ...
        && ~isempty(ms.vera_mat_minimal.tala.electrodes)
    tf = true;
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


function tf = has_usable_clinical_info(sn_obj)
tf = false;
if ~isprop(sn_obj, 'filt_ops') || isempty(sn_obj.filt_ops)
    return;
end
ops = sn_obj.filt_ops;
tf = isstruct(ops) && isfield(ops, 'elec_clinic_info') && ~isempty(ops.elec_clinic_info);
end


function files = auxiliary_mni_candidate_files(subject, taskType, searchDirs)
files = {};
patterns = {
    '%s_%s_crunched_defaultSEEGorBOTHBroadBand.mat'
    '%s_%s_crunched_*BroadBand*.mat'
    '%s_*%s*BroadBand*.mat'
    '%s_%s_crunched.mat'
    };

for d = 1:numel(searchDirs)
    if isempty(searchDirs{d}) || ~isfolder(searchDirs{d})
        continue;
    end
    for p = 1:numel(patterns)
        hits = dir(fullfile(searchDirs{d}, sprintf(patterns{p}, subject, taskType)));
        for h = 1:numel(hits)
            files{end+1} = fullfile(hits(h).folder, hits(h).name); %#ok<AGROW>
        end
    end
end

[~, ia] = unique(files, 'stable');
files = files(ia);
end
