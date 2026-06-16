function sn_obj = load_sn_obj_from_crunched(crunchedFile)
% LOAD_SN_OBJ_FROM_CRUNCHED  Load sn_obj from a pipeline crunched .mat file.
%
%   sn_obj = load_sn_obj_from_crunched(crunchedFile)
%
%   Prefers sn_obj; falls back to obj if sn_obj was not saved.

if ~isfile(crunchedFile)
    error('Crunched file not found:\n  %s', crunchedFile);
end

vars = whos('-file', crunchedFile);
varNames = {vars.name};
if ismember('sn_obj', varNames)
    load(crunchedFile, 'sn_obj');
elseif ismember('obj', varNames)
    warning('sn_obj not found in file — using obj instead.');
    load(crunchedFile, 'obj');
    sn_obj = obj;
else
    error('Expected sn_obj or obj in:\n  %s', crunchedFile);
end

if ~isa(sn_obj, 'ecog_sn_data_seeg') && ~isa(sn_obj, 'ecog_data_seeg')
    warning('Loaded object is %s, not ecog_sn_data_seeg.', class(sn_obj));
end

if ~istable(sn_obj.s_vs_n_sig)
    warning(['s_vs_n_sig is missing on the loaded object. ' ...
        'Run test_s_vs_n or apply_roi before cross-checking.']);
end
end
