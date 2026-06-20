function obj = load_broadband_crunched(matFile, repoRoot)
% LOAD_BROADBAND_CRUNCHED  Load a merged-pipeline crunched .mat as ieeg ecog_data.
%
%   obj = load_broadband_crunched(matFile)
%   obj = load_broadband_crunched(matFile, repoRoot)
%
%   Crunched files from the ieeg_pipeline are saved as class ecog_data
%   (@ecog_data). The legacy ecog_data.m in the repo root has the SAME class
%   name and will shadow @ecog_data when the repo root is on the MATLAB path
%   (including when the current folder is the repo). Loading with the wrong
%   class drops elec_data / condition / for_preproc and breaks downstream steps.
%
%   This helper loads from a neutral working directory with only the
%   brainstorm + ieeg_pipeline paths, so @ecog_data is used.

if nargin < 2 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end

if ~isfile(matFile)
    error('Crunched file not found:\n  %s', matFile);
end

oldPwd = pwd;
cleanupPwd = onCleanup(@() cd(oldPwd)); %#ok<NASGU>
cd(tempdir);

oldPath = path;
cleanupPath = onCleanup(@() path(oldPath)); %#ok<NASGU>

restoredefaultpath;
brainstormDir = fullfile(repoRoot, 'brainstorm_pipeline');
ieegRoot      = fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master');
utilsRoot     = fullfile(ieegRoot, 'utils', 'kumar_ieeg_utils');

addpath(genpath(brainstormDir));
addpath(ieegRoot);
if isfolder(utilsRoot)
    addpath(genpath(utilsRoot));
end

ieegClass = which('ecog_data');
if contains(ieegClass, [filesep 'evlab_ecog_tools' filesep 'ecog_data.m']) ...
        && ~contains(ieegClass, '@ecog_data')
    error(['Legacy ecog_data.m is still shadowing ieeg_pipeline @ecog_data.\n' ...
        '  Resolved: %s\n  Remove the repo root from the path before loading.'], ieegClass);
end

S = load(matFile);
if ~isfield(S, 'obj')
    error('Expected variable ''obj'' in:\n  %s', matFile);
end
obj = S.obj;

if ~isa(obj, 'ecog_data') && ~isa(obj, 'ecog_data_ieeg') && ~isa(obj, 'ecog_data_seeg')
    warning('Loaded object class is %s (expected ieeg_pipeline ecog_data).', class(obj));
end

if ~isprop(obj, 'elec_data') || isempty(obj.elec_data)
    error(['No broadband elec_data found after load. The file may have been read ' ...
        'with the wrong ecog_data class (legacy vs ieeg_pipeline @ecog_data).']);
end

end
