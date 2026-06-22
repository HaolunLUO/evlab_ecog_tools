function [sn_obj, meta] = load_sn_obj_for_comparison(subject, crunchedDir, repoRoot, taskType, auxSearchDirs, cohortKind)
% LOAD_SN_OBJ_FOR_COMPARISON  Load sn_obj from a finished analysis .mat file.
%
%   cohortKind:
%     'brainstorm'  -> <subject>_<taskType>_crunched.mat
%                      or legacy <subject>_MITLangloc_crunched.mat
%     'broadband'   -> <subject>_<taskType>_crunched_HG_ZScore.mat (other lab)

if nargin < 3 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 4 || isempty(taskType)
    taskType = 'MITLangloc';
end
if nargin < 5 || isempty(auxSearchDirs)
    auxSearchDirs = {crunchedDir, fileparts(crunchedDir), 'F:\iEEG_evlab'};
end
if nargin < 6 || isempty(cohortKind)
    cohortKind = 'brainstorm';
end

setup_comparison_path(repoRoot);

candidates = analysis_file_candidates(subject, crunchedDir, taskType, cohortKind);

matFile = '';
for i = 1:numel(candidates)
    if isfile(candidates{i}) && has_sn_obj_file(candidates{i})
        matFile = candidates{i};
        break;
    end
end

if isempty(matFile)
    error(['No analysis file with sn_obj found for %s in\n  %s\n' ...
        'Tried:\n  %s'], subject, crunchedDir, strjoin(candidates, '\n  '));
end

S = load(matFile, 'sn_obj', 'obj');
sn_obj = S.sn_obj;

if isfield(S, 'obj') && ~isempty(S.obj)
    sn_obj = enrich_sn_obj_from_obj(sn_obj, S.obj);
end

sn_obj = attach_auxiliary_mni(sn_obj, subject, taskType, auxSearchDirs);

meta = struct();
meta.file = matFile;
meta.subject = subject;
meta.taskType = taskType;
meta.cohortKind = cohortKind;
if isprop(sn_obj, 'subject') && ~isempty(sn_obj.subject)
    meta.subject = sn_obj.subject;
end
if isprop(sn_obj, 'experiment')
    meta.experiment = sn_obj.experiment;
end

end


function candidates = analysis_file_candidates(subject, crunchedDir, taskType, cohortKind)
switch lower(cohortKind)
    case 'brainstorm'
        candidates = {
            fullfile(crunchedDir, sprintf('%s_%s_crunched.mat', subject, taskType))
            fullfile(crunchedDir, sprintf('%s_MITLangloc_crunched.mat', subject))
        };
    case 'broadband'
        candidates = {
            fullfile(crunchedDir, sprintf('%s_%s_crunched_HG_ZScore.mat', subject, taskType))
            fullfile(crunchedDir, sprintf('%s_%s_crunched_HG.mat', subject, taskType))
            fullfile(crunchedDir, sprintf('%s_MITLangloc_crunched_HG_ZScore.mat', subject))
            fullfile(crunchedDir, sprintf('%s_MITLangloc_crunched_HG.mat', subject))
        };
    otherwise
        error('Unknown cohortKind: %s', cohortKind);
end
end


function ok = has_sn_obj_file(matFile)
ok = false;
try
    vars = whos('-file', matFile);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function setup_comparison_path(repoRoot)
persistent pathReady lastRepo
if ~isempty(pathReady) && strcmp(lastRepo, repoRoot)
    return;
end
scriptDir = fullfile(repoRoot, 'brainstorm_pipeline');
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));
pathReady = true;
lastRepo = repoRoot;
end
