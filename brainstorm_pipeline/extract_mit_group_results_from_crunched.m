%% ========================================================================
% EXTRACT MIT GROUP RESULTS FROM CRUNCHED FILES
% extract_mit_group_results_from_crunched.m
% =========================================================================
% PURPOSE
%   Batch-write MIT_multi_single-compatible outputs from finished crunched
%   .mat files (after complete_mit_pipeline_brainstorm.m has already run):
%
%     output/<taskType>/<Subject>_<taskType>_groupResult.mat
%     output/<taskType>/<Subject>_<taskType>_allChanGammaPower.mat
%
%   Use this when localization/plots are already in the crunched file but
%   groupResult / allChanGammaPower were not saved (or need to be refreshed).
%
% INPUT
%   <workingDir>/<Subject>_<taskType>_crunched.mat  (must contain sn_obj
%   with s_vs_n_sig populated)
%
% USAGE
%   1) Edit workingDir / taskType / subjects below
%   2) Run this script once to process all subjects
%% ========================================================================

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(scriptDir);
addpath(genpath(fullfile(repoRoot, 'brainstorm_pipeline')));
addpath(genpath(repoRoot));

%% USER SETTINGS
workingDir  = 'F:\seeg\luohong\analysisEV';
taskType    = 'MITSWJNTask';
crunchedDir = workingDir;
outputDir   = fullfile(workingDir, 'output', taskType);

% Leave empty to auto-discover *_<taskType>_crunched.mat in crunchedDir.
subjects = {};

overwrite            = true;
useOddForInference   = false;
includeEffectSizes   = true;
includeWordwise      = true;
includeWordBoundaries = false;

%% SETUP
taskConfig = get_mit_task_config(taskType);

if isempty(subjects)
    subjects = discover_subjects_from_crunched(crunchedDir, taskType);
end
if isempty(subjects)
    error('No subjects found. Set subjects manually or check crunchedDir:\n  %s', crunchedDir);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

saveOpts = struct();
saveOpts.includeEffectSizes    = includeEffectSizes;
saveOpts.includeWordwise       = includeWordwise;
saveOpts.includeWordBoundaries = includeWordBoundaries;
saveOpts.useOddForInference    = useOddForInference;
saveOpts.overwrite             = overwrite;
saveOpts.saveGammaPower        = true;

fprintf('\n=== EXTRACT GROUP RESULTS FROM CRUNCHED ===\n');
fprintf('Task: %s\n', taskType);
fprintf('Crunched dir: %s\n', crunchedDir);
fprintf('Output dir: %s\n', outputDir);
fprintf('Subjects (%d): %s\n\n', numel(subjects), strjoin(subjects, ', '));

okCount = 0;
failCount = 0;

for i = 1:numel(subjects)
    subjectName = subjects{i};
    crunchedFile = fullfile(crunchedDir, sprintf('%s_%s_crunched.mat', subjectName, taskType));

    fprintf('--- %s ---\n', subjectName);
    if ~isfile(crunchedFile)
        warning('Crunched file not found:\n  %s', crunchedFile);
        failCount = failCount + 1;
        continue;
    end

    try
        sn_obj = load_sn_obj_from_crunched(crunchedFile);
        if ~istable(sn_obj.s_vs_n_sig)
            error(['s_vs_n_sig missing on loaded object. Re-run localization ' ...
                '(complete_mit_pipeline_brainstorm Step 8) for this subject.']);
        end

        saveOpts.sourceFile = crunchedFile;
        outFiles = save_mit_group_outputs(sn_obj, subjectName, taskType, taskConfig, outputDir, saveOpts);

        nSig = sum(logical(sn_obj.s_vs_n_sig.elec_data{1}));
        nSigBip = 0;
        if ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
            nSigBip = sum(logical(sn_obj.s_vs_n_sig.bip_elec_data{1}));
        end
        fprintf('  sig uni=%d  sig bip=%d\n', nSig, nSigBip);
        fprintf('  -> %s\n', outFiles.groupResultFile);
        if ~isempty(outFiles.gammaPowerFile)
            fprintf('  -> %s\n', outFiles.gammaPowerFile);
        end
        okCount = okCount + 1;
    catch ME
        warning('Failed for %s: %s', subjectName, ME.message);
        failCount = failCount + 1;
    end
end

fprintf('\n=== DONE ===\n');
fprintf('Succeeded: %d | Failed: %d\n', okCount, failCount);


function subjects = discover_subjects_from_crunched(crunchedDir, taskType)
d = dir(fullfile(crunchedDir, sprintf('*%s_crunched.mat', taskType)));
subjects = cell(numel(d), 1);
suffix = sprintf('_%s_crunched.mat', taskType);
for i = 1:numel(d)
    name = d(i).name;
    if endsWith(name, suffix)
        subjects{i} = erase(name, suffix);
    else
        subjects{i} = '';
    end
end
subjects = subjects(~cellfun(@isempty, subjects));
subjects = unique(subjects, 'stable');
end
