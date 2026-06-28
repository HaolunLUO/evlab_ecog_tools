function outFiles = save_mit_group_outputs(sn_obj, subjectName, taskType, taskConfig, outputDir, opts)
% SAVE_MIT_GROUP_OUTPUTS  Write groupResult + allChanGammaPower .mat files.
%
%   outFiles = save_mit_group_outputs(sn_obj, 'Subject01', 'MITSWJNTask', ...
%       taskConfig, outputDir, struct('sourceFile', crunchedFile));
%
%   Files written (MIT_multi_single-compatible naming):
%     <outputDir>/<Subject>_<Task>_groupResult.mat
%     <outputDir>/<Subject>_<Task>_allChanGammaPower.mat

if nargin < 6 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'saveGammaPower')
    opts.saveGammaPower = true;
end
if ~isfield(opts, 'overwrite')
    opts.overwrite = true;
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

groupResultFile = fullfile(outputDir, sprintf('%s_%s_groupResult.mat', subjectName, taskType));
gammaPowerFile  = fullfile(outputDir, sprintf('%s_%s_allChanGammaPower.mat', subjectName, taskType));

if opts.overwrite || ~isfile(groupResultFile)
    groupResult = build_group_result_from_sn_obj(sn_obj, subjectName, taskType, taskConfig, opts);
    save(groupResultFile, 'groupResult', '-v7.3');
    fprintf('Saved group result: %s\n', groupResultFile);
else
    fprintf('Skipping existing group result: %s\n', groupResultFile);
end

if opts.saveGammaPower && ~isempty(sn_obj.bip_elec_data)
    if opts.overwrite || ~isfile(gammaPowerFile)
        try
            gammaPower = build_gamma_power_from_sn_obj(sn_obj, taskConfig, opts);
            save(gammaPowerFile, 'gammaPower', '-v7.3');
            fprintf('Saved gamma power: %s\n', gammaPowerFile);
        catch ME
            warning('Gamma power export failed for %s: %s', subjectName, ME.message);
            gammaPowerFile = '';
        end
    else
        fprintf('Skipping existing gamma power: %s\n', gammaPowerFile);
    end
else
    gammaPowerFile = '';
end

outFiles = struct();
outFiles.groupResultFile = groupResultFile;
outFiles.gammaPowerFile  = gammaPowerFile;
end
