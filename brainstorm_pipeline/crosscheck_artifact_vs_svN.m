%% ========================================================================
% CROSS-CHECK SHARP ARTIFACTS vs S-vs-N SIGNIFICANT CHANNELS
% ========================================================================
% Loads a crunched .mat file produced by complete_mit_pipeline_brainstorm.m
% and reports overlap between:
%   - sn_obj.stats.artifact_stats_*  (sharp HG transients)
%   - sn_obj.s_vs_n_sig              (language localization significance)
%
% PREREQUISITES
%   addpath(genpath('<repo>/brainstorm_pipeline'));
%   Pipeline must have reached Step 7.5+ (test_s_vs_n or ROI applied).
%   artifact_stats are created when detectSharpArtifacts = true (default).
%
% OUTPUT
%   results  - struct with .unipolar / .bipolar tables and overlap indices
%   Optional CSV export and overlap-channel plots (see USER SETTINGS).
%
% See also: artifact_svN_crosscheck, load_sn_obj_from_crunched,
%           plot_crosscheck_channels

%% USER SETTINGS — edit paths before running
subject    = 'Subject06';
taskType   = 'MITSWJNTask';
workingDir = 'F:\seeg\luohong\analysisEV\';

montage = 'bipolar';              % 'unipolar', 'bipolar', or 'both'
runDetectionIfMissing = true;  % call detect_sharp_artifacts if stats absent
min_amplitude = 15;            % z-score threshold for detection
min_slope     = 10;            % z-score/sample slope threshold

doPlotOverlap = true;          % plot HG traces for sig+artifact channels
maxOverlapPlots = 6;

doSaveCsv = true;
csvOutputDir = workingDir;     % where to write summary tables


%% BUILD FILE PATH
crunchedFile = fullfile(workingDir, sprintf('%s_%s_crunched.mat', subject, taskType));
fprintf('Loading: %s\n', crunchedFile);


%% LOAD + CROSS-CHECK
sn_obj = load_sn_obj_from_crunched(crunchedFile);

results = artifact_svN_crosscheck(sn_obj, ...
    'montage', montage, ...
    'runDetectionIfMissing', runDetectionIfMissing, ...
    'min_amplitude', min_amplitude, ...
    'min_slope', min_slope);


%% OPTIONAL: PLOT OVERLAPPING CHANNELS
if doPlotOverlap
    if isfield(results, 'unipolar') && results.unipolar.n_overlap > 0
        plot_crosscheck_channels(sn_obj, results.unipolar.overlap_idx, 'unipolar', maxOverlapPlots);
    end
    if isfield(results, 'bipolar') && isfield(results.bipolar, 'overlap_idx') ...
            && results.bipolar.n_overlap > 0
        plot_crosscheck_channels(sn_obj, results.bipolar.overlap_idx, 'bipolar', maxOverlapPlots);
    end
end


%% OPTIONAL: EXPORT TABLES
if doSaveCsv
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    if isfield(results, 'unipolar') && ~isempty(results.unipolar.table)
        uniFile = fullfile(csvOutputDir, ...
            sprintf('%s_%s_artifact_svN_unipolar_%s.csv', subject, taskType, stamp));
        writetable(results.unipolar.table, uniFile);
        fprintf('Saved: %s\n', uniFile);

        if results.unipolar.n_overlap > 0
            ovFile = fullfile(csvOutputDir, ...
                sprintf('%s_%s_artifact_svN_overlap_uni_%s.csv', subject, taskType, stamp));
            writetable(results.unipolar.overlap, ovFile);
            fprintf('Saved: %s\n', ovFile);
        end
    end

    if isfield(results, 'bipolar') && isfield(results.bipolar, 'table') ...
            && ~isempty(results.bipolar.table)
        bipFile = fullfile(csvOutputDir, ...
            sprintf('%s_%s_artifact_svN_bipolar_%s.csv', subject, taskType, stamp));
        writetable(results.bipolar.table, bipFile);
        fprintf('Saved: %s\n', bipFile);

        if results.bipolar.n_overlap > 0
            ovFile = fullfile(csvOutputDir, ...
                sprintf('%s_%s_artifact_svN_overlap_bip_%s.csv', subject, taskType, stamp));
            writetable(results.bipolar.overlap, ovFile);
            fprintf('Saved: %s\n', ovFile);
        end
    end
end

fprintf('Done. Inspect results.unipolar.overlap and results.bipolar.overlap in the workspace.\n');
