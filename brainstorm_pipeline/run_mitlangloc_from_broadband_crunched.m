%% ========================================================================
% MITLangloc FROM BROADBAND-CRUNCHED DATA
% run_mitlangloc_from_broadband_crunched.m
% =========================================================================
% PURPOSE
%   Continue analysis from an external lab crunched file that already has
%   broadband preprocessing (defaultSEEGorBOTHBroadBand), e.g.
%   AMC092_MITLangloc_crunched_defaultSEEGorBOTHBroadBand.mat
%
%   This script:
%     1) Loads the broadband crunched .mat (variable: obj)
%     2) Extracts high-gamma (NAPLAB), downsamples
%     3) Baseline z-scores the HG envelope
%     4) Runs MITLangloc S-vs-N localization, plots, and group results
%
%   Does NOT modify complete_mit_pipeline_brainstorm.m.
%
% INPUT
%   A merged ieeg_pipeline crunched file with broadband referencing already
%   applied (obj.for_preproc.order present, sample_freq still ~500 Hz).
%
% MATLAB PATH
%   addpath(genpath('<repo>'))
%% ========================================================================

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
% Do NOT addpath(genpath(repoRoot)) here — legacy ecog_data.m shadows
% ieeg_pipeline @ecog_data and breaks loading of broadband crunched files.
addpath(scriptDir);
addpath(genpath(fullfile(repoRoot, 'brainstorm_pipeline')));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));
addpath(genpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master', 'utils', 'kumar_ieeg_utils')));

%% ========================================================================
% USER SETTINGS  (EDIT THESE)
%% ========================================================================
taskType = 'MITLangloc';

% Full path to broadband-referenced crunched file from the other lab:
broadbandCrunchedFile = 'F:\iEEG_evlab\BJH011_MITLangloc_crunched_defaultSEEGorBOTHBroadBand.mat';  % e.g. 'F:\...\AMC092_MITLangloc_crunched_defaultSEEGorBOTHBroadBand.mat'

workingDir  = 'F:\seeg\luohong\analysisEV';
anatomyPath = fullfile(workingDir, 'anatomy\');

% Where to save the HG+z-scored working copy (original file is not overwritten):
saveProcessedCopy = true;
processedCopyDir    = fullfile(workingDir, 'crunched', 'MITLangloc');

% Re-run HG extraction even if sample_freq already equals decimationFreq:
forceReprocessHG = false;

decimationFreq = 200;
detectSharpArtifacts = true;

% S-vs-N contrast for MITLangloc:
%   'Jabberwocky'    -> MGH / ieeg_pipeline MITLangloc (Sentences vs Jabberwocky)
%   'Nonword-lists'  -> only if your file includes that condition
N_condition = 'Jabberwocky';

% --- Language channel selection / statistics ---
useOddForInference = false;
computeEffectSizes = true;
doWordwiseLangloc  = true;
wordwiseMinConsec  = 3;
doWordBoundaries   = false;
wordBoundaryEpoch  = [-0.25 0.25];

% --- Optional PDF report (needs MATLAB Report Generator) ---
generateLanglocReport = true;
useLanglocReportV2    = true;

% --- Optional behavioral log ---
behaviorFile    = '';
behaviorAlignBy = 'order';

%% ========================================================================
% TASK CONFIG (MITLangloc)
%% ========================================================================
taskConfig = struct();
taskConfig.nWordPositions = 12;
taskConfig.wordDuration   = 0.45;
taskConfig.eventPattern   = '%s_tp%d';
taskConfig.conditionMap = containers.Map(...
    {'sentence', 'word lists', 'Jabberwocky', 'non-word lists'}, ...
    {'Sentences', 'Word-lists', 'Jabberwocky', 'Nonword-lists'});
taskConfig.S_condition = 'Sentences';
taskConfig.N_condition = N_condition;
taskConfig.W_condition = 'Word-lists';
taskConfig.J_condition = 'Jabberwocky';
taskConfig.testWords   = 1:12;
taskConfig.subAverage  = true;

%% ========================================================================
% SETUP
%% ========================================================================
fprintf('\n=== SETUP ===\n');

if isempty(broadbandCrunchedFile) || ~isfile(broadbandCrunchedFile)
    error(['Set broadbandCrunchedFile to a valid .mat file, e.g.\n' ...
        '  AMC092_MITLangloc_crunched_defaultSEEGorBOTHBroadBand.mat']);
end

cd(workingDir);

outputDir = fullfile(workingDir, 'output', taskType);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end
if saveProcessedCopy && ~exist(processedCopyDir, 'dir'); mkdir(processedCopyDir); end

fprintf('Input:  %s\n', broadbandCrunchedFile);
fprintf('Output: %s\n', outputDir);

%% ========================================================================
% STEP 1: LOAD (ieeg_pipeline ecog_data, not legacy ecog_data.m)
%% ========================================================================
fprintf('\n=== STEP 1: LOAD BROADBAND CRUNCHED FILE ===\n');

obj = load_broadband_crunched(broadbandCrunchedFile, repoRoot);

if ~isa(obj, 'ecog_data_seeg')
    fprintf('Wrapping %s as ecog_data_seeg ...\n', class(obj));
    obj = wrap_crunched_obj_as_seeg(obj);
end

subjectName = obj.subject;
if isempty(subjectName)
    [~, baseName] = fileparts(broadbandCrunchedFile);
    subjectName = regexp(baseName, '^([A-Za-z0-9]+)', 'match', 'once');
end
if isempty(subjectName)
    error('Could not determine subject name from obj.subject or filename.');
end
fprintf('Subject: %s | Fs: %.1f Hz\n', subjectName, obj.sample_freq);

if ~isempty(behaviorFile)
    if ~isfile(behaviorFile)
        error('behaviorFile not found:\n  %s', behaviorFile);
    end
    obj = attach_behavior_to_obj(obj, behaviorFile, 'alignBy', behaviorAlignBy);
end

%% ========================================================================
% STEP 2: ANATOMY
%% ========================================================================
fprintf('\n=== STEP 2: ANATOMY ===\n');

if (~isfield(obj, 'anatomy') || isempty(obj.anatomy)) && isfolder(anatomyPath)
    try
        obj.add_anatomy(anatomyPath);
        fprintf('Anatomy loaded from %s\n', anatomyPath);
    catch ME
        warning('Could not load anatomy: %s', ME.message);
    end
else
    fprintf('Using anatomy already in obj (or anatomyPath not set).\n');
end

%% ========================================================================
% STEP 3: HIGH-GAMMA + DOWNSAMPLE  (skip broadband)
%% ========================================================================
fprintf('\n=== STEP 3: HIGH-GAMMA + DOWNSAMPLE ===\n');

needHGExtraction = forceReprocessHG || obj.sample_freq > decimationFreq;

if needHGExtraction
    fprintf('Extracting high-gamma (NAPLAB) and downsampling to %d Hz ...\n', decimationFreq);
    obj.extract_high_gamma('doNapLabFilterExtraction', true);
    obj.downsample_signal('decimationFreq', decimationFreq);
else
    fprintf('High-gamma + downsampling already done (Fs=%.1f Hz); skipping.\n', obj.sample_freq);
end

fprintf('Unipolar: [%d x %d] | Bipolar: [%d x %d] | Fs: %.1f Hz\n', ...
    size(obj.elec_data,1), size(obj.elec_data,2), ...
    size(obj.bip_elec_data,1), size(obj.bip_elec_data,2), obj.sample_freq);

%% ========================================================================
% STEP 4: CONDITIONS, EVENTS, TRIALS
%% ========================================================================
fprintf('\n=== STEP 4: CONDITIONS + EVENTS + TRIALS ===\n');

[obj, taskConfig, eventReport] = normalize_mitlangloc_events(obj, taskConfig);
fprintf('Using contrast: %s vs %s\n', taskConfig.S_condition, taskConfig.N_condition);

if isempty(obj.trial_data)
    obj.make_trials();
    fprintf('Created %d trial segments.\n', numel(obj.trial_data));
else
    fprintf('Trials already present: %d\n', numel(obj.trial_data));
end

if saveProcessedCopy
    crunchedFile = fullfile(processedCopyDir, ...
        sprintf('%s_%s_crunched_HG.mat', subjectName, taskType));
else
    crunchedFile = broadbandCrunchedFile;
end
save(crunchedFile, 'obj', '-v7.3');
fprintf('Saved pre-zscore object: %s\n', crunchedFile);

%% ========================================================================
% STEP 5: ANALYSIS OBJECT
%% ========================================================================
fprintf('\n=== STEP 5: ANALYSIS OBJECT ===\n');

sn_obj = ecog_sn_data_seeg(...
    outputDir, ...
    crunchedFile, ...
    workingDir, ...
    'ecog_data_seeg.m', ...
    workingDir);
sn_obj.anatomy = obj.anatomy;

if ~exist(sn_obj.langloc_save_path, 'dir')
    mkdir(sn_obj.langloc_save_path);
end

%% ========================================================================
% STEP 6: SIGNIFICANCE + BASELINE Z-SCORE
%% ========================================================================
fprintf('\n=== STEP 6: SIGNIFICANCE + BASELINE Z-SCORE ===\n');

if isempty(sn_obj.stats) || ~isstruct(sn_obj.stats)
    sn_obj.stats = struct();
end

sn_obj.extract_significant_channel();
sn_obj.extract_time_significance();
sn_obj.extract_normalization_metrics();
sn_obj.normalize_signal('normtype', 'z-score');
sn_obj.make_trials();

if detectSharpArtifacts
    sn_obj.detect_sharp_artifacts();
end

fprintf('Baseline z-score applied.\n');

%% ========================================================================
% STEP 7: S vs N LOCALIZATION
%% ========================================================================
fprintf('\n=== STEP 7: S vs N LOCALIZATION ===\n');
fprintf('Contrast: %s vs %s\n', taskConfig.S_condition, taskConfig.N_condition);

sn_obj.test_s_vs_n('words', taskConfig.testWords, ...
    'S_condition_flag', taskConfig.S_condition, ...
    'N_condition_flag', taskConfig.N_condition, ...
    'n_rep', 10000, ...
    'threshold', 0.05, ...
    'side', 'right', ...
    'use_odd_for_inference', useOddForInference, ...
    'do_plot', true);

if computeEffectSizes
    sn_obj.compute_hg_power_diff_s_vs_n('words', taskConfig.testWords, ...
        'S_condition_flag', taskConfig.S_condition, ...
        'N_condition_flag', taskConfig.N_condition, ...
        'use_odd_for_inference', useOddForInference);
    sn_obj.compute_hg_sn_corr('words', taskConfig.testWords, ...
        'S_condition_flag', taskConfig.S_condition, ...
        'N_condition_flag', taskConfig.N_condition);
end

if doWordwiseLangloc
    sn_obj.test_s_vs_n_wordwise('words', taskConfig.testWords, ...
        'S_condition_flag', taskConfig.S_condition, ...
        'N_condition_flag', taskConfig.N_condition, ...
        'n_rep', 1000, ...
        'corr_type', 'Spearman', ...
        'threshold', 0.05, ...
        'side', 'right', ...
        'min_consecutive', wordwiseMinConsec, ...
        'consecutiveness', true, ...
        'use_odd_for_inference', useOddForInference);
end

if doWordBoundaries
    sn_obj.test_s_vs_n_wordboundaries(...
        'S_condition_flag', taskConfig.S_condition, ...
        'N_condition_flag', taskConfig.N_condition, ...
        'n_rep', 1000, ...
        'threshold', 0.05, ...
        'epoch_range', wordBoundaryEpoch, ...
        'num_words', numel(taskConfig.testWords), ...
        'use_odd_for_inference', useOddForInference, ...
        'do_plot', false);
end

if ~istable(sn_obj.s_vs_n_sig)
    error('sn_obj.s_vs_n_sig is missing or invalid.');
end

%% ========================================================================
% STEP 8: PLOTS + SUMMARY
%% ========================================================================
fprintf('\n=== STEP 8: PLOTS ===\n');

plotArgs = {'words', taskConfig.testWords, ...
    'S_condition_flag', taskConfig.S_condition, ...
    'N_condition_flag', taskConfig.N_condition, ...
    'subAverage', taskConfig.subAverage};
if ~isempty(taskConfig.W_condition)
    plotArgs = [plotArgs, {'W_condition_flag', taskConfig.W_condition}];
end
if ~isempty(taskConfig.J_condition)
    plotArgs = [plotArgs, {'J_condition_flag', taskConfig.J_condition}];
end
sn_obj.lang_resp_plots(plotArgs{:});

fprintf('\n=== SUMMARY STATISTICS ===\n');
try
    disp(sn_obj.get_summary_statistics());
catch ME
    warning('Summary statistics failed: %s', ME.message);
end

roiUni = find(sn_obj.s_vs_n_sig.elec_data{1});
fprintf('Significant unipolar channels: %d (%.1f%% of clean)\n', ...
    numel(roiUni), 100*numel(roiUni)/numel(sn_obj.elec_ch_clean));

%% ========================================================================
% STEP 9: SAVE RESULTS
%% ========================================================================
fprintf('\n=== STEP 9: SAVE RESULTS ===\n');

groupResultFile = fullfile(outputDir, sprintf('%s_%s_groupResult.mat', subjectName, taskType));
groupResult = struct();
groupResult.subject  = subjectName;
groupResult.taskType = taskType;
groupResult.sourceFile = broadbandCrunchedFile;
groupResult.sig_uni     = logical(sn_obj.s_vs_n_sig.elec_data{1});
groupResult.p_ratio_uni = sn_obj.s_vs_n_p_ratio.elec_data{1};
groupResult.elec_labels = sn_obj.elec_ch_label;
groupResult.nClean      = numel(sn_obj.elec_ch_clean);
groupResult.nSig        = sum(groupResult.sig_uni);
groupResult.sample_freq = sn_obj.sample_freq;
groupResult.nWords      = numel(taskConfig.testWords);
save(groupResultFile, 'groupResult', '-v7.3');
fprintf('Saved group result: %s\n', groupResultFile);

obj = sync_obj_from_sn_obj(obj, sn_obj);
if saveProcessedCopy
    hgZscoreFile = fullfile(processedCopyDir, ...
        sprintf('%s_%s_crunched_HG_ZScore.mat', subjectName, taskType));
    save(hgZscoreFile, 'obj', 'sn_obj', '-v7.3');
    fprintf('Saved HG+z-scored object: %s\n', hgZscoreFile);
else
    save(crunchedFile, 'obj', 'sn_obj', '-v7.3');
end

%% ========================================================================
% STEP 10: OPTIONAL PDF REPORT
%% ========================================================================
if generateLanglocReport
    fprintf('\n=== STEP 10: LANGLOC PDF REPORT ===\n');
    try
        run_langloc_report(obj, taskType, taskConfig, outputDir, ...
            subjectName, useLanglocReportV2);
    catch ME
        warning('Langloc report generation failed: %s', ME.message);
    end
end

fprintf('\n========================================\n');
fprintf('MITLangloc analysis complete\n');
fprintf('Subject: %s\n', subjectName);
fprintf('Results: %s\n', sn_obj.langloc_save_path);
fprintf('========================================\n');
