%% ========================================================================
% BRAINSTORM -> MIT PIPELINE  (complete_mit_pipeline_brainstorm.m)
% =========================================================================
% PURPOSE
%   End-to-end iEEG/SEEG analysis from Brainstorm-exported data through the
%   MIT high-gamma pipeline, including cross-task ROI support.
%
% RELATIONSHIP TO EXISTING evlab_ecog_tools PIPELINE
%   The existing evlab_ecog_tools pipeline (crunch_subject_ALBANY.m) starts
%   from BCI2000 .dat files and produces crunched .mat files.  This script
%   starts from Brainstorm-exported .mat files (converted by
%   brainstorm_to_mit_crunched_new.m) and produces the same crunched .mat
%   format, making the two pipelines interoperable at the crunched-file
%   boundary.  See INTEGRATION_GUIDE.md for details.
%
% WORKFLOW
%   1) Convert Brainstorm SEEG data -> MIT crunched .mat file
%   2) (optional) Attach behavioral log -> obj.events_table
%   3) Preprocess: highpass, notch, CAR/bipolar, high-gamma extraction
%   4) (roiSourceTask) Run S vs N localization and save ROI
%      (other tasks)   Load ROI and apply to current task
%   5) Generate timecourse + barplots + anatomy plots
%   6) Save group-compatible result .mat files (groupResult + allChanGammaPower); sync obj; optional langloc PDF
%
% RECOMMENDED USAGE
%   A) Run once with taskType='MITSWJNTask' -> defines and saves ROI
%   B) Run with taskType='WM'        -> plots MITSWJNTask ROI on WM data
%   C) Run with taskType='Auditory'  -> plots MITSWJNTask ROI on Auditory
%
% MATLAB PATH REQUIREMENTS
%   addpath(genpath('<repo>'))                      % evlab_ecog_tools (incl. brainstorm_pipeline)
%   brainstorm_to_mit_crunched_new.m must be on path (provided separately)
%   JancaCodePapers/ must be on path (in this repo) for IED removal
%
%   NOTE: The SEEG classes are uniquely named (ecog_data_seeg /
%   ecog_sn_data_seeg). ecog_data_seeg now subclasses the advanced ieeg_pipeline
%   engine, vendored here as brainstorm_pipeline/@ecog_data_ieeg (uniquely named
%   to avoid the collision with the legacy ecog_data classes), so path ORDER
%   still does not matter.
%
%   Because of this re-base the inherited engine behaviour is in effect: e.g.
%   normalize_signal now Gaussian-smooths the high-gamma envelope. Optional
%   advanced steps are exposed as the opt-in flags in the USER SETTINGS below
%   (along-shank Laplacian referencing, sharp-artifact detection). See
%   INTEGRATION_GUIDE.md ("Re-based onto the ieeg_pipeline engine") for details.
%% ========================================================================

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(scriptDir);
addpath(genpath(fullfile(repoRoot, 'brainstorm_pipeline')));
addpath(genpath(repoRoot));
addpath(genpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master', 'utils', 'kumar_ieeg_utils')));
%% ========================================================================
% USER SETTINGS
%% ========================================================================
taskType = 'MITSWJNTask';  % 'MITSWJNTask' | 'WM' | 'vWM' | 'Math' | 'MSIT' | 'vMSIT' | 'Auditory' | 'MITLangloc'

% --- Cross-task ROI ---
useROIfromSource = false;
roiSourceTask    = 'MITSWJNTask';

% --- Re-run control ---
forceRebuildCrunched = false;
forceReprocess       = false;

% --- UI control ---
isPlotVisible        = false;
doneVisualInspection = true;

% --- Preprocessing / feature extraction (ieeg_pipeline engine) ---
% Broadband referencing -> NAPLAB high-gamma envelope -> downsample, followed
% (in STEP 7.5) by significance + baseline z-score. This chain is the input to
% the language channel selection.
% NOTE: the inherited engine normalize_signal also Gaussian-smooths the HG
% envelope.
preprocOrder   = 'defaultSEEGorBOTHBroadBand';  % highpass,notch,IED,CAR,Laplacian,bipolar (no envelope)
decimationFreq = 200;                            % Hz, for downsample_signal
detectSharpArtifacts = true;                    % optional sharp-transient QC after z-scoring

% --- Language channel selection / statistics (ecog_MITLangloc-master features) ---
% Split-half cross-validation: odd trials drive selection (inference), even
% trials are reserved for plotting and held-out effect sizes (avoids
% double-dipping). Applies to test_s_vs_n, lang_resp_plots, and the effect
% size / word-level analyses below.
useOddForInference   = false;
% Held-out effect sizes (computed on even trials):
computeEffectSizes   = true;   % compute_hg_power_diff_s_vs_n + compute_hg_sn_corr
% Word-level selection (per-word permutation + consecutiveness criterion):
doWordwiseLangloc    = true;   % test_s_vs_n_wordwise
wordwiseMinConsec    = 3;      % require this many (consecutive) significant words
% Word-boundaries cluster-based time-series permutation test:
doWordBoundaries     = false;  % test_s_vs_n_wordboundaries (slower; needs timePermCluster)
wordBoundaryEpoch    = [-0.25 0.25];

% --- Langloc PDF report (MITSWJNTask only) ---
% Requires MATLAB Report Generator (mlreportgen). Saves PDF to output/MITSWJNTask/.
generateLanglocReport = true;
useLanglocReportV2    = true;  % v2: includes LangLoc Responsive Electrodes chapter

% --- Behavioral log (MITSWJNTask only; for PDF performance metrics) ---
% CSV/table with accuracy + RT (seconds), or response + probe columns.
% Leave empty ('') to skip. Row order must match neural trials (default),
% or set behaviorAlignBy = 'session_trial' with session/trial columns.
workingDir  = 'F:\seeg\luohong\analysisEV';
behaviorFile    = fullfile(workingDir, 'behavior', 'Subject12_MITSWJNTask_behavior.csv');   % e.g. fullfile(workingDir, 'behavior', 'Subject01_MITSWJNTask.csv')
behaviorAlignBy = 'order';  % 'order' | 'session_trial'

% --- Paths (EDIT THESE) ---
workingDir  = 'F:\seeg\luohong\analysisEV';
anatomyPath = fullfile(workingDir, 'anatomy\');

% --- Subject/protocol ---
params = struct();
params.SubjectName  = 'Subject12';
params.ProtocolName = 'analysis';
params.outputPath   = workingDir;
params.taskType     = taskType;

% --- Data files (EDIT PER TASK) ---
switch taskType
    case 'MITSWJNTask'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011XQ/data_block001.mat'
        };
    case 'WM'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011X8/data_block001.mat'
        };
    case 'vWM'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011XQ/data_block001_05.mat'  % EDIT: vWM task data
        };
    case 'Math'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011XQ/data_block001_02.mat'  % EDIT: Math task data
        };
    case 'MSIT'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011XQ/data_block001_04.mat'  % EDIT: MSIT task data
        };
    case 'vMSIT'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject12/DA0011XQ/data_block001_03.mat'  % EDIT: vMSIT task data
        };
    case 'Auditory'
        allDataFiles = {
            'D:\seeg\analysis\data\Subject06\DA010035\data_block001_12.mat'
            'D:\seeg\analysis\data\Subject06\DA010035\data_block001_13.mat'
        };
    otherwise
        allDataFiles = {
            'D:\seeg\analysis\data\Subject03\DA01001N\data_block001_08.mat'
            'D:\seeg\analysis\data\Subject03\DA01001N\data_block001_07.mat'
        };
end


%% ========================================================================
% STEP 0: TASK CONFIGURATION
%% ========================================================================
taskConfig = get_mit_task_config(taskType);

params.nWordPositions = taskConfig.nWordPositions;
params.wordDuration   = taskConfig.wordDuration;
params.taskConfig     = taskConfig;


%% ========================================================================
% STEP 1: SETUP
%% ========================================================================
fprintf('\n=== STEP 1: SETUP ===\n');

cd(workingDir);
addpath(genpath(workingDir));

applyROI = useROIfromSource && ~strcmpi(taskType, roiSourceTask);

outputDir = fullfile(workingDir, 'output', taskType);
if applyROI
    outputDir = fullfile(workingDir, 'output', taskType, sprintf('ROIfrom_%s', roiSourceTask));
end
if ~exist(outputDir,'dir'); mkdir(outputDir); end

roiFile = fullfile(workingDir, 'output', roiSourceTask, ...
    sprintf('%s_ROI_from_%s.mat', params.SubjectName, roiSourceTask));

fprintf('Subject: %s\n', params.SubjectName);
fprintf('Task: %s\n', taskType);
fprintf('Apply ROI from %s? %d\n', roiSourceTask, applyROI);
fprintf('Data files: %d\n', numel(allDataFiles));
fprintf('OutputDir: %s\n', outputDir);


%% ========================================================================
% STEP 2: CONVERT BRAINSTORM -> MIT CRUNCHED FILE
%% ========================================================================
fprintf('\n=== STEP 2: CONVERTING BRAINSTORM -> MIT FORMAT ===\n');

defaultCrunchedFile = fullfile(workingDir, sprintf('%s_MITLangloc_crunched.mat', params.SubjectName));
taskCrunchedFile    = fullfile(workingDir, sprintf('%s_%s_crunched.mat', params.SubjectName, taskType));

if forceRebuildCrunched && exist(taskCrunchedFile,'file')
    fprintf('Deleting existing crunched file (forceRebuildCrunched=true)\n');
    delete(taskCrunchedFile);
end

if exist(taskCrunchedFile,'file')
    % Check that the saved object is the SEEG pipeline class (ecog_data_seeg).
    % A plain engine/ecog_data_v2 object uses the 60 Hz notch defaults; if we
    % find one here it means the file was created before the SEEG pipeline
    % existed and must be rebuilt with ecog_data_seeg (50 Hz notch, no CAR
    % before bipolar).
    tmp = load(taskCrunchedFile, 'obj');
    if ~isa(tmp.obj, 'ecog_data_seeg')
        fprintf(['WARNING: crunched file contains a ''%s'' object instead of ' ...
            '''ecog_data_seeg''.\n  Deleting and rebuilding with the SEEG ' ...
            'pipeline (50 Hz notch).\n'], class(tmp.obj));
        delete(taskCrunchedFile);
    else
        fprintf('Found existing crunched file. Skipping conversion.\n');
    end
    clear tmp;
end

if ~exist(taskCrunchedFile,'file')
    fprintf('Running brainstorm_to_mit_crunched_new...\n');
    brainstorm_to_mit_crunched_new(allDataFiles, params);

    if exist(defaultCrunchedFile,'file') && ~exist(taskCrunchedFile,'file')
        movefile(defaultCrunchedFile, taskCrunchedFile, 'f');
    end

    if ~exist(taskCrunchedFile,'file')
        error('Conversion did not produce expected crunched file:\n  %s', taskCrunchedFile);
    end
end

crunchedFile = taskCrunchedFile;
fprintf('Using crunched file: %s\n', crunchedFile);


%% ========================================================================
% STEP 2.5: ATTACH BEHAVIORAL LOG (MITSWJNTask only)
%% ========================================================================
if strcmpi(taskType, 'MITSWJNTask')
    fprintf('\n=== STEP 2.5: BEHAVIORAL DATA ===\n');

    load(crunchedFile, 'obj');

    if ~isempty(behaviorFile)
        if ~isfile(behaviorFile)
            error('behaviorFile not found:\n  %s', behaviorFile);
        end
        fprintf('Loading behavior: %s\n', behaviorFile);
        obj = attach_behavior_to_obj(obj, behaviorFile, 'alignBy', behaviorAlignBy);
        save(crunchedFile, 'obj', '-v7.3');
    else
        if isprop(obj, 'events_table') && istable(obj.events_table) ...
                && height(obj.events_table) == numel(obj.condition)
            fprintf('No behaviorFile set; using existing obj.events_table (%d trials).\n', ...
                height(obj.events_table));
        else
            fprintf('No behaviorFile set; obj.events_table not populated (report will use defaults).\n');
        end
    end
else
    fprintf('\n=== STEP 2.5: BEHAVIORAL DATA (skipped; MITSWJNTask only) ===\n');
end


%% ========================================================================
% STEP 3: FIX CHANNEL TYPES + CLEAN CHANNELS
%% ========================================================================
fprintf('\n=== STEP 3: FIXING CHANNEL TYPES ===\n');

% obj already loaded in Step 2.5 (or load if Step 2.5 was skipped in a prior edit)
if ~exist('obj', 'var') || ~isa(obj, 'ecog_data_seeg')
    load(crunchedFile, 'obj');
end

for i = 1:numel(obj.elec_ch_type)
    obj.elec_ch_type{i} = lower(obj.elec_ch_type{i});
end

validSEEG = false(numel(obj.elec_ch_label), 1);
for i = 1:numel(obj.elec_ch_label)
    label = obj.elec_ch_label{i};
    if strcmp(obj.elec_ch_type{i}, 'seeg') && ~isempty(regexp(label, '\d', 'once'))
        validSEEG(i) = true;
    end
end

invalidChans = find(~validSEEG);
obj.elec_ch_prelim_deselect = union(obj.elec_ch_prelim_deselect, invalidChans);
obj.define_clean_channels();

fprintf('Total channels: %d | Valid SEEG: %d | Clean: %d | Excluded: %d\n', ...
    numel(obj.elec_ch), sum(validSEEG), numel(obj.elec_ch_clean), numel(obj.elec_ch_prelim_deselect));


%% ========================================================================
% STEP 3.5: ADD ANATOMY
%% ========================================================================
fprintf('\n=== STEP 3.5: ANATOMY ===\n');

if ~isfield(obj, 'anatomy') || isempty(obj.anatomy)
    try
        obj.add_anatomy(anatomyPath);
        fprintf('Anatomy loaded.\n');
    catch ME
        warning('Could not load anatomy: %s', ME.message);
    end
else
    fprintf('Anatomy already present.\n');
end


%% ========================================================================
% STEP 4: PREPROCESS (HIGH GAMMA EXTRACTION)
%% ========================================================================
fprintf('\n=== STEP 4: PREPROCESSING ===\n');

needPreproc = forceReprocess || ~isfield(obj.for_preproc, 'order') || isempty(obj.for_preproc.order);

if needPreproc
    % 1) Broadband referencing (engine order, delegated to @ecog_data_ieeg):
    %    highpass -> notch -> IED -> CAR -> Laplacian -> bipolar, with NO
    %    envelope extraction and NO downsampling yet.
    fprintf('Preprocessing order: %s\n', preprocOrder);
    obj.preprocess_signal('order', preprocOrder, ...
        'isPlotVisible', isPlotVisible, ...
        'doneVisualInspection', doneVisualInspection);

    % 2) High-gamma envelope via the NAPLAB filterbank.
    obj.extract_high_gamma('doNapLabFilterExtraction', true);

    % 3) Downsample the HG envelope.
    obj.downsample_signal('decimationFreq', decimationFreq);
else
    fprintf('Preprocessing already present; skipping.\n');
end

fprintf('Unipolar: [%d x %d] | Bipolar: [%d x %d] | Fs: %.1f Hz\n', ...
    size(obj.elec_data,1), size(obj.elec_data,2), ...
    size(obj.bip_elec_data,1), size(obj.bip_elec_data,2), obj.sample_freq);


%% ========================================================================
% STEP 5: MAP CONDITION NAMES
%% ========================================================================
fprintf('\n=== STEP 5: CONDITION MAPPING ===\n');

condMap = taskConfig.conditionMap;
fprintf('Original conditions: %s\n', strjoin(unique(obj.condition), ', '));
for i = 1:numel(obj.condition)
    if isKey(condMap, obj.condition{i})
        obj.condition{i} = condMap(obj.condition{i});
    end
end
fprintf('Mapped conditions:   %s\n', strjoin(unique(obj.condition), ', '));


%% ========================================================================
% STEP 6: SEGMENT INTO TRIALS
%% ========================================================================
fprintf('\n=== STEP 6: TRIAL SEGMENTATION ===\n');

if isempty(obj.trial_data)
    obj.make_trials();
    fprintf('Created %d trial segments.\n', numel(obj.trial_data));
else
    fprintf('Trials already present: %d\n', numel(obj.trial_data));
end

save(crunchedFile, 'obj', '-v7.3');


%% ========================================================================
% STEP 7: CREATE ANALYSIS OBJECT
%% ========================================================================
fprintf('\n=== STEP 7: ANALYSIS OBJECT ===\n');

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
fprintf('Output directory: %s\n', sn_obj.langloc_save_path);


%% ========================================================================
% STEP 7.5: SIGNIFICANCE + BASELINE Z-SCORE
%% ========================================================================
fprintf('\n=== STEP 7.5: SIGNIFICANCE + BASELINE Z-SCORE ===\n');

if isempty(sn_obj.stats) || ~isstruct(sn_obj.stats)
    sn_obj.stats = struct();
end

% This is the ieeg_pipeline feature chain that feeds the language channel
% selection, run on the analysis object:
%   extract_significant_channel -> extract_time_significance
%   -> 


%-> normalize_signal('z-score')
%
% The significance steps run on the (un-normalized) high-gamma envelope and
% store their results in sn_obj.stats. The baseline for both significance and
% normalization is anchored to probe_key = 1 (fixation row when present, else
% first event); default baseTimeRange [-0.5 0] is the pre-fixation baseline.
sn_obj.extract_significant_channel();
sn_obj.extract_time_significance();
sn_obj.extract_normalization_metrics();

% normalize_signal is inherited from the ieeg_pipeline engine and additionally
% Gaussian-smooths the high-gamma envelope (so make_trials below captures the
% smoothed, z-scored signal).
sn_obj.normalize_signal('normtype', 'z-score');
sn_obj.make_trials();

fprintf('Significance computed; baseline z-score (+ engine smoothing) applied. Trials rebuilt.\n');

% Optional sharp-artifact detection on the z-scored high-gamma envelope.
% Results are stored in sn_obj.stats.artifact_stats_unipolar/_bipolar.
if detectSharpArtifacts
    sn_obj.detect_sharp_artifacts();
end


%% ========================================================================
% STEP 8: LOCALIZE (source task) OR APPLY ROI (other tasks)
%% ========================================================================
fprintf('\n=== STEP 8: S vs N / ROI ===\n');

if ~applyROI
    fprintf('Running localization: %s vs %s\n', taskConfig.S_condition, taskConfig.N_condition);

    sn_obj.test_s_vs_n('words', taskConfig.testWords, ...
        'S_condition_flag', taskConfig.S_condition, ...
        'N_condition_flag', taskConfig.N_condition, ...
        'n_rep', 10000, ...
        'threshold', 0.05, ...
        'side', 'right', ...
        'use_odd_for_inference', useOddForInference, ...
        'do_plot', true);

    % --- Held-out effect sizes (even trials) ---
    if computeEffectSizes
        fprintf('Computing held-out HG effect sizes (S - N power diff + S/N correlation)\n');
        sn_obj.compute_hg_power_diff_s_vs_n('words', taskConfig.testWords, ...
            'S_condition_flag', taskConfig.S_condition, ...
            'N_condition_flag', taskConfig.N_condition, ...
            'use_odd_for_inference', useOddForInference);
        sn_obj.compute_hg_sn_corr('words', taskConfig.testWords, ...
            'S_condition_flag', taskConfig.S_condition, ...
            'N_condition_flag', taskConfig.N_condition);
    end

    % --- Word-level selection (per-word permutation + consecutiveness) ---
    if doWordwiseLangloc
        fprintf('Running word-wise langloc (min_consecutive=%d)\n', wordwiseMinConsec);
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

    % --- Word-boundaries cluster-based time-series permutation test ---
    if doWordBoundaries
        fprintf('Running word-boundaries time-series cluster permutation test\n');
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

    if useROIfromSource && strcmpi(taskType, roiSourceTask)
        save_roi_from_sn_obj(sn_obj, roiFile, roiSourceTask);
        fprintf('Saved ROI: %s\n', roiFile);
    end

else
    fprintf('Applying ROI from: %s\n', roiSourceTask);
    if ~exist(roiFile, 'file')
        error('ROI file not found. Run taskType=%s first.\nMissing: %s', roiSourceTask, roiFile);
    end
    tmp = load(roiFile, 'roi');
    sn_obj = apply_roi_to_sn_obj(sn_obj, tmp.roi);
    fprintf('ROI applied. Unipolar ROI channels: %d\n', sum(sn_obj.s_vs_n_sig.elec_data{1}));
end

if ~istable(sn_obj.s_vs_n_sig)
    error('sn_obj.s_vs_n_sig is missing or invalid.');
end


%% ========================================================================
% STEP 9: GENERATE PLOTS
%% ========================================================================
fprintf('\n=== STEP 9: PLOTTING ===\n');

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


%% ========================================================================
% STEP 10: SUMMARY STATISTICS
%% ========================================================================
fprintf('\n=== STEP 10: SUMMARY STATISTICS ===\n');
try
    summary = sn_obj.get_summary_statistics();
    disp(summary);
catch ME
    warning('Summary statistics failed (non-fatal): %s', ME.message);
end


%% ========================================================================
% STEP 11: CHANNEL LIST
%% ========================================================================
fprintf('\n=== STEP 11: CHANNEL LIST ===\n');

roiUni = find(sn_obj.s_vs_n_sig.elec_data{1});
fprintf('Significant unipolar channels: %d (%.1f%% of clean)\n', ...
    numel(roiUni), 100*numel(roiUni)/numel(sn_obj.elec_ch_clean));

if ~isempty(roiUni)
    fprintf('  Idx  Label        p-ratio\n');
    for i = 1:numel(roiUni)
        pval = sn_obj.s_vs_n_p_ratio.elec_data{1}(roiUni(i));
        fprintf('  %3d  %-12s %.4f\n', roiUni(i), sn_obj.elec_ch_label{roiUni(i)}, pval);
    end
end


%% ========================================================================
% STEP 12: ANATOMY PLOT (SUBJECT SPACE)
%% ========================================================================
fprintf('\n=== STEP 12: ANATOMY (SUBJECT SPACE) ===\n');

if isfield(sn_obj, 'anatomy') && isfield(sn_obj.anatomy, 'subject_space') ...
        && isfield(sn_obj.anatomy.subject_space, 'tala') ...
        && isfield(sn_obj.anatomy.subject_space.tala, 'electrodes')

    electrodes = sn_obj.anatomy.subject_space.tala.electrodes;
    figure('Position', [100 100 1200 800], 'Color', 'w');
    scatter3(electrodes(:,1), electrodes(:,2), electrodes(:,3), ...
        50, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k');
    hold on;

    if ~isempty(roiUni)
        mappedCells = sn_obj.anatomy.mapping(roiUni);
        ok = ~cellfun(@isempty, mappedCells);
        mappedIdx = cell2mat(mappedCells(ok));
        mappedIdx = mappedIdx(mappedIdx >= 1 & mappedIdx <= size(electrodes,1));
        mappedIdx = unique(mappedIdx(:));
        roiCoords = electrodes(mappedIdx, :);
        scatter3(roiCoords(:,1), roiCoords(:,2), roiCoords(:,3), ...
            150, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end

    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
    title(sprintf('%s - Significant Electrodes (%s)', params.SubjectName, taskType));
    legend({'All Electrodes', 'Significant'}, 'Location', 'best');
    grid on; axis equal; view(3); rotate3d on;

    saveas(gcf, fullfile(outputDir, sprintf('%s_anatomy_subject_3D.png', params.SubjectName)));
    fprintf('Saved subject-space anatomy plot.\n');
else
    warning('No subject-space anatomy found. Skipping Step 12.');
end


%% ========================================================================
% STEP 13: ANATOMY PLOT (MNI)
%% ========================================================================
fprintf('\n=== STEP 13: ANATOMY (MNI) ===\n');

try
    plot_sig_electrodes_anatomy(sn_obj, outputDir, ...
        'subjectName', params.SubjectName, ...
        'space', 'mni', ...
        'plotAllElectrodes', false, ...
        'signalType', 'bipolar', ...
        'doGIF', true, ...
        'angle', 270, ...
        'brainAlpha', 0.35, ...
        'sigColor', [1 0 0], ...
        'sigSize', 120);
    fprintf('Saved MNI anatomy plot.\n');
catch ME
    warning('MNI anatomy plotting failed (non-fatal): %s', ME.message);
end


%% ========================================================================
% STEP 14: SAVE GROUP-COMPATIBLE RESULTS (MIT_multi_single format)
%% ========================================================================
fprintf('\n=== STEP 14: SAVING GROUP RESULTS ===\n');

saveOpts = struct();
saveOpts.sourceFile            = crunchedFile;
saveOpts.includeEffectSizes    = computeEffectSizes;
saveOpts.includeWordwise       = doWordwiseLangloc;
saveOpts.includeWordBoundaries = doWordBoundaries;
saveOpts.useOddForInference    = useOddForInference;
saveOpts.overwrite             = true;
saveOpts.saveGammaPower        = true;
outFiles = save_mit_group_outputs(sn_obj, params.SubjectName, taskType, taskConfig, outputDir, saveOpts);
groupResultFile = outFiles.groupResultFile;


%% ========================================================================
% STEP 14.5: SYNC obj FROM sn_obj
%% ========================================================================
fprintf('\n=== STEP 14.5: SYNCING obj FROM sn_obj ===\n');

obj = sync_obj_from_sn_obj(obj, sn_obj);

if isstruct(obj.stats) && ~isempty(fieldnames(obj.stats))
    fprintf('obj.stats fields: %s\n', strjoin(fieldnames(obj.stats), ', '));
else
    warning('obj.stats is still empty after sync — Step 7.5 may not have completed.');
end

save(crunchedFile, 'obj', 'sn_obj', '-v7.3');
fprintf('Saved synced obj and sn_obj: %s\n', crunchedFile);


%% ========================================================================
% STEP 15: LANGLOC PDF REPORT (MITSWJNTask only)
%% ========================================================================
if generateLanglocReport && strcmpi(taskType, 'MITSWJNTask')
    fprintf('\n=== STEP 15: LANGLOC PDF REPORT ===\n');
    try
        run_langloc_report(obj, taskType, taskConfig, outputDir, ...
            params.SubjectName, useLanglocReportV2);
    catch ME
        warning('Langloc report generation failed (non-fatal): %s', ME.message);
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        if contains(ME.message, 'mlreportgen') || contains(ME.identifier, 'MATLAB:UndefinedFunction')
            fprintf(['  Tip: install MATLAB Report Generator, or set ' ...
                'generateLanglocReport = false.\n']);
        end
    end
else
    if generateLanglocReport && ~strcmpi(taskType, 'MITSWJNTask')
        fprintf('\n=== STEP 15: LANGLOC PDF REPORT (skipped; MITSWJNTask only) ===\n');
    else
        fprintf('\n=== STEP 15: LANGLOC PDF REPORT (skipped) ===\n');
    end
end

%% ========================================================================
% DONE
%% ========================================================================
fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');
fprintf('Task: %s\n', taskType);
if applyROI
    fprintf('ROI applied from: %s\n', roiSourceTask);
else
    fprintf('Localization performed on: %s\n', taskType);
end
fprintf('Results: %s\n\n', sn_obj.langloc_save_path);
