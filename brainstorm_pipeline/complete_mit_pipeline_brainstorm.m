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
%   2) Preprocess: highpass, notch, CAR/bipolar, high-gamma extraction
%   3) (roiSourceTask) Run S vs N localization and save ROI
%      (other tasks)   Load ROI and apply to current task
%   4) Generate timecourse + barplots + anatomy plots
%   5) Save group-compatible result .mat file
%
% RECOMMENDED USAGE
%   A) Run once with taskType='MITSWJNTask' -> defines and saves ROI
%   B) Run with taskType='WM'        -> plots MITSWJNTask ROI on WM data
%   C) Run with taskType='Auditory'  -> plots MITSWJNTask ROI on Auditory
%
% MATLAB PATH REQUIREMENTS
%   addpath(genpath('<repo>/brainstorm_pipeline'))  % MUST come FIRST
%   addpath(genpath('<repo>'))                      % rest of evlab_ecog_tools
%   brainstorm_to_mit_crunched_new.m must be on path (provided separately)
%   JancaCodePapers/ must be on path (in this repo) for IED removal
%% ========================================================================

clear; clc; close all;
addpath(genpath('F:\seeg\luohong\analysisEV\v2_piepeline\evlab_ecog_tools\brainstorm_pipeline'))  % MUST come FIRST
addpath(genpath('F:\seeg\luohong\analysisEV\v2_piepeline\evlab_ecog_tools\')) 
%% ========================================================================
% USER SETTINGS
%% ========================================================================
taskType = 'MITSWJNTask';  % 'MITSWJNTask' | 'WM' | 'Auditory' | 'MITLangloc'

% --- Cross-task ROI ---
useROIfromSource = false;
roiSourceTask    = 'MITSWJNTask';

% --- Re-run control ---
forceRebuildCrunched = false;
forceReprocess       = true;

% --- UI control ---
isPlotVisible        = true;
doneVisualInspection = true;

% --- Paths (EDIT THESE) ---
workingDir  = 'F:\seeg\luohong\analysisEV';
anatomyPath = fullfile(workingDir, 'anatomy\');

% --- Subject/protocol ---
params = struct();
params.SubjectName  = 'Subject11';
params.ProtocolName = 'analysis';
params.outputPath   = workingDir;
params.taskType     = taskType;

% --- Data files (EDIT PER TASK) ---
switch taskType
    case 'MITSWJNTask'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject11\DA01008R\data_block001_02.mat'
        };
    case 'WM'
        allDataFiles = {
            'F:\seeg\analysis\data\Subject01\DA0011IJ\data_block001_07.mat'
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
switch taskType
    case 'MITLangloc'
        taskConfig = struct();
        taskConfig.nWordPositions = 12;
        taskConfig.wordDuration   = 0.45;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map(...
            {'sentence', 'word lists', 'Jabberwocky', 'non-word lists'}, ...
            {'Sentences', 'Word-lists', 'Jabberwocky', 'Nonword-lists'});
        taskConfig.S_condition = 'Sentences';
        taskConfig.N_condition = 'Nonword-lists';
        taskConfig.W_condition = 'Word-lists';
        taskConfig.J_condition = 'Jabberwocky';
        taskConfig.testWords   = 1:12;
        taskConfig.subAverage  = true;

    case 'MITSWJNTask'
        taskConfig = struct();
        taskConfig.nWordPositions = 8;
        taskConfig.wordDuration   = 0.70;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map(...
            {'sentence', 'word lists', 'Jabberwocky', 'non-word lists'}, ...
            {'SENTENCES', 'WORDS', 'JABBERWOCKY', 'NONWORDS'});
        taskConfig.S_condition = 'SENTENCES';
        taskConfig.N_condition = 'NONWORDS';
        taskConfig.W_condition = 'WORDS';
        taskConfig.J_condition = 'JABBERWOCKY';
        taskConfig.testWords   = 1:8;
        taskConfig.subAverage  = false;

    case 'Auditory'
        taskConfig = struct();
        taskConfig.nWordPositions = 1;
        taskConfig.wordDuration   = 18;
        taskConfig.eventPattern   = '%s';
        taskConfig.conditionMap = containers.Map({'Intact', 'Degraded'}, {'Intact', 'Degraded'});
        taskConfig.S_condition = 'Intact';
        taskConfig.N_condition = 'Degraded';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1;
        taskConfig.subAverage  = true;

    case 'WM'
        taskConfig = struct();
        taskConfig.nWordPositions = 5;
        taskConfig.wordDuration   = 1;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map({'hard', 'easy'}, {'Hard', 'Easy'});
        taskConfig.S_condition = 'Hard';
        taskConfig.N_condition = 'Easy';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1:5;
        taskConfig.subAverage  = true;

    otherwise
        error('Unknown task type: %s', taskType);
end

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
    fprintf('Found existing crunched file. Skipping conversion.\n');
else
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
% STEP 3: FIX CHANNEL TYPES + CLEAN CHANNELS
%% ========================================================================
fprintf('\n=== STEP 3: FIXING CHANNEL TYPES ===\n');

load(crunchedFile, 'obj');

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
    obj.preprocess_signal('order', 'defaultSEEGorBOTH', ...
        'isPlotVisible', isPlotVisible, ...
        'doneVisualInspection', doneVisualInspection);
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

sn_obj = ecog_sn_data(...
    outputDir, ...
    crunchedFile, ...
    workingDir, ...
    'ecog_data.m', ...
    workingDir);

sn_obj.anatomy = obj.anatomy;

if ~exist(sn_obj.langloc_save_path, 'dir')
    mkdir(sn_obj.langloc_save_path);
end
fprintf('Output directory: %s\n', sn_obj.langloc_save_path);


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
        'do_plot', true);

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
% STEP 14: SAVE GROUP-COMPATIBLE RESULTS
%% ========================================================================
fprintf('\n=== STEP 14: SAVING GROUP RESULTS ===\n');

groupResultFile = fullfile(outputDir, sprintf('%s_%s_groupResult.mat', params.SubjectName, taskType));

groupResult = struct();
groupResult.subject  = params.SubjectName;
groupResult.taskType = taskType;

groupResult.sig_uni     = logical(sn_obj.s_vs_n_sig.elec_data{1});
groupResult.p_ratio_uni = sn_obj.s_vs_n_p_ratio.elec_data{1};
groupResult.elec_labels = sn_obj.elec_ch_label;
groupResult.elec_valid  = sn_obj.elec_ch_valid;
groupResult.elec_clean  = sn_obj.elec_ch_clean;
groupResult.nClean      = numel(sn_obj.elec_ch_clean);
groupResult.nSig        = sum(groupResult.sig_uni);

if ~isempty(sn_obj.bip_elec_data) && ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
    groupResult.sig_bip     = logical(sn_obj.s_vs_n_sig.bip_elec_data{1});
    groupResult.p_ratio_bip = sn_obj.s_vs_n_p_ratio.bip_elec_data{1};
    groupResult.bip_labels  = sn_obj.bip_ch_label;
    groupResult.nSigBip     = sum(groupResult.sig_bip);
else
    groupResult.sig_bip  = [];
    groupResult.bip_labels = {};
    groupResult.nSigBip  = 0;
end

if isfield(sn_obj, 'anatomy') && isfield(sn_obj.anatomy, 'mni_space') ...
        && isfield(sn_obj.anatomy.mni_space, 'tala')
    allMNI = sn_obj.anatomy.mni_space.tala.electrodes;
    groupResult.mni_coords = allMNI;
    sigIdx = find(groupResult.sig_uni);
    mappedCells = sn_obj.anatomy.mapping(sigIdx);
    ok = ~cellfun(@isempty, mappedCells);
    sigAnatomyIdx = cell2mat(mappedCells(ok));
    sigAnatomyIdx = sigAnatomyIdx(sigAnatomyIdx >= 1 & sigAnatomyIdx <= size(allMNI,1));
    sigMask = false(size(allMNI,1), 1);
    sigMask(sigAnatomyIdx) = true;
    groupResult.mni_sig_mask = sigMask;
end

groupResult.sample_freq = sn_obj.sample_freq;
groupResult.nWords      = numel(taskConfig.testWords);

try
    [S_tc, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, 'condition', taskConfig.S_condition, 'signalType', 'unipolar');
    [N_tc, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, 'condition', taskConfig.N_condition, 'signalType', 'unipolar');
    groupResult.S_timecourse_mean = mean(S_tc, 1);
    groupResult.S_timecourse_sem  = std(S_tc, [], 1) / sqrt(size(S_tc,1));
    groupResult.N_timecourse_mean = mean(N_tc, 1);
    groupResult.N_timecourse_sem  = std(N_tc, [], 1) / sqrt(size(N_tc,1));
    groupResult.nSigElecs_tc      = size(S_tc, 1);
catch ME
    warning('Monopolar timecourses failed: %s', ME.message);
end

if ~isempty(sn_obj.bip_elec_data) && ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames) ...
        && sum(sn_obj.s_vs_n_sig.bip_elec_data{1}) > 0
    try
        [S_tc_bip, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, 'condition', taskConfig.S_condition, 'signalType', 'bipolar');
        [N_tc_bip, ~] = sn_obj.get_timecourses('words', taskConfig.testWords, 'condition', taskConfig.N_condition, 'signalType', 'bipolar');
        groupResult.S_timecourse_mean_bip = mean(S_tc_bip, 1);
        groupResult.S_timecourse_sem_bip  = std(S_tc_bip, [], 1) / sqrt(size(S_tc_bip,1));
        groupResult.N_timecourse_mean_bip = mean(N_tc_bip, 1);
        groupResult.N_timecourse_sem_bip  = std(N_tc_bip, [], 1) / sqrt(size(N_tc_bip,1));
        groupResult.nSigElecs_tc_bip      = size(S_tc_bip, 1);
    catch ME
        warning('Bipolar timecourses failed: %s', ME.message);
    end
end

save(groupResultFile, 'groupResult', '-v7.3');
fprintf('Saved group result: %s\n', groupResultFile);

save(crunchedFile, 'sn_obj', '-append');

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
