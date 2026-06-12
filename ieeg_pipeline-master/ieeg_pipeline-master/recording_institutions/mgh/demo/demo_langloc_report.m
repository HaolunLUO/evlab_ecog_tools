%% MGH iEEG — LangLoc Audio: Full Pipeline Demo with Report Generation
%
% This demo walks through the complete LangLoc Audio preprocessing workflow
% and generates a PDF report using generateReportLangloc_v2.
%
% The report includes:
%   - Experiment metadata table (subject, trial counts, electrode counts)
%   - Behavioral performance (accuracy, RT distributions by condition)
%   - Audio duration distributions (Natus vs behavioral)
%   - LangLoc-responsive electrode plots (Dysoc-style, per-word boundaries)
%   - High-gamma time series across all channels
%   - Summary of significant channels (unipolar + bipolar)
%
% Prerequisites:
%   1. Run langloc_audio_template.m first to produce the crunched .mat file.
%   2. Ensure utils/kumar_ieeg_utils/ is on the MATLAB path (it contains
%      generateReportLangloc_v2 and its helper functions).
%   3. MATLAB Report Generator toolbox must be installed.
%
% Kumar Duraivel, EvLab @ MIT

%% EDIT: Point to the crunched object produced by langloc_audio_template.m
SUBJECT   = 'sub-XXXX';
SESSION   = 'LangLocAudio';
DATAPATH  = '/path/to/LangLoc/data';
order     = 'defaultSEEGorBOTHBroadBand';

crunched_path = fullfile(DATAPATH, 'derivatives', SUBJECT, 'preproc', 'crunched', ...
    [SUBJECT '_' SESSION '_crunched_' order '.mat']);

%% Add utilities to path
addpath(genpath('../../../utils'));  % ieeg_pipeline/utils (includes kumar_ieeg_utils)

%% Load the preprocessed object
fprintf('Loading crunched object: %s\n', crunched_path);
load(crunched_path, 'obj');

%% (Optional) Re-run post-processing steps if not saved in object
% Uncomment and run if obj.stats is missing the required fields:
%
% obj.extract_high_gamma('doNapLabFilterExtraction', true);
% obj.downsample_signal('decimationFreq', 100);
% obj.extract_significant_channel();
% obj.extract_time_significance();
% obj.extract_normalization_metrics();
% obj.normalize_signal('normtype', 'z-score');

%% Generate the LangLoc report
% The report PDF will be saved to:
%   derivatives/SUBJECT/preproc/crunched/SUBJECT_SESSION.pdf
%
% generateReportLangloc_v2 produces:
%   Page 1  — Experiment metadata table
%   Page 2  — Behavioral accuracy + RT distributions (sentence vs nonword)
%   Page 3  — Audio duration distributions
%   Pages 4+ — LangLoc-responsive electrode plots (split-half scatter,
%               condition bar chart, word-boundary time series)
%   Pages N+ — High-gamma time series (10 channels per page)
%   Final   — Summary of significant channels

reportName = [obj.subject '_' obj.experiment];
fprintf('Generating report: %s.pdf\n', reportName);
generateReportLangloc_v2(obj, reportName);
fprintf('Report saved to: %s\n', fullfile(obj.crunched_file_path, [reportName '.pdf']));
