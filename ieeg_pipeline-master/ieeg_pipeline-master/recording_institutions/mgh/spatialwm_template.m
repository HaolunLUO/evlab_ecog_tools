%% MGH iEEG — Spatial Working Memory Preprocessing Template
%
% Template for preprocessing Spatial WM data recorded at Massachusetts
% General Hospital using the Natus EEG system.
%
% Full pipeline (subject-specific scripts, task utilities) lives in the
% companion repository:
%   https://github.com/YOURUSERNAME/MGH_IEEG_preproc
%
% Dependencies (add to MATLAB path before running):
%   - @ecog_data class (this repository)
%   - utils/ folder (this repository)
%   - spatialWM_utils/ folder from MGH_IEEG_preproc
%   - utils/ folder from MGH_IEEG_preproc
%
% Kumar Duraivel, EvLab @ MIT

clear all
close all

%% EDIT: Subject / session info
SUBJECT = 'sub-XXXX';           % BIDS subject ID (e.g., 'sub-EM1233')
SESSION = 'SpatialWM';

%% EDIT: Data paths
DATAPATH = '/path/to/SpatialWM/data';

%% EDIT: Path to MGH_IEEG_preproc utilities
MGH_PREPROC_REPO = '/path/to/MGH_IEEG_preproc';

%% Add utilities to path
addpath(genpath(fullfile(MGH_PREPROC_REPO, 'spatialWM', 'spatialWM_utils')));
addpath(genpath(fullfile(MGH_PREPROC_REPO, 'utils')));
addpath(genpath('../../utils'));

%% Define data paths (BIDS-like structure)
PATH_DATA    = fullfile(DATAPATH, 'raw_data', SUBJECT);
PATH_SESSION = fullfile(PATH_DATA, ['ses-' SESSION]);
PATH_EDF     = fullfile(PATH_SESSION, 'natus');
PATH_EVENTS  = fullfile(PATH_SESSION, 'tasks');
PATH_DER     = fullfile(DATAPATH, 'derivatives');
PATH_ANNOT   = fullfile(PATH_DER, SUBJECT, 'annot');
PATH_SAVE    = fullfile(PATH_DER, SUBJECT, 'preproc');

if ~exist(PATH_ANNOT, 'dir'), mkdir(PATH_ANNOT); end
if ~exist(PATH_SAVE,  'dir'), mkdir(PATH_SAVE);  end

%% Load neural data (Natus EDF)
edflist = dir(fullfile(PATH_EDF, '*.EDF'));
edfname = edflist(1).name;

[hdr, record] = edfread(fullfile(PATH_EDF, edfname));
info = edfinfo(fullfile(PATH_EDF, edfname));
sampling_frequency = hdr.frequency;

%% Read or create BIDS channels table
channels_table = create_channels_table_bids(info, PATH_ANNOT, SUBJECT, SESSION);

%% Parse Natus trigger channel
chan_insp = {'TRIG'};
DC_files  = cell2mat(cellfun(@(x) find(strcmp(hdr.label, x)), chan_insp, 'uni', false));
TrigMat1  = record(DC_files, :)';

filteredEventTimes = processAndPlotTriggerEventsSpatialWM(TrigMat1);

%% Load behavioral CSV files
d_events = dir(fullfile(PATH_EVENTS, '*.csv'));
if isempty(d_events)
    d_events = dir(fullfile(PATH_SESSION, 'task', '*.csv'));
end

% EDIT: select which runs to include
task_files_to_pick = 1:length(d_events);
d_events = d_events(task_files_to_pick);

[events_table] = extract_behavioral_events_for_spatialWM( ...
    'behavior_files', d_events, ...
    'sampling', unique(sampling_frequency), ...
    'filteredEventTimes', filteredEventTimes);
% EDIT: update expected trial count (typically 72 per session)
assert(numel(events_table) == 72, ...
    'Expected 72 trials; got %d. Check d_events.', numel(events_table));

%% Align Natus triggers with behavioral data
time2save = (events_table(1).trial_onset - 15) * sampling_frequency(1) : ...
            (events_table(end).trial_onset + 15) * sampling_frequency(1);
timeStart = time2save(1);

natusTrialStart = [events_table.trial_onset];
behTimingOnset  = [events_table.beh_onset];

isiNatus = [diff(natusTrialStart(1:36)); diff(natusTrialStart(37:72))];
isiBeh   = [diff(behTimingOnset(1:36));  diff(behTimingOnset(37:72))];
figure; histogram(isiNatus - isiBeh, 50);
xlabel('Discrepancy (s)'); ylabel('Trials');

%% Build trial timing
[trial_timing, events_table] = make_trials_swm( ...
    events_table, unique(sampling_frequency), timeStart);

%% Write ecog_data object
subject    = SUBJECT;
experiment = SESSION;
order      = 'defaultSEEGorBOTHBroadBand';

save_filename     = [subject '_' experiment '_crunched_' order '.mat'];
save_path         = fullfile(PATH_SAVE, 'crunched');
formattedDateTime = datestr(datetime('now'), 'yyyymmdd_HHMM');
log_filename      = fullfile(PATH_SAVE, 'logs', ...
    [subject '_' experiment '_' order '_' formattedDateTime '.txt']);

d_files   = fullfile(PATH_EDF, edfname);
ch_labels = channels_table.name;
ch_type   = string(channels_table.type);
ch_select = find(contains(ch_type, 'seeg'));

for_preproc = struct;
for_preproc.elec_data_raw    = single(record(ch_select, time2save));
for_preproc.event_table      = events_table;
for_preproc.stitch_index_raw = 1;
for_preproc.sample_freq_raw  = sampling_frequency(1);
for_preproc.log_file_name    = log_filename;
for_preproc.decimation_freq  = sampling_frequency(1) / 4;

obj = ecog_data(for_preproc, subject, experiment, save_filename, save_path, ...
    d_files, PATH_EDF, ch_labels(ch_select), 1:length(ch_select), [], ch_type(ch_select));

obj.preprocess_signal('order', order, 'isPlotVisible', false, 'doneVisualInspection', false);

obj.events_table = struct2table(obj.for_preproc.event_table);
obj.condition    = [obj.for_preproc.event_table.condition];
obj.session      = [obj.for_preproc.event_table.session];
obj.trial_timing = trial_timing';

if ~isfolder(save_path), mkdir(save_path); end
save(fullfile(save_path, save_filename), 'obj', '-v7.3');

%% High-gamma extraction, downsampling, normalization
obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 200);
obj.extract_significant_channel();
obj.extract_time_significance();
obj.extract_normalization_metrics();
obj.normalize_signal('normtype', 'z-score');
