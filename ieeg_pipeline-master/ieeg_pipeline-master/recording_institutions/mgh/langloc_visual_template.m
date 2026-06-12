%% MGH iEEG — Language Localizer (Visual) Preprocessing Template
%
% Template for preprocessing LangLoc Visual data recorded at Massachusetts
% General Hospital using the Natus EEG system.
%
% Full pipeline (subject-specific scripts, task utilities) lives in the
% companion repository:
%   https://github.com/YOURUSERNAME/MGH_IEEG_preproc
%
% This template covers the core preprocessing steps. For the complete
% workflow including report generation, refer to the langloc/ folder of
% the above repository.
%
% Dependencies (add to MATLAB path before running):
%   - @ecog_data class (this repository)
%   - utils/ folder (this repository)
%   - langloc_utils/ folder from MGH_IEEG_preproc
%   - utils/ folder from MGH_IEEG_preproc
%
% Kumar Duraivel, EvLab @ MIT

clear all
close all

%% EDIT: Subject / session info
SUBJECT  = 'sub-XXXX';          % BIDS subject ID (e.g., 'sub-EM1220')
SESSION  = 'LangLocVideo';
MODALITY = 'video';

%% EDIT: Data paths
DATAPATH = '/path/to/LangLoc/data';

%% EDIT: Path to MGH_IEEG_preproc utilities
MGH_PREPROC_REPO = '/path/to/MGH_IEEG_preproc';

%% Add utilities to path
addpath(genpath(fullfile(MGH_PREPROC_REPO, 'langloc', 'langloc_utils')));
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
% LangLoc Visual uses TRIG; bit1 = trial onset
chan_insp = {'TRIG'};
DC_files  = cell2mat(cellfun(@(x) find(strcmp(hdr.label, x)), chan_insp, 'uni', false));
TrigMat1  = record(DC_files, :)';

filteredEventTimes = processAndPlotTriggerEventsLangLocVisual(TrigMat1);
% filteredEventTimes{1}  = trial onset times
% filteredEventTimes{10} = trial offset times

%% Load behavioral CSV files
d_events = dir(fullfile(PATH_EVENTS, '*.csv'));
if isempty(d_events)
    d_events = dir(fullfile(PATH_SESSION, 'task', '*.csv'));
end

% EDIT: select which runs to include
task_files_to_pick = 1:length(d_events);
d_events = d_events(task_files_to_pick);

[events_table] = extract_behavioral_events_for_langloc_visual( ...
    'behavior_files', d_events, 'sampling', unique(sampling_frequency));
% EDIT: update expected trial count
assert(size(events_table, 1) == 80, ...
    'Expected 80 trials; got %d. Check d_events.', size(events_table, 1));

%% Align Natus triggers with behavioral data
time2save = filteredEventTimes{1}(1) - 30*sampling_frequency(1) : ...
            filteredEventTimes{1}(end) + 30*sampling_frequency(1);
timeStart = time2save(1);

natusTrialStart = (filteredEventTimes{1}  - timeStart) ./ sampling_frequency(1);
natusTrialEnd   = (filteredEventTimes{10} - timeStart) ./ sampling_frequency(1);

% Sanity check: compare ISI between Natus and behavioral
isiNatus = [diff(natusTrialStart(1:40)); diff(natusTrialStart(41:80))];
isiBeh   = [diff(events_table.actual_onset(1:40)); diff(events_table.actual_onset(41:80))];
figure; histogram(isiNatus - isiBeh, 50);
xlabel('Discrepancy (s)'); ylabel('Trials');

%% Build trial timing
[trial_timing] = get_timing_for_LangLocVideo(filteredEventTimes, events_table, timeStart);

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

obj.events_table = obj.for_preproc.event_table;
obj.condition    = cellfun(@(x) replace(x, {'S','N'}, {'sentence','nonword'}), ...
    obj.for_preproc.event_table.condition, 'UniformOutput', false);
obj.session      = obj.for_preproc.event_table.list;
obj.trial_timing = trial_timing(:, 1);

if ~isfolder(save_path), mkdir(save_path); end
save(fullfile(save_path, save_filename), 'obj', '-v7.3');

%% High-gamma extraction, downsampling, normalization
obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 100);
obj.extract_significant_channel();
obj.extract_time_significance();
obj.extract_normalization_metrics();
obj.normalize_signal('normtype', 'z-score');
