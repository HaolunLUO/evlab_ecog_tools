%% MGH iEEG — Speech Language Localizer Preprocessing Template
%
% Template for preprocessing Speech LangLoc data recorded at Massachusetts
% General Hospital using the Natus EEG system.
%
% Full pipeline (subject-specific scripts, task utilities) lives in the
% companion repository:
%   https://github.com/YOURUSERNAME/MGH_IEEG_preproc
%
% Paradigm: 4 runs x 36 trials = 144 total trials.
%   Conditions per run: sentence (12), nonword (12), quilt (12).
%   Passive listening — no behavioral probe.
%   Trigger bits: bit2 = audio onset, bit6 = audio offset, bit7 = space press.
%
% Dependencies (add to MATLAB path before running):
%   - @ecog_data class (this repository)
%   - utils/ folder (this repository)
%   - speechlangloc_utils/ folder from MGH_IEEG_preproc
%   - utils/ folder from MGH_IEEG_preproc
%
% Kumar Duraivel, EvLab @ MIT

clear all
close all

%% EDIT: Subject / session info
SUBJECT = 'sub-XXXX';           % BIDS subject ID (e.g., 'sub-EM1367')
SESSION = 'SpeechLangLocAudio';

%% EDIT: Data paths (raw data not yet moved to project BIDS tree)
PATH_RAW    = '/path/to/sub-XXXX_ses-SpeechLangLocAudio_iEEG';
PATH_EDF    = fullfile(PATH_RAW, 'natus');
PATH_EVENTS = fullfile(PATH_RAW, 'task');

%% EDIT: Derivatives output location
PATH_DER   = '/path/to/speechlangloc/derivatives';
PATH_ANNOT = fullfile(PATH_DER, SUBJECT, 'annot');
PATH_SAVE  = fullfile(PATH_DER, SUBJECT, 'preproc');

%% EDIT: Path to MGH_IEEG_preproc utilities
MGH_PREPROC_REPO = '/path/to/MGH_IEEG_preproc';

%% Add utilities to path
addpath(genpath(fullfile(MGH_PREPROC_REPO, 'speechLangLoc', 'speechlangloc_utils')));
addpath(genpath(fullfile(MGH_PREPROC_REPO, 'utils')));
addpath(genpath('../../utils'));

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
% Bit assignments:
%   bit1 = start_expt  bit2 = audio onset   bit3 = end_expt
%   bit4 = event       bit6 = audio offset  bit7 = space press
chan_insp = {'TRIG'};
DC_files  = cell2mat(cellfun(@(x) find(strcmp(hdr.label, x)), chan_insp, 'uni', false));
TrigMat1  = record(DC_files, :)';

filteredEventTimes = processAndPlotTriggerEventsSpeechLangLocAudio(TrigMat1);
% filteredEventTimes{2} should contain N_trials rising edges (audio onset)
% filteredEventTimes{6} should contain N_trials rising edges (audio offset)

%% Load behavioral CSV files
d_events = dir(fullfile(PATH_EVENTS, '*.csv'));
% EDIT: select which runs to include
task_files_to_pick = 1:length(d_events);
d_events = d_events(task_files_to_pick);

events_table = extract_behavioral_events_for_speechlangloc('behavior_files', d_events);
% EDIT: update expected trial count (4 runs x 36 trials = 144)
assert(size(events_table, 1) == 144, ...
    'Expected 144 trial rows; got %d. Check d_events.', size(events_table, 1));

%% Align Natus triggers with behavioral data
trialId   = 2;  % filteredEventTimes index for audio onset
time2save = filteredEventTimes{trialId}(1)  - 30*sampling_frequency(1) : ...
            filteredEventTimes{trialId}(end) + 30*sampling_frequency(1);
timeStart = time2save(1);

natusAudioStart = (filteredEventTimes{trialId} - timeStart) ./ sampling_frequency(1);
natusAudioEnd   = (filteredEventTimes{10}       - timeStart) ./ sampling_frequency(1);

assert(length(natusAudioStart) == size(events_table, 1));
assert(length(natusAudioEnd)   == size(events_table, 1));

audioDurNatus = natusAudioEnd - natusAudioStart;
audioDurBeh   = events_table.trial_offset - events_table.trial_onset;
figure; scatter(audioDurNatus, audioDurBeh, 20, 'black', 'filled');
figure; histogram(audioDurNatus - audioDurBeh, 50);
xlabel('Discrepancy (s)'); ylabel('Trials');

events_table.trial_onset_natus  = natusAudioStart - 0.2;
events_table.audio_onset_natus  = natusAudioStart;
events_table.audio_ended_natus  = natusAudioEnd;
events_table.trial_ended_natus  = natusAudioEnd + 0.2;

dt_post_ti = events_table.post_ti_onset - events_table.trial_onset;
events_table.post_ti_onset_natus = events_table.audio_onset_natus + dt_post_ti;

events_table.response_onset_natus = NaN(height(events_table), 1);
has_resp = ~isnan(events_table.response_RT);
events_table.response_onset_natus(has_resp) = ...
    events_table.post_ti_onset_natus(has_resp) + events_table.response_RT(has_resp);

%% Build trial timing
trial_timing = make_trials_mit_speechlangloc(events_table, sampling_frequency(1));

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
obj.condition    = obj.for_preproc.event_table.condition;
obj.session      = obj.for_preproc.event_table.run;
obj.trial_timing = trial_timing;

if ~isfolder(save_path), mkdir(save_path); end
save(fullfile(save_path, save_filename), 'obj', '-v7.3');

%% High-gamma extraction, downsampling, normalization
obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 200);
obj.extract_normalization_metrics();
obj.normalize_signal('normtype', 'z-score');
