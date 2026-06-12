# MGH iEEG Preprocessing Pipeline

Preprocessing pipeline for ECoG/SEEG data recorded at **Massachusetts General Hospital (MGH)** under the care of **Mark Richardson**. Data is acquired using the Natus EEG system and stored as EDF files with a BIDS-adjacent folder structure.

---

## Contact

| Person | Role | Email |
|---|---|---|
| **Kumar Duraivel** | Primary maintainer | dsuseendar@gmail.com |
| **Colton Casto** | Original pipeline author | ccasto@mit.edu |

---

## Companion Repository

The full MGH pipeline — including subject-specific preprocessing scripts, audio stimulus alignment files, and task utility functions — lives in a separate repository:

> **[MGH_IEEG_preproc](https://github.com/YOURUSERNAME/MGH_IEEG_preproc)**

The templates in this folder are parameterized starting points. For the complete workflow with audio alignment, visual inspection utilities, and detailed per-subject notes, refer to the companion repository.

---

## Workflow Overview

Each task follows the same five-stage pattern:

```
1. Load raw EDF (Natus)
       ↓
2. Parse trigger channel → filteredEventTimes
       ↓
3. Load behavioral CSVs → events_table
       ↓
4. Align Natus triggers with behavioral data (sanity-check plots)
       ↓
5. Build ecog_data object → preprocess → extract HG → normalize → save
```

---

## Data Structure

Data follows a BIDS-like layout:

```
DATAPATH/
├── raw_data/
│   └── sub-XXXX/
│       └── ses-SESSION/
│           ├── natus/           *.EDF  (Natus recording)
│           └── tasks/           *.csv  (behavioral log per run)
└── derivatives/
    └── sub-XXXX/
        ├── annot/               *_channels.tsv  (BIDS channel table)
        └── preproc/
            ├── crunched/        *_crunched_ORDER.mat
            └── logs/            preprocessing logs
```

Subject IDs follow the BIDS convention: `sub-EMXXXX` (e.g., `sub-EM1296`).

---

## Channels Table

`create_channels_table_bids()` auto-generates a TSV from the EDF header on first run and saves it to `derivatives/sub-XXXX/annot/`. On subsequent runs it loads the existing file, allowing manual edits (e.g., correcting channel types) to be preserved.

Channel type labels used:
- `seeg` — intracranial SEEG contacts
- `eeg` — scalp EEG channels
- `eog`, `ecg`, `emg` — physiological channels
- `MISC`, `OTHER` — other channels (e.g., trigger, DC)

Only channels typed `seeg` are passed to the `ecog_data` constructor.

---

## Trigger Parsing

Each task has a dedicated trigger parser in the MGH_IEEG_preproc `utils/` folder:

| Task | Parser | Trigger channel |
|---|---|---|
| LangLoc Audio | `processAndPlotTriggerEventsLangLocAudio` | `DC1` |
| LangLoc Visual | `processAndPlotTriggerEventsLangLocVisual` | `TRIG` |
| Spatial WM | `processAndPlotTriggerEventsSpatialWM` | `TRIG` |
| Sentences | `processAndPlotTriggerEventsLangLocAudio` | `TRIG` + `DC1` (microphone) |
| Speech LangLoc | `processAndPlotTriggerEventsSpeechLangLocAudio` | `TRIG` |

The parsers return `filteredEventTimes`, a cell array where each cell contains the sample indices of one event type (audio onset, audio offset, probe, etc.). Bit assignments are documented at the top of each template.

---

## Task Templates

Each template covers the full preprocessing pipeline for one task. Edit the `% EDIT:` variables at the top before running.

| Template | Task | Conditions | Typical trial count |
|---|---|---|---|
| `langloc_audio_template.m` | Language Localizer (Audio) | sentence, nonword | 80 (2 runs × 40) |
| `langloc_visual_template.m` | Language Localizer (Visual) | sentence, nonword | 80 (2 runs × 40) |
| `spatialwm_template.m` | Spatial Working Memory | various | 72 (2 runs × 36) |
| `sentences_template.m` | MIT Sentences | sentence, nonword | 120 (3 runs × 40) |
| `speechlangloc_template.m` | Speech Language Localizer | sentence, nonword, quilt | 144 (4 runs × 36) |

### Variables to edit in every template

```matlab
SUBJECT          = 'sub-XXXX';           % BIDS subject ID
DATAPATH         = '/path/to/data';       % Root data directory
MGH_PREPROC_REPO = '/path/to/MGH_IEEG_preproc';  % Companion repo

% Also review:
task_files_to_pick = 1:length(d_events);  % Exclude aborted/incomplete runs
assert(size(events_table,1) == N, ...);   % Update expected trial count
```

---

## Preprocessing Order

All MGH templates use `'defaultSEEGorBOTHBroadBand'`, which applies broadband preprocessing without inline high-gamma extraction. High-gamma is extracted separately after saving:

```matlab
obj.preprocess_signal('order', 'defaultSEEGorBOTHBroadBand');
% ... set trial_timing, condition, session ...
save(fullfile(save_path, save_filename), 'obj', '-v7.3');

obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 100);   % 200 for SpatialWM / SpeechLangLoc
obj.extract_normalization_metrics();
obj.normalize_signal('normtype', 'z-score');
```

---

## Demo

`demo/demo_langloc_report.m` walks through loading a crunched LangLoc Audio object and generating a full PDF report using `generateReportLangloc_v2` (in `utils/kumar_ieeg_utils/`).

The report includes:
- Experiment metadata table
- Behavioral accuracy and RT distributions
- Audio duration distributions (Natus vs behavioral)
- LangLoc-responsive electrode plots (split-half reliability, condition bar chart, word-boundary time series)
- High-gamma time series across all channels
- Summary of significant channels (unipolar and bipolar)

Requires the MATLAB Report Generator toolbox.
