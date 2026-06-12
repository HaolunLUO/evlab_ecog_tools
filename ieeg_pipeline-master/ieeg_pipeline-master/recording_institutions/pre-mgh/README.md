# Pre-MGH iEEG Preprocessing Pipeline

Preprocessing pipeline for ECoG/SEEG data recorded at **Albany Medical College (AMC)**, **Mayo Clinic Jacksonville (MCJ)**, **Barnes-Jewish Hospital (BJH)**, and **St. Louis Children's Hospital (SLCH)**. Data is acquired using BCI2000 and stored as `.dat` files.

This pipeline was constructed in December 2021 – January 2022 by **Colton Casto** (ccasto@mit.edu), modeled closely on earlier separate pipelines for Albany and MGH data built by Eghbal Hosseini and Hannah Small in 2020–2021.

---

## Contact

| Person | Role | Email |
|---|---|---|
| **Kumar Duraivel** | Primary maintainer (2023–present) | dsuseendar@gmail.com |
| **Colton Casto** | Original pipeline author | ccasto@mit.edu |

---

## Workflow Overview

Preprocessing consists of three steps:

1. **Construct the object** — load raw `.dat` files, subject op info, and build the `ecog_data` object
2. **Preprocess the signal** — apply a preprocessing order (filtering, referencing, high-gamma extraction, normalization)
3. **Extract trial/behavioral information** — parse BCI2000 state variables to build trial timing and events tables

Steps 1 and 3 require access to raw data files. Step 2 can be re-run at any time on a saved object without returning to raw files — useful for trying different preprocessing orders.

---

## Entry Points

### `crunch.m` — Full pipeline in one call

```matlab
% Edit REPO_PATH, MINDHIVE_PATH, and NESE_PATH at the top of crunch.m first

crunch('AMC082', 'MITLangloc', ...
    'fromScratch',          true, ...
    'order',                'defaultSEEGorBOTH', ...
    'doneVisualInspection', true, ...
    'isPlotVisible',        false);
```

| Parameter | Default | Description |
|---|---|---|
| `subject` | required | Subject ID (e.g., `'AMC082'`) |
| `experiment` | required | Experiment name (e.g., `'MITLangloc'`) |
| `fromScratch` | `false` | Build object from raw `.dat` files |
| `order` | `'defaultECOG'` | Preprocessing order to apply |
| `doneVisualInspection` | `false` | Apply previously saved visual inspection markings |
| `preEnvelopeExtraction` | `false` | Save a pre-envelope-extraction version of the signal |
| `addAnatomy` | `false` | Add electrode anatomical labels (VERA/ANT format) |
| `decimation_factor` | `4` | Downsampling factor |
| `deselectElecsPrelim` | `true` | Remove GND/REF/non-cortical channels at construction |
| `isPlotVisible` | `true` | Show figures during preprocessing |

### `crunch.sh` — SLURM batch wrapper

Runs `crunch.m` as a job array, reading subjects and orders from text files:

```bash
sbatch --array=1-$(wc -l < subjects_working.txt) crunch.sh \
    subjects_working.txt orders_working.txt MITLangloc
```

The working text files (`subjects_working.txt`, `orders_working.txt`, etc.) list one entry per line corresponding to each array index.

---

## Subject Operation Info (`info/`)

Each subject requires an `*_op_info.mat` file in `info/`. The file contains a struct named `SUBJECTID_op_info` (e.g., `AMC082_op_info`) with the following fields:

| Field | Description |
|---|---|
| `channel_label` | Cell array of channel names from the amplifier |
| `channel_type` | Cell array of channel types (`'ecog_grid'`, `'ecog_strip'`, `'seeg'`, `'ground'`, etc.) |
| `GND` | Ground channel indices |
| `REF` | Reference channel indices |
| `bad_channels` | Clinician-marked bad channels |
| `skull_eeg_channels` | Skull EEG channel indices |
| `microphone_channels` | Microphone channel indices |
| `EMG_channels` | EMG channel indices |
| `visual_trigger` | Visual trigger channel index or name |
| `button_trigger` | Button trigger channel index or name |
| `buzzer_trigger` | Buzzer trigger channel index or name |
| `audio_trigger` | Audio trigger channel index or name |

Use `create_info.m` at the repository root to create new op info files.

---

## Constructing the `ecog_data` Object

The constructor takes raw signal and metadata and initializes the object for preprocessing:

```matlab
obj = ecog_data(for_preproc, subject, experiment, save_filename, save_path, ...
                d_files, file_path, elec_ch_label, elec_ch, ...
                elec_ch_prelim_deselect, elec_ch_type);
```

### `for_preproc` struct fields

| Field | Description |
|---|---|
| `elec_data_raw` | Raw signal combined across all data files (nChans × nSamples) |
| `stitch_index_raw` | Sample index in raw fs where each data file starts (nFiles × 1) |
| `stitch_index_dec` | Sample index in decimated fs where each data file starts (nFiles × 1) |
| `sample_freq_raw` | Original sampling frequency (Hz) |
| `decimation_freq` | Target downsampled frequency (Hz) |
| `decimation_factor` | `sample_freq_raw / decimation_freq` |
| `elecs_per_amp` | Electrodes per amplifier — used for Common Average Referencing (CAR) |

> **Note:** Even if you do not plan to downsample, `stitch_index_dec`, `decimation_freq`, and `decimation_factor` must still be provided. The pipeline stores both raw and decimated trial timing so that you can move between sampling rates without rebuilding from scratch. If you later want a different decimation frequency you must rebuild the object.

---

## Preprocessing the Signal

### From raw signal

```matlab
obj.preprocess_signal('order', 'defaultSEEGorBOTH');
```

`preprocess_signal()` accepts only named orders — it does not accept arbitrary sequences. To add a new order, add it directly to the method. This was intentional: it prevents accidental step omissions and ensures all orders are reviewed before use.

To modify step-level parameters (e.g., filter cutoffs, outlier thresholds), edit `define_parameters()` in the class.

### From an existing saved object

Load the `.mat` file and call any preprocessing method directly on the object. This overwrites `obj.elec_data` (and `obj.bip_elec_data` if applicable):

```matlab
load(save_filename); % loads obj
obj.zscore_signal();
```

### Preprocessing orders

#### `defaultECOG`
For subjects with only ECoG electrodes. Uses Common Average Referencing (CAR) — preferred over global mean removal when electrode groupings per amplifier are known.

1. Highpass filter
2. Notch filter
3. IED removal
4. Visual inspection
5. **Common average referencing (CAR)**
6. High gamma extraction (gaussian filters)
7. Z-score
8. Downsample
9. Remove outliers

#### `defaultSEEGorBOTH`
For SEEG-only or mixed ECoG+SEEG subjects. Uses global mean removal (preferred when amplifier groupings are uncertain) and adds bipolar referencing for SEEG electrodes.

1. Highpass filter
2. Notch filter
3. IED removal
4. Visual inspection
5. **Global mean removal**
6. **Bipolar referencing** *(SEEG only; stored in `obj.bip_elec_data`)*
7. High gamma extraction
8. Z-score
9. Downsample
10. Remove outliers

#### `defaultMCJandBJH`
For MCJ and BJH subjects. Global mean removal is applied *before* IED removal because global "blips" unique to these sites cause excessive channel rejection during IED removal without prior mean removal.

1. Highpass filter
2. Notch filter
3. **Global mean removal** *(before IED removal — MCJ/BJH specific)*
4. IED removal
5. Visual inspection
6. Bipolar referencing
7. High gamma extraction
8. Z-score
9. Downsample
10. Remove outliers

#### Pre-envelope extraction variants
Produce a broadband version of the signal before high-gamma extraction. Useful when you want to run envelope extraction separately or with different parameters later.

- `preEnvelopeExtractionECOG` — ECoG version (stops after CAR + downsample)
- `preEnvelopeExtractionSEEGorBOTH` — SEEG/mixed version (stops after global mean + bipolar + downsample)
- `preEnvelopeExtractionMCJandBJH` — MCJ/BJH version

---

## Main Signal Preprocessing Methods

These are called internally by `preprocess_signal()` but can also be called directly on a loaded object.

| Method | Description |
|---|---|
| `highpass_filter()` | Highpass filters signal using MATLAB `filtfilt()` |
| `notch_filter()` | Removes 60 Hz line noise and harmonics |
| `remove_IED()` | Marks electrodes with significant IEDs using the Janca protocol. Channels with IEDs are excluded from subsequent steps |
| `visual_inspection()` | Interactive plot for manual channel rejection. Saved markings can be reloaded via `visual_inspection_working.csv` |
| `reference_signal()` | Applies referencing. Args: `'doGlobalMeanRemoval'`, `'doCAR'`, `'doShankCSR'`, `'doBipolarReferencing'` |
| `extract_high_gamma()` | Extracts high-gamma envelope. Args: `'doGaussianFilterExtraction'`, `'doBandpassExtraction'`, `'doNapLabFilterExtraction'` |
| `zscore_signal()` | Z-scores signal |
| `downsample_signal()` | Downsamples signal |
| `remove_outliers()` | Removes outlier timepoints |

All parameters are set in `define_parameters()`. Preprocessing is applied **per file** for experiments with multiple runs.

---

## Post-preprocessing Steps

These are applied after signal preprocessing and can be called on a loaded object:

```matlab
obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 100);
obj.extract_significant_channel();        % windowed permutation test (channel level)
obj.extract_time_significance();          % cluster-corrected time-series test
obj.extract_normalization_metrics();      % compute baseline statistics
obj.normalize_signal('normtype', 'z-score');
```

---

## Trial Epoch Extraction

```matlab
[epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', [-0.5 2.0]);
```

Returns `nChans × nTrials × nSamples`. `epochData_bip` is empty if bipolar referencing was not applied.

---

## Object Output Structure

| Property | Description |
|---|---|
| `elec_data` | Current unipolar signal (nChans × nSamples) |
| `bip_elec_data` | Current bipolar signal (nBipChans × nSamples) |
| `sample_freq` | Current sampling frequency |
| `stitch_index` | File start samples at current sampling rate |
| `elec_ch_label` | Channel labels |
| `elec_ch_clean` | Indices of clean (non-excluded) channels |
| `elec_ch_valid` | Indices of channels used in analysis |
| `elec_ch_with_IED` | Channels excluded by IED removal |
| `elec_ch_with_noise` | Channels excluded by visual inspection |
| `bip_ch_label` | Bipolar channel labels |
| `trial_timing` | Per-trial timing table at current sampling rate |
| `events_table` | Trial metadata (condition, stimulus, probe, RT, accuracy) |
| `condition` | Condition label per trial |
| `session` | Run/session number per trial |
| `for_preproc` | Struct with raw signal, raw timing, and all preprocessing parameters |

---

## Supported Experiments

| Experiment | Description |
|---|---|
| `MITLangloc` | Language localizer (visual presentation, BCI2000) |
| `MITNaturalisticStoriesTask` | Naturalistic stories (continuous audio) |
| `MITConstituentBounds` | Constituent boundary task |
| `MITNLengthSentences` | N-length sentence task |

Stimulus timing files for shared experiments are in `expts/` at the repository root.
