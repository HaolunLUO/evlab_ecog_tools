# Brainstorm iEEG Pipeline — Integration Guide

## Overview

This folder integrates a Brainstorm-to-MIT iEEG/SEEG analysis pipeline into
`evlab_ecog_tools`.  Both pipelines produce identical **crunched `.mat` files**
containing an `ecog_data` object, making them interoperable from Step 2 onward.

---

## How the Two Pipelines Relate

```
┌─────────────────────────────────────────────────────────────────────┐
│                      RAW DATA                                        │
│  BCI2000 .dat (MIT EEG lab)    │  Brainstorm .mat (SEEG/iEEG)        │
└──────────────┬──────────────────┴──────────────┬────────────────────┘
               │                                 │
        load_bcidat.m                 brainstorm_to_mit_crunched_new.m
               │                                 │
               ▼                                 ▼
┌────────────────────────────────────────────────────────────────────┐
│            CRUNCHED .mat  ← COMMON INTEGRATION POINT              │
│   Variable: obj  (ecog_data object)                                │
│   Contains: raw signal, trial timing, channel labels, anatomy      │
└──────────────────────────────────┬─────────────────────────────────┘
                                   │
                         ┌─────────▼──────────┐
                         │  ecog_data methods  │
                         │  (preprocess_signal,│
                         │   make_trials, etc) │
                         └─────────┬───────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │  ecog_sn_data methods        │
                    │  (test_s_vs_n, lang_resp_    │
                    │   plots, get_summary_stats)  │
                    └──────────────────────────────┘
```

---

## File Reference

| File | Purpose |
|------|---------|
| `ecog_data.m` | Monolithic version of `ecog_data_v2`; identical interface, adds SEEG-specific improvements |
| `ecog_sn_data.m` | Full analysis class: configurable S/N flags, timecourse/barplots, summary stats, cross-task ROI |
| `roi_utils.m` | `save_roi_from_sn_obj`, `apply_roi_to_sn_obj`, `normalize_labels` |
| `plot_anatomy_utils.m` | `plot_sig_electrodes_anatomy`, `get_electrode_coords`, `plot_template_cortex` |
| `complete_mit_pipeline_brainstorm.m` | End-to-end pipeline script (main entry point) |

---

## MATLAB Path Setup

> **Critical:** `brainstorm_pipeline/` must be added **before** `MGH_utils/` to
> ensure MATLAB resolves `ecog_data` and `ecog_sn_data` to the new classes.

```matlab
% Recommended path setup for Brainstorm pipeline
addpath(genpath('/path/to/evlab_ecog_tools/brainstorm_pipeline'));  % FIRST
addpath(genpath('/path/to/evlab_ecog_tools'));                       % rest of repo

% brainstorm_to_mit_crunched_new.m must also be on path (provided separately)
```

When running the **BCI2000 pipeline** (crunch_subject_ALBANY.m), do NOT add
`brainstorm_pipeline/` to the path — the old `ecog_data` / `ecog_sn_data` in
`MGH_utils/` will be used instead.

---

## Class Compatibility

### `ecog_data` vs `ecog_data_v2`

The new `brainstorm_pipeline/ecog_data.m` and `MGH_utils/@ecog_data_v2/` share:
- **Identical constructor signature** (11 arguments: `for_preproc`, `subject`, etc.)
- **Identical properties** (`elec_data`, `bip_elec_data`, `anatomy`, etc.)

This means:
- Crunched files written by `ecog_data_v2` can be loaded by `ecog_sn_data`
  (the new Brainstorm-pipeline class)
- Crunched files written by the Brainstorm pipeline can be used with
  `ecog_data_v2` methods

### `ecog_sn_data` (new) vs `MGH_utils/ecog_sn_data.m` (old)

| Feature | Old (`MGH_utils/`) | New (`brainstorm_pipeline/`) |
|---------|-------------------|------------------------------|
| Condition flags | Hardcoded `'S'` / `'N'` | Configurable via `S_condition_flag`, `N_condition_flag` |
| Preprocessing methods | Inherited from old `ecog_data` | Inherited from new `ecog_data` |
| Timecourse extraction | Not present | `get_timecourses()` |
| Word averages | Not present | `get_word_averages()` |
| Summary stats | Not present | `get_summary_statistics()` |
| Cross-task ROI | Not present | Via `apply_roi_to_sn_obj()` |
| Plot methods | Separate `ecog_sn_analysis.m` | Embedded (`plot_timecourse`, `plot_barplot`) |
| Anatomy optional | No | Yes (gracefully degrades without anatomy) |

### Features ported FROM `ecog_data_v2` INTO `brainstorm_pipeline/ecog_data.m`

These methods existed in `ecog_data_v2` but were absent in the original shared
iEEG pipeline code. They have all been added to `brainstorm_pipeline/ecog_data.m`:

| Method | Purpose |
|--------|---------|
| `stats` property | Struct that accumulates all analysis outputs |
| `get_cond_id()` | Returns a logical row vector of trial indices for a condition |
| `extract_trial_epochs()` | 3-D `[nChans × nTrials × nSamples]` epoch extraction (used for normalization and significance) |
| `extract_normalization_metrics()` | Computes per-channel `[mean, std]` from fixation/baseline epochs; stores in `obj.stats.normMetrics` |
| `normalize_signal()` | Normalizes `elec_data` and `bip_elec_data` using those metrics; 6 methods: `z-score`, `mean-sub`, `perc-change`, `ratio`, `log-ratio`, `norm` |
| `extract_time_significance()` | Cluster-permutation test at every time point; results in `obj.stats.time_series.pSigChan` |
| `extract_significant_channel()` | Per-channel permutation test (epoch > baseline power) + FDR correction; results in `obj.stats.sig_hg_channel` |
| `doNapLabFilterExtraction` | Third high-gamma extraction method via NAPLAB Columbia filterbank (`naplab_filterbank` static method); alternative to Chang-lab Gaussian |

External dependencies for the stats methods (already in this repo):
- `remove_bad_trials`, `extractCommonTrials`, `extendTimeEpoch`, `timePermCluster` — `kumar_ieeg_utils/`
- `fdr_bh` — `fdr_bh/`

### Other improvements in `brainstorm_pipeline/ecog_data.m` vs `ecog_data_v2`

1. `extract_shanks()` — operates only on clean channels (avoids label-parse failures on reference/excluded contacts)
2. `combine_data_files()` — explicit per-segment sample offset (fixes multi-block alignment)
3. `reference_signal()` — full inline bipolar referencing (no separate script needed)
4. `get_summary_statistics()` (in `ecog_sn_data`) — anatomy is optional

### 50 Hz vs 60 Hz line noise

| Setting | `brainstorm_pipeline/ecog_data.m` | `ecog_data_v2` |
|---------|----------------------------------|----------------|
| Notch filter harmonics | 50, 100, 150, 200, 250 Hz (European/Chinese standard) | 60, 120, 180, 240 Hz (US standard) |
| Peak filter (noise QC) | 45, 50, 55 Hz | 55, 60, 65 Hz |

**If your recordings were made in a 60 Hz country (US), edit `define_parameters()` in
`brainstorm_pipeline/ecog_data.m`:**

```matlab
param.notch.fcenter = [60, 120, 180, 240];  % US line noise
param.peak.fcenter  = [55, 60, 65];
```

---

## Anatomy Format

The Brainstorm pipeline uses the VERA anatomy format:

```
obj.anatomy.subject_space.tala.electrodes  [nElec × 3]  subject-space XYZ
obj.anatomy.mni_space.tala.electrodes      [nElec × 3]  MNI-space XYZ
obj.anatomy.template_brain.cortex          struct        brain mesh
obj.anatomy.mapping                        {nChan × 1}  chan → anatomy index
obj.anatomy.hemisphere                     {nElec × 1}  'left' | 'right'
```

The **old** pipeline stores coordinates directly on the object:

```
obj.elec_ch_pos_mni   [nChan × 3]   MNI
obj.elec_ch_pos_anat  [nChan × 3]   subject-space
```

These are **not interchangeable** — anatomy plotting code must match the format
of the object being analyzed.

---

## Cross-Task ROI Workflow

```matlab
% 1. Run MITSWJNTask to define ROI
taskType = 'MITSWJNTask';
useROIfromSource = true;
roiSourceTask    = 'MITSWJNTask';
% ... run pipeline ... ROI saved to output/MITSWJNTask/<subject>_ROI_from_MITSWJNTask.mat

% 2. Apply that ROI to a WM task
taskType = 'WM';
useROIfromSource = true;
roiSourceTask    = 'MITSWJNTask';
% ... run pipeline ... ROI loaded and applied; plots generated for WM data
```

Label matching in `apply_roi_to_sn_obj` is case/punctuation insensitive, so
`'LF12'`, `'lf12'`, and `'LF-12'` all match.

---

## Quick-Start Checklist

- [ ] `brainstorm_to_mit_crunched_new.m` on MATLAB path
- [ ] `brainstorm_pipeline/` folder added to MATLAB path **before** `MGH_utils/`
- [ ] `JancaCodePapers/` on path (for IED removal; already in this repo)
- [ ] Anatomy files present at `anatomyPath` (or set `obj.anatomy = []` to skip)
- [ ] Edit `workingDir`, `params.SubjectName`, and `allDataFiles` in the pipeline script
- [ ] For cross-task ROI: run source task (`useROIfromSource=true`, `taskType=roiSourceTask`) first

---

## What the Existing Pipeline is NOT affected by

- `crunch_subject_ALBANY.m` and `general_crunch_script.m` are unchanged
- `MGH_utils/@ecog_data_v2/` is unchanged  
- `MGH_utils/ecog_sn_data.m` (old) is unchanged
- All filter scripts in `ecog-filters/` are unchanged
- BCI2000 I/O in `mex/` and `albany_mex_files/` is unchanged

---

## Decision Summary: Why this folder structure?

The key trade-off was between:

**Option A — Rename the new classes** (`ecog_data_brainstorm` etc.)  
→ Avoids all MATLAB path conflicts but requires changes to `ecog_sn_data.m`'s
  `extends` declaration and the main pipeline script.

**Option B — Isolate in a subfolder** (this implementation)  
→ Zero risk to existing code. MATLAB resolves classes by path order, so adding
  `brainstorm_pipeline/` first gives priority to the new classes for Brainstorm
  runs, while the BCI2000 pipeline uses the old classes as before. The crunched
  `.mat` format is the single integration point.

Option B was chosen because it requires no edits to existing files and the
shared crunched-file format already provides the correct integration boundary.
