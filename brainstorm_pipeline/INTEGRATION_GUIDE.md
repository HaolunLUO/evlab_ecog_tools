# Brainstorm iEEG Pipeline — Integration Guide

## Overview

This folder integrates a Brainstorm-to-MIT iEEG/SEEG analysis pipeline into
`evlab_ecog_tools`.  The SEEG classes are **subclasses of the canonical
`MGH_utils/@ecog_data_v2`** class, so they stay aligned with `ecog_data_v2` by
inheritance and only override the handful of methods that must differ for SEEG.
Crunched `.mat` files produced here contain an `ecog_data_seeg` object (an
`ecog_data_v2` subclass), so they remain compatible with `ecog_data_v2`-based
analysis from Step 2 onward.

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
│   Variable: obj  (ecog_data_seeg object, an ecog_data_v2 subclass) │
│   Contains: raw signal, trial timing, channel labels, anatomy      │
└──────────────────────────────────┬─────────────────────────────────┘
                                   │
                    ┌──────────────▼───────────────┐
                    │  ecog_data_seeg methods        │
                    │  (preprocess_signal, etc.;     │
                    │   inherits ecog_data_v2)       │
                    └──────────────┬─────────────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │  ecog_sn_data_seeg methods   │
                    │  (test_s_vs_n, lang_resp_    │
                    │   plots, get_summary_stats)  │
                    └──────────────────────────────┘
```

---

## File Reference

| File | Purpose |
|------|---------|
| `ecog_data_seeg.m` | `ecog_data_v2` subclass; inherits everything and overrides only the SEEG-specific methods |
| `ecog_sn_data_seeg.m` | `ecog_data_seeg` subclass: configurable S/N flags, timecourse/barplots, summary stats, cross-task ROI |
| `roi_utils.m` | `save_roi_from_sn_obj`, `apply_roi_to_sn_obj`, `normalize_labels` |
| `plot_anatomy_utils.m` | `plot_sig_electrodes_anatomy`, `get_electrode_coords`, `plot_template_cortex` |
| `complete_mit_pipeline_brainstorm.m` | End-to-end pipeline script (main entry point) |

---

## MATLAB Path Setup

> The SEEG classes are uniquely named (`ecog_data_seeg` / `ecog_sn_data_seeg`),
> so path **order no longer matters** — there is no name collision with the
> legacy `ecog_data` / `ecog_sn_data` classes elsewhere in the repo. Just make
> sure the whole repo is on the path.

```matlab
% Path setup for the SEEG (Brainstorm) pipeline
addpath(genpath('/path/to/evlab_ecog_tools'));   % includes brainstorm_pipeline + MGH_utils

% brainstorm_to_mit_crunched_new.m must also be on path (provided separately)
```

The **BCI2000 pipeline** (`crunch_subject_ALBANY.m`) is completely unaffected:
it uses the legacy `ecog_data` / `ecog_sn_data` classes, which no longer share a
name with anything in `brainstorm_pipeline/`.

---

## Class Compatibility

### `ecog_data_seeg < ecog_data_v2`

`ecog_data_seeg` **is** an `ecog_data_v2` (subclass), so:
- It has the **identical constructor signature** (11 arguments: `for_preproc`,
  `subject`, etc.) and **identical properties** (`elec_data`, `bip_elec_data`,
  `anatomy`, `stats`, etc.) — all inherited.
- Every `ecog_data_v2` method is available unless explicitly overridden below.
- Crunched files written here load as `ecog_data_seeg` objects and work anywhere
  an `ecog_data_v2` object is expected (`isa(obj,'ecog_data_v2')` is true).

**SEEG-specific overrides (the only methods that differ from `ecog_data_v2`):**

| Override | Why it differs for SEEG |
|----------|--------------------------|
| `define_parameters()` | 50 Hz line-noise standard (notch 50/100/150/200/250 Hz, peak 45/50/55 Hz) and excludes Gaussian high-gamma bands that land on 50 Hz harmonics |
| `notch_filter()` | Reports the actual line-noise frequency (reads `for_preproc.filter_params.line_noise_hz`) |
| `extract_shanks()` | Operates only on clean channels, so excluded contacts with unparseable labels can't break shank parsing |
| `reference_signal()` | Derives bipolar pairs directly from channel labels, keeping it consistent with the clean-only `extract_shanks`; also adds optional along-shank **Laplacian** referencing (`doLaplacianReferencing`) ported from `ieeg_pipeline` |
| `preprocess_signal()` | SEEG preprocessing orders that do **NOT** apply CAR before bipolar referencing; adds Laplacian-based orders (`'defaultSEEGLaplacian'`, `'preEnvelopeExtractionSEEGLaplacian'`) |

Everything else used by the pipeline — `extract_high_gamma`, `make_trials`,
`extract_trial_epochs`, `extract_normalization_metrics`, `normalize_signal`,
`extract_time_significance`, `extract_significant_channel`, `combine_data_files`,
`define_clean_channels`, `get_cond_id`, the `stats` property, etc. — is
**inherited unchanged** from `ecog_data_v2`.

---

## Advanced preprocessing (ported from `ieeg_pipeline`)

The repo also ships a newer, more advanced preprocessing pipeline
(`ieeg_pipeline-master/@ecog_data`, Duraivel/Casto, EvLab). Its most useful
SEEG preprocessing capabilities have been brought into `ecog_data_seeg` as
**additive, opt-in** features. The existing default flow (`'defaultSEEG'`) is
**unchanged**; nothing happens unless you explicitly enable these.

| Capability | How to use | Notes |
|------------|-----------|-------|
| **Along-shank Laplacian referencing** | `preprocess_signal('order','defaultSEEGLaplacian')` or `reference_signal('doLaplacianReferencing',true)` | Local spatial reference: endpoints subtract the single adjacent contact, interior contacts subtract the mean of both neighbours. Uses clean channels only (consistent with the SEEG `extract_shanks`). Writes to `elec_data` (no bipolar produced). New orders: `'defaultSEEGLaplacian'`, `'preEnvelopeExtractionSEEGLaplacian'`. |
| **Sharp-artifact detection** | `obj.detect_sharp_artifacts()` (after z-scoring) | Flags sharp transients on the z-scored high-gamma envelope using OR logic over amplitude (`'min_amplitude'`, default 15) and slope (`'min_slope'`, default 10) criteria. Results in `obj.stats.artifact_stats_unipolar` / `.artifact_stats_bipolar`. For `MITNaturalisticStoriesTask` it restricts analysis to `story_*` epochs. |
| **High-gamma smoothing** | `obj.smooth_high_gamma()` (before `make_trials`) | Gaussian smoothing of the (normalized) high-gamma envelope (default 100 ms window, `'window_s'`). Matches the smoothing the `ieeg_pipeline` performs inside `normalize_signal`; kept opt-in here so the inherited `normalize_signal` behaviour is preserved by default. |
| **Bandpass extraction** | `[uni,bip] = obj.extract_bandpass_signal(lo,hi)` | Stitch-aware segment-wise bandpass filter. **Requires `eegfilt` (EEGLAB)** on the path; errors with a clear message otherwise. |
| **Save processed object** | `obj.saveUpdatedObject()` | Saves to `<crunched_file_path>/<subject>_<experiment>_crunched_HG_ZScore.mat`. |

In `complete_mit_pipeline_brainstorm.m` these are exposed as the opt-in flags
`useLaplacianReferencing`, `smoothHighGamma`, and `detectSharpArtifacts` in the
USER SETTINGS block (all default `false`).

> **Why additive instead of re-basing onto `ieeg_pipeline/@ecog_data`?** That
> class is named `ecog_data`, which already exists twice elsewhere in the repo
> (`./ecog_data.m`, `MGH_utils/ecog_data.m`). Subclassing it would reintroduce
> the exact class-name collision the SEEG classes were designed to avoid (the
> superclass would resolve by path order). Keeping the robust, uniquely-named
> `ecog_data_v2` base and porting the advanced steps as overrides/additions
> avoids that collision and keeps the existing flow working, while still making
> the advanced preprocessing available.

External dependencies for the inherited stats methods (already in this repo):
- `remove_bad_trials`, `extractCommonTrials`, `extendTimeEpoch`, `timePermCluster` — `kumar_ieeg_utils/`
- `fdr_bh` — `fdr_bh/`

### `ecog_sn_data_seeg` vs `MGH_utils/ecog_sn_data.m` (old)

`ecog_sn_data_seeg < ecog_data_seeg` is the S-vs-N analysis layer. The old
`MGH_utils/ecog_sn_data.m` extends the legacy `ecog_data` and is independent.

| Feature | Old (`MGH_utils/`) | New (`brainstorm_pipeline/`) |
|---------|-------------------|------------------------------|
| Base class | legacy `ecog_data` | `ecog_data_seeg` (→ `ecog_data_v2`) |
| Condition flags | Hardcoded `'S'` / `'N'` | Configurable via `S_condition_flag`, `N_condition_flag` |
| Timecourse extraction | Not present | `get_timecourses()` |
| Word averages | Not present | `get_word_averages()` |
| Summary stats | Not present | `get_summary_statistics()` (anatomy optional) |
| Cross-task ROI | Not present | Via `apply_roi_to_sn_obj()` |
| Plot methods | Separate `ecog_sn_analysis.m` | Embedded (`plot_timecourse`, `plot_barplot`) |

### CAR is not applied before bipolar referencing (SEEG)

The SEEG preprocessing orders defined in `ecog_data_seeg.preprocess_signal`
(`'defaultSEEG'`, `'defaultSEEGbyShank'`, `'preEnvelopeExtractionSEEG'`) omit the
CAR step. Bipolar referencing already removes shared/common signal between
adjacent contacts, so a prior common-average step is unnecessary and can distort
the local bipolar estimate. (The parent `ecog_data_v2` orders that do use CAR,
e.g. `'defaultECOG'` / `'defaultSEEGorBOTH'`, are still reachable — any order not
recognized by the subclass is delegated to `ecog_data_v2`.)

### 50 Hz vs 60 Hz line noise

| Setting | `ecog_data_seeg` | `ecog_data_v2` |
|---------|------------------|----------------|
| Notch filter harmonics | 50, 100, 150, 200, 250 Hz (European/Chinese standard) | 60, 120, 180, 240 Hz (US standard) |
| Peak filter (noise QC) | 45, 50, 55 Hz | 55, 60, 65 Hz |

**If your recordings were made in a 60 Hz country (US), edit `define_parameters()` in
`brainstorm_pipeline/ecog_data_seeg.m`:**

```matlab
param.line_noise_hz = 60;                   % single source of truth
param.peak.fcenter  = param.line_noise_hz + [-5, 0, 5];   % [55 60 65]
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
- [ ] Whole repo on the MATLAB path (path order no longer matters)
- [ ] `JancaCodePapers/` on path (for IED removal; already in this repo)
- [ ] Anatomy files present at `anatomyPath` (or set `obj.anatomy = []` to skip)
- [ ] Edit `workingDir`, `params.SubjectName`, and `allDataFiles` in the pipeline script
- [ ] For cross-task ROI: run source task (`useROIfromSource=true`, `taskType=roiSourceTask`) first

---

## What the Existing Pipeline is NOT affected by

- `crunch_subject_ALBANY.m` and `general_crunch_script.m` are unchanged
- `MGH_utils/@ecog_data_v2/` is unchanged (the SEEG classes subclass it; they do
  not modify it)
- `MGH_utils/ecog_sn_data.m` (old) is unchanged
- All filter scripts in `ecog-filters/` are unchanged
- BCI2000 I/O in `mex/` and `albany_mex_files/` is unchanged

---

## Decision Summary: Why uniquely-named subclasses?

Earlier, `brainstorm_pipeline/` shipped its own monolithic `ecog_data.m` and
`ecog_sn_data.m`. Because those reused the names `ecog_data` / `ecog_sn_data`
(which also exist in the repo root and in `MGH_utils/`), MATLAB could only ever
resolve one class per name by **path order**. In practice the brainstorm copies
were silently shadowed, so the pipeline never actually ran on them.

The fix:

1. **Unique names** (`ecog_data_seeg`, `ecog_sn_data_seeg`) eliminate the
   collision, so path order is irrelevant.
2. **Subclass `ecog_data_v2`** so the SEEG pipeline stays aligned with the
   canonical class by inheritance and only the genuinely SEEG-specific methods
   are overridden — no large monolithic copy to drift out of sync.
3. **No edits to `ecog_data_v2`** itself, so the ECoG (grid/strip) pipelines that
   depend on it are unaffected.
