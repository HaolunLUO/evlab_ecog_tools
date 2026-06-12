# Brainstorm iEEG Pipeline — Integration Guide

## Overview

This folder integrates a Brainstorm-to-MIT iEEG/SEEG analysis pipeline into
`evlab_ecog_tools`.  The SEEG classes are now **re-based onto the advanced EvLab
`ieeg_pipeline` engine**, which is vendored here as
**`brainstorm_pipeline/@ecog_data_ieeg`** (a uniquely-named copy of
`ieeg_pipeline-master/@ecog_data`). `ecog_data_seeg` subclasses that engine,
inheriting all of its advanced preprocessing methods and overriding only the
handful that must differ for SEEG. Crunched `.mat` files produced here contain
an `ecog_data_seeg` object (an `ecog_data_ieeg` subclass).

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
│   Variable: obj  (ecog_data_seeg object, an ecog_data_ieeg subclass)│
│   Contains: raw signal, trial timing, channel labels, anatomy      │
└──────────────────────────────────┬─────────────────────────────────┘
                                   │
                    ┌──────────────▼───────────────┐
                    │  ecog_data_seeg methods        │
                    │  (preprocess_signal, etc.;     │
                    │   inherits @ecog_data_ieeg)    │
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
| `@ecog_data_ieeg/` | Advanced iEEG/SEEG preprocessing engine, vendored from `ieeg_pipeline-master/@ecog_data` (renamed to avoid the `ecog_data` name collision). Kept verbatim apart from the class name and `arguments obj` type validators, for easy re-sync. |
| `ecog_data_seeg.m` | `ecog_data_ieeg` subclass; inherits the engine and overrides only the SEEG-specific methods |
| `ecog_sn_data_seeg.m` | `ecog_data_seeg` subclass: configurable S/N flags, timecourse/barplots, summary stats, cross-task ROI |
| `roi_utils.m` | `save_roi_from_sn_obj`, `apply_roi_to_sn_obj`, `normalize_labels` |
| `plot_anatomy_utils.m` | `plot_sig_electrodes_anatomy`, `get_electrode_coords`, `plot_template_cortex` |
| `complete_mit_pipeline_brainstorm.m` | End-to-end pipeline script (main entry point) |

---

## MATLAB Path Setup

> The SEEG classes are uniquely named (`ecog_data_seeg` / `ecog_sn_data_seeg`)
> and the engine they subclass is the uniquely-named, vendored
> `@ecog_data_ieeg`, so path **order no longer matters** — there is no name
> collision with the legacy `ecog_data` / `ecog_sn_data` classes elsewhere in
> the repo. Just make sure the whole repo is on the path.

```matlab
% Path setup for the SEEG (Brainstorm) pipeline
addpath(genpath('/path/to/evlab_ecog_tools'));   % includes brainstorm_pipeline + MGH_utils

% brainstorm_to_mit_crunched_new.m must also be on path (provided separately)
```

`@ecog_data_ieeg` is self-contained for the core preprocessing, but a few
inherited engine methods call shared utilities that live elsewhere in the repo
(all picked up by `genpath`): `remove_bad_trials`, `extractCommonTrials`,
`extendTimeEpoch`, `timePermCluster` (`kumar_ieeg_utils/`), `fdr_bh` (`fdr_bh/`),
the HDF5 helpers used by `output_xarray*` (`ieeg_pipeline-master/.../utils/`),
and `eegfilt` (EEGLAB; only for `extract_bandpass_signal`).

The **BCI2000 pipeline** (`crunch_subject_ALBANY.m`) is completely unaffected:
it uses the legacy `ecog_data` / `ecog_sn_data` classes, which no longer share a
name with anything in `brainstorm_pipeline/`.

---

## Class Compatibility

### `ecog_data_seeg < ecog_data_ieeg`

`ecog_data_seeg` **is** an `ecog_data_ieeg` (subclass of the advanced engine), so:
- It has the **identical constructor signature** (11 arguments: `for_preproc`,
  `subject`, etc.) and **identical properties** (`elec_data`, `bip_elec_data`,
  `anatomy`, `stats`, etc.) — all inherited from the engine.
- Every `ecog_data_ieeg` (engine) method is available unless explicitly
  overridden below.
- `isa(obj,'ecog_data_ieeg')` is true.

**Inherited from the advanced engine (`@ecog_data_ieeg`), unchanged:**
`extract_high_gamma`, `normalize_signal` (now Gaussian-smooths the HG envelope),
`downsample_signal`, `make_trials`, `highpass_filter`,
`remove_IED`, `visual_inspection`, `extract_significant_channel`,
`extract_time_significance`, `extract_normalization_metrics` (baseline anchored
to `probe_key = 1` via the `extract_trial_epochs` override), `combine_data_files`,
`define_clean_channels`, `first_step`, `output_data_structures`, `output_xarray`,
`output_xarray_minimal`, `plus` (recording concatenation), `zscore_signal`, the
`stats` property, etc.

**SEEG-specific overrides (the only methods that differ from the engine):**

| Override | Why it differs for SEEG |
|----------|--------------------------|
| `define_parameters()` | 50 Hz line-noise standard (notch 50/100/150/200/250 Hz, peak 45/50/55 Hz) and excludes Gaussian high-gamma bands that land on 50 Hz harmonics (the engine defaults to the 60 Hz US standard) |
| `measure_line_noise()` | Measures noise power at the configured line-noise frequency (50 Hz, via the peak filters from `define_parameters`) and reports it, instead of the engine's hardcoded `"Measuring 60Hz noise power"` message |
| `notch_filter()` | Reports the actual line-noise frequency (reads `for_preproc.filter_params.line_noise_hz`) and runs the interactive noisy-channel review (the engine notch is clean-channel-only and non-interactive) |
| `plot_line_noise()` | 50 Hz axis labels; figure name from `crunched_file_name` |
| `extract_shanks()` | Operates only on clean channels, so excluded contacts with unparseable labels can't break shank parsing |
| `reference_signal()` | Derives bipolar pairs directly from channel labels, keeping it consistent with the clean-only `extract_shanks`; also adds along-shank **Laplacian** referencing (`doLaplacianReferencing`) |
| `preprocess_signal()` | SEEG preprocessing orders that do **NOT** apply CAR before bipolar referencing; adds Laplacian-based orders (`'defaultSEEGLaplacian'`, `'preEnvelopeExtractionSEEGLaplacian'`); unrecognized orders are delegated to the engine |
| `extract_trial_epochs()` | Supports **both** a string `key` (`'key','word_1'`) and the engine's numeric `probe_key` (default 1, used when `key` is empty), with the window rounded to whole samples. This lets the brainstorm scripts and the inherited engine methods (`extract_significant_channel`, `extract_time_significance`, `extract_normalization_metrics`) both anchor epochs. |
| `get_cond_id` / `get_cond_resp` / `get_value` / `get_ave_cond_trial` | Kept for the string condition flags and the exact table layout the S-vs-N analysis layer consumes |

**SEEG-tuned helpers also on the subclass** (thin variants of, or extras
alongside, the engine equivalents): `detect_sharp_artifacts` (inputParser-based;
avoids `parfor`), `smooth_high_gamma` (explicit HG smoothing), and
`extract_bandpass_signal` (errors clearly if `eegfilt`/EEGLAB is missing).

---

## What the re-base onto the engine changes

Because the base class is now the advanced engine rather than `ecog_data_v2`,
the inherited methods are the **advanced** versions. The most visible behavioural
differences vs. the previous `ecog_data_v2`-based pipeline:

| Area | Now (engine) | Before (`ecog_data_v2`) |
|------|--------------|--------------------------|
| `normalize_signal` | Z-scores **and Gaussian-smooths** the HG envelope (100 ms window) | Z-score only |
| `extract_normalization_metrics` | Engine version (inherited): baseline `[-0.5 0]` anchored to `probe_key = 1` (first word onset), over all/sampled trials | `key='fix'`/`'word_1'`, padded window, common good trials |
| New engine methods | `output_data_structures`, `output_xarray`, `output_xarray_minimal` (Python/HDF5 export), `plus` (concatenate recordings), `detect_sharp_artifacts`, `extract_bandpass_signal`, `extract_significant_channel`, `extract_time_significance` | not present / not used |
| Referencing options | adds along-shank **Laplacian** plus broadband orders (`'defaultSEEGorBOTHBroadBand'`, `'defaultSEEGLaplacian'`) | bipolar / CAR / CSR only |

### Default pipeline (`complete_mit_pipeline_brainstorm.m`)

The script runs the ieeg_pipeline feature chain that feeds the language channel
selection:

```matlab
% Preprocess (STEP 4, on the preproc object):
obj.preprocess_signal('order','defaultSEEGorBOTHBroadBand');  % highpass,notch,IED,CAR,Laplacian,bipolar (broadband)
obj.extract_high_gamma('doNapLabFilterExtraction', true);     % NAPLAB high-gamma envelope
obj.downsample_signal('decimationFreq', 200);

% Significance + baseline z-score (STEP 7.5, on the analysis object):
sn_obj.extract_significant_channel();    % -> sn_obj.stats.sig_hg_channel
sn_obj.extract_time_significance();      % -> sn_obj.stats.time_series
sn_obj.extract_normalization_metrics();  % -> sn_obj.stats.normMetrics (baseline = probe_key 1)
sn_obj.normalize_signal('normtype','z-score');   % z-score (+ engine smoothing)
sn_obj.make_trials();

% -> language channel selection
sn_obj.test_s_vs_n(...);
```

Configured via `USER SETTINGS`: `preprocOrder` (default `'defaultSEEGorBOTHBroadBand'`),
`decimationFreq` (default `200`), and the optional `detectSharpArtifacts` flag.
HG smoothing is not a separate flag because the inherited `normalize_signal`
already performs it. Note `extract_significant_channel` / `extract_time_significance`
run permutation tests and can be slow on high channel counts.

> **Why a vendored engine (`@ecog_data_ieeg`) rather than subclassing
> `ieeg_pipeline-master/@ecog_data` directly?** The upstream class is named
> `ecog_data`, which already exists twice elsewhere in the repo (`./ecog_data.m`,
> `MGH_utils/ecog_data.m`). Subclassing the bare name would resolve the
> superclass by MATLAB path order and could silently pick a legacy class. The
> uniquely-named vendored copy removes that collision while keeping the engine
> byte-for-byte re-syncable from upstream (only the class name and the
> `arguments obj` type validators were renamed).

### `ecog_sn_data_seeg` vs `MGH_utils/ecog_sn_data.m` (old)

`ecog_sn_data_seeg < ecog_data_seeg` is the S-vs-N analysis layer. The old
`MGH_utils/ecog_sn_data.m` extends the legacy `ecog_data` and is independent.

| Feature | Old (`MGH_utils/`) | New (`brainstorm_pipeline/`) |
|---------|-------------------|------------------------------|
| Base class | legacy `ecog_data` | `ecog_data_seeg` (→ `ecog_data_ieeg`) |
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
the local bipolar estimate. (The engine orders that do use CAR, e.g.
`'defaultECOG'` / `'defaultSEEGorBOTH'`, are still reachable — any order not
recognized by the subclass is delegated to `@ecog_data_ieeg`.)

### 50 Hz vs 60 Hz line noise

| Setting | `ecog_data_seeg` | `@ecog_data_ieeg` (engine) |
|---------|------------------|----------------------------|
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
- `MGH_utils/@ecog_data_v2/` is unchanged (no longer used by the SEEG pipeline,
  but still present for any other consumers)
- `ieeg_pipeline-master/@ecog_data` (upstream) is unchanged — `@ecog_data_ieeg`
  is a renamed vendored copy of it
- `MGH_utils/ecog_sn_data.m` (old) is unchanged
- All filter scripts in `ecog-filters/` are unchanged
- BCI2000 I/O in `mex/` and `albany_mex_files/` is unchanged

---

## Decision Summary: re-basing onto the advanced engine

Earlier, `brainstorm_pipeline/` shipped its own monolithic `ecog_data.m` /
`ecog_sn_data.m`, which collided by name with the repo-root and `MGH_utils/`
classes and were silently shadowed. That was fixed by giving the SEEG classes
**unique names** (`ecog_data_seeg`, `ecog_sn_data_seeg`) and subclassing
`ecog_data_v2`.

This update re-bases the SEEG classes onto the newer, more advanced EvLab
`ieeg_pipeline` engine:

1. **Vendor the engine under a unique name** (`@ecog_data_ieeg`, copied from
   `ieeg_pipeline-master/@ecog_data`). The upstream class is named `ecog_data`
   and would collide with the two legacy `ecog_data` classes; the unique name
   keeps path order irrelevant. Only the class name and the `arguments obj` type
   validators were renamed, so it stays re-syncable from upstream.
2. **`ecog_data_seeg < ecog_data_ieeg`** inherits every advanced engine method
   and overrides only the genuinely SEEG-specific ones.
3. **No edits to the upstream `ieeg_pipeline-master/@ecog_data`** or to
   `MGH_utils/@ecog_data_v2`, so other consumers are unaffected.
