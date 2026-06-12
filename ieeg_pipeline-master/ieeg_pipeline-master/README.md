# iEEG Preprocessing Pipeline

A MATLAB pipeline for preprocessing intracranial EEG (ECoG/SEEG) data, built around the `ecog_data` class. The pipeline handles signal filtering, artifact rejection, referencing, high-gamma extraction, normalization, and trial epoch extraction.

Developed in the [EvLab](https://evlab.mit.edu) at MIT. Data has been collected across multiple recording institutions; institution-specific scripts and details live in `recording_institutions/`.

---

## Contact

| Person | Role | Email |
|---|---|---|
| **Kumar Duraivel** | Primary maintainer (2023–present) | dsuseendar@gmail.com |
| **Colton Casto** | Original pipeline author | ccasto@mit.edu |

---

## Repository Structure

```
ieeg_pipeline/
├── @ecog_data/                  Core preprocessing class (all methods)
├── recording_institutions/
│   ├── pre-mgh/                 Albany Medical, Mayo Jacksonville, Barnes-Jewish, SLCH
│   │   ├── README.md            ← detailed pre-MGH workflow
│   │   ├── crunch.m             Main entry point (BCI2000 .dat files)
│   │   ├── crunch.sh            SLURM batch wrapper
│   │   ├── info/                Subject operation info files (*_op_info.mat)
│   │   └── demo/                Example preprocessing script
│   └── mgh/                     Massachusetts General Hospital (Natus EDF)
│       ├── README.md            ← detailed MGH workflow
│       ├── langloc_audio_template.m
│       ├── langloc_visual_template.m
│       ├── spatialwm_template.m
│       ├── sentences_template.m
│       ├── speechlangloc_template.m
│       └── demo/                LangLoc Audio full pipeline + report generation
├── expts/                       Experiment stimulus files (shared across sites)
├── utils/                       Statistical utilities, colormaps, IED detection,
│                                Albany MEX I/O, analysis helpers
├── create_info.m                Utility to create op_info files (pre-MGH)
└── README.md
```

---

## Dependencies

- MATLAB R2020b or later
- Signal Processing Toolbox
- Parallel Computing Toolbox (optional, speeds up permutation tests)
- MATLAB Report Generator (optional, required for PDF report generation)

No SPM, FieldTrip, or MNE required.

---

## Quick Start

### Pre-MGH sites (Albany Medical / Mayo Jacksonville / Barnes-Jewish)

Data is stored in BCI2000 `.dat` format.

```matlab
addpath(genpath('/path/to/ieeg_pipeline'));

crunch('AMC082', 'MITLangloc', ...
    'fromScratch', true, ...
    'order',       'defaultSEEGorBOTH', ...
    'doneVisualInspection', true);
```

Edit the path variables at the top of `recording_institutions/pre-mgh/crunch.m` before running.

**→ See [`recording_institutions/pre-mgh/README.md`](recording_institutions/pre-mgh/README.md) for the full workflow, all `crunch.m` parameters, op info format, preprocessing orders, and output structure.**

---

### MGH (Massachusetts General Hospital)

Data is stored in Natus EDF format. Each task has a template in `recording_institutions/mgh/`.

1. Copy the relevant template, edit `SUBJECT`, `DATAPATH`, and `MGH_PREPROC_REPO` at the top.
2. Run section by section, reviewing the trigger alignment plots.
3. For a full walkthrough with report generation: `recording_institutions/mgh/demo/demo_langloc_report.m`.

A companion repository contains subject-specific scripts and audio alignment files:
> **[MGH_IEEG_preproc](https://github.com/YOURUSERNAME/MGH_IEEG_preproc)**

**→ See [`recording_institutions/mgh/README.md`](recording_institutions/mgh/README.md) for the full workflow, data structure, trigger parsing, and template reference.**

---

## The `ecog_data` Class

Constructing the object and preprocessing are separate steps. Once an object is saved, preprocessing can be re-run with a different order without returning to raw data files.

### Constructor

```matlab
obj = ecog_data(for_preproc, subject, experiment, save_filename, save_path, ...
                d_files, file_path, elec_ch_label, elec_ch, ...
                elec_ch_prelim_deselect, elec_ch_type);
```

Key `for_preproc` fields: `elec_data_raw`, `stitch_index_raw`, `sample_freq_raw`, `decimation_freq`, `decimation_factor`, `elecs_per_amp`. See [`recording_institutions/pre-mgh/README.md`](recording_institutions/pre-mgh/README.md) for full field descriptions.

### Preprocessing orders

| Order | Use case |
|---|---|
| `defaultECOG` | ECoG-only (CAR referencing) |
| `defaultSEEGorBOTH` | SEEG or mixed ECoG+SEEG (global mean + bipolar) |
| `defaultMCJandBJH` | MCJ / BJH sites |
| `defaultSEEGorBOTHBroadBand` | MGH SEEG, broadband |
| `preEnvelopeExtractionECOG` | ECoG before envelope extraction |
| `preEnvelopeExtractionSEEGorBOTH` | SEEG/mixed before envelope extraction |
| `preEnvelopeExtractionMCJandBJH` | MCJ/BJH before envelope extraction |

### Post-preprocessing

```matlab
obj.extract_high_gamma('doNapLabFilterExtraction', true);
obj.downsample_signal('decimationFreq', 100);
obj.extract_significant_channel();
obj.extract_time_significance();
obj.extract_normalization_metrics();
obj.normalize_signal('normtype', 'z-score');

[epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', [-0.5 2.0]);
```

---

## Utilities

| Location | Contents |
|---|---|
| `utils/JancaCodePapers/` | IED detection (Janca & Podvalny methods) |
| `utils/albany_mex_files/` | BCI2000 `.dat` file I/O (MEX binaries, multiple platforms) |
| `utils/Colormaps/` | Perceptually uniform colormaps |
| `utils/fdr_bh/` | Benjamini-Hochberg FDR correction |
| `utils/kumar_ieeg_utils/` | Analysis helpers: permutation tests, HDF5 I/O, report generation |
