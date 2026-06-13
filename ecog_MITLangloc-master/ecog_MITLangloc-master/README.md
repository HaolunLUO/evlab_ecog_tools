# MITLangloc

## Overview

This repository is used to identify and analyze language-responsive ECoG/SEEG channels collected in one of three similar experimental paradigms. Language-responsive channels are defined as channels that respond significantly more to sentences (+meaning, +form) than a perceptually-matched control (often sequences of nonwords like "florp" or "flickit"; -meaning, -form). See **Experiments** for more on the paradigms and **Identification of Language-Responsive Channels** for more on the selection criteria. 

All data was collected by Peter Brunner and colleagues across four collection sites. See the repository below for more information on how these data were collected and preprocessed: 

https://github.mit.edu/ccasto/ecog_pipeline_merged


## Experiments 

The following experiments can be analyzed using this repository. 

### MITLangloc
 
17 subjects to date (w/ usable signal) read sentences and lists of nonwords. Each trial consisted of **12** words/nonwords presented one at a time, and at the end of the sequence, a memory probe was presented. A fixation cross was presented both before and after the probe, and the presentation speed was either fast (450 ms per word), medium (600 ms per word), or slow (750 ms per word) depending on the subjects' preference based on a speed test. Subjects were shown some combination of speed tests and slow/medium/fast langloc. The number of runs varied by subject and presentation speed, but it was most common for subjects to see 2 runs of medium langloc (72 trials per run, 36 trials per condition per run). 

### MITSWJNTask

6 subjects read sentences, lists of words, Jabberwock sentences, and lists of nonwords. Each trial consisted of **8** words/nonwords presented one at a time, and at the end of the sequence, a memory probe was presented. A fixation cross was presented both before and after the probe, and the presentation speed was either fast (450 ms per word) or slow (700 ms per word) depending on the subjects' preference. Subjects who opted for the fast presentation speed saw a total of 80 trials per condition (320 trials total, 10 runs, 32 trials per run, 8 trials per condition per run), and those who opted for the slow presentation speed saw a total of 60 trials per condition (240 trials total, 10 runs, 24 trials per run, 6 trials per condition per run). 

### MIT_SWJN_audio_visual_task

Forthcoming.


## Identification of Language-Responsive Channels

The following procedure is applied to identify language responsive channels: 

1. Take the average response for each word, electrode, and trial
2. Average the response over all 12 words to get 1 value per trial per electrode (`Cond_R` vector)
3. Create a true condition vector where 1 corresponds to sentences and -1 corresponds to nonwords (`Cond` vector)
4. Take the correlation between `Cond_R` and `Cond` (`Corr_SN`)
5. Permutate the `Cond` vector over 10k iterations and correlate with `Cond_R` (`Corr_rand`)
6. Compate `Corr_SN` with `Corr_rand` using a t-test (0.05 threshold)


## Usage 

`analyze_MITLangloc.m`

### Arguments
  * **Type of Analysis**
    * `'fromScratch'` (default : `false`) :
      * Constructs `ecog_sn_data` object from scratch and identifies electrodes with significant Sentences>Nonwords activity. Also produces S>N plots for **all** electrodes (averaged by words) and a histogram of correlations between the `Cond_R` and permutted `Cond` vectors (null distirbution, 10k iterations, see **Identification of Language-Responsive Channels**) that is saved in `output/all_electrodes`. 
      * This step must be done before any subsequent analyses can be performed. 
      * ***NB: This analysis will overwrite the `ecog_data` object in the crunched directory with a new `ecog_sn_data` object.***
    * `'averageSub'` (default : `false`) :
      * Produces plots for all language-responsive channels and a plot of all language-responsive channels for a given subject averaged together with word-averaged and timecourse signal (saved in `output/all_lang_electrodes` and `output/sub_average`).
    * `'summaryStatistics'` (default : `false`) :
      * Compiles list of significant (unipolar and bipolar) channels for all specified subjects. 
    * `'toStruct'` (default : `false`) :
      * Loads existing `ecog_sn_data` object and saves it as a structure in the `structures/` subdirectory of the crunched directory. 
  * **Subjects**
    * `'doOneSub'` (default : `[]`) :
      * Only runs analysis on specified subject (e.g., `'AMC091'`).
    * `'doOneSub'` (default : `[]`) :
      * Only runs analysis on specified subject **group** (e.g., `'BJH'`).
  * **Other**
    * `'words'` (default : `[1:12]`) :
      * Number of words presented to subjects in a given experiment to be used for the purposes of word-level averaging and plotting. 
    * `'SConditionName'` (default : `'Sentences'`) : 
      * Name of S condition in `ecog_data` objects (e.g., `'Sentences'` for MITLangloc and `'SENTENCES'` for MITSWJNTask).
    * `'NConditionName'` (default : `'Jabberwocky'`) : 
      * Name of N condition in `ecog_data` objects (e.g., `'Jabberwocky'` for MITLangloc (not actually Jabberwocky, just mislabeled in raw data files) and `'NONWORDS'` for MITSWJNTask).
    * `'experiment'` (default : `'MITLangloc'`) : 
      * Experiment name (alternatives : `'MITSWJNTask'`)
    
### Example

`analyze_MITLangloc('averageSub',true,'SConditionName','SENTENCES','NConditionName','NONWORDS','experiment','MITSWJNTask')`

Produces plots for all lang-resp chans and average plots by subject for all MITSWJNTask subjects (AMC026, AMC029, AMC031, AMC037, AMC038, AMC044) from pre-existing `ecog_sn_data` objects.


## TODO 

* Split trials in half for identification of lang-resp chans (half to identify, half to produce S>N plots / measure effect sizes)
* Review S>N selection criteria w/ Ev (and compared to PNAS method)
