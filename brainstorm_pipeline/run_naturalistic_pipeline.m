function summary = run_naturalistic_pipeline(obj, allDataFiles, taskConfig, params, outputDir, opts)
% RUN_NATURALISTIC_PIPELINE  Preprocess and segment long-form naturalistic SEEG.
%
%   Mirrors the workflow in MITNatural.m:
%     1) Extract NA/NA2/NA3 markers from Brainstorm exports
%     2) Broadband preprocess (defaultSEEGorBOTHBroadBand), NAPLAB HG, downsample
%     3) Build a parallel broadband-referenced signal for raw segment export
%     4) Cut HG + broadband segments and write summary files
%
%   summary = run_naturalistic_pipeline(obj, allDataFiles, taskConfig, params, outputDir, opts)
%
%   Required opts fields:
%     .crunchedFile
%
%   Optional opts fields:
%     .forceReprocess         (default false)
%     .isPlotVisible          (default false)
%     .doneVisualInspection   (default true)
%     .preprocOrder           (default 'defaultSEEGorBOTHBroadBand')
%     .decimationFreq         (default 500)
%     .saveBroadbandSegments  (default true)

if nargin < 6 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'crunchedFile') || isempty(opts.crunchedFile)
    error('opts.crunchedFile is required.');
end

crunchedFile = char(opts.crunchedFile);
forceReprocess = get_naturalistic_opt(opts, 'forceReprocess', false);
isPlotVisible = get_naturalistic_opt(opts, 'isPlotVisible', false);
doneVisualInspection = get_naturalistic_opt(opts, 'doneVisualInspection', true);
preprocOrder = char(get_naturalistic_opt(opts, 'preprocOrder', 'defaultSEEGorBOTHBroadBand'));
decimationFreq = get_naturalistic_opt(opts, 'decimationFreq', 500);
saveBroadbandSegments = get_naturalistic_opt(opts, 'saveBroadbandSegments', true);

outputDir = char(outputDir);
segmentDir = fullfile(outputDir, 'segments');
rawSegmentDir = fullfile(outputDir, 'segments_raw');
if ~exist(segmentDir, 'dir'); mkdir(segmentDir); end
if saveBroadbandSegments && ~exist(rawSegmentDir, 'dir')
    mkdir(rawSegmentDir);
end

%% Extract marker events from Brainstorm exports
fprintf('\n=== NATURALISTIC: EXTRACT EVENTS ===\n');
naturalEvents = extract_naturalistic_events(allDataFiles, taskConfig, obj);
save(crunchedFile, 'obj', 'naturalEvents', '-v7.3');

%% Preprocess: broadband referencing, HG extraction, downsample (same as STEP 4)
fprintf('\n=== NATURALISTIC: PREPROCESS ===\n');
needPreproc = forceReprocess || ~isfield(obj.for_preproc, 'order') ...
    || isempty(obj.for_preproc.order);

if needPreproc
    fprintf('Preprocessing order: %s\n', preprocOrder);
    fprintf('This may take a while for long naturalistic recordings...\n');
    obj.preprocess_signal('order', preprocOrder, ...
        'isPlotVisible', isPlotVisible, ...
        'doneVisualInspection', doneVisualInspection);

    fprintf('Extracting high-gamma envelope (NAPLAB filterbank)...\n');
    obj.extract_high_gamma('doNapLabFilterExtraction', true);

    fprintf('Downsampling HG to %d Hz...\n', decimationFreq);
    obj.downsample_signal('decimationFreq', decimationFreq);
else
    fprintf('Preprocessing already present; ensuring HG + downsample...\n');
    obj = ensure_naturalistic_hg_ready(obj, decimationFreq);
end

fprintf('Unipolar: [%d x %d] | Bipolar: [%d x %d] | Fs: %.1f Hz | Duration: %.1f s\n', ...
    size(obj.elec_data, 1), size(obj.elec_data, 2), ...
    size(obj.bip_elec_data, 1), size(obj.bip_elec_data, 2), ...
    obj.sample_freq, size(obj.elec_data, 2) / obj.sample_freq);

save(crunchedFile, 'obj', 'naturalEvents', '-v7.3');

%% Optional broadband signal for raw segment export
bb_elec_data = [];
bb_bip_elec_data = [];
bb_sample_freq = obj.sample_freq;
hg_sample_freq = obj.sample_freq;
if saveBroadbandSegments
    fprintf('\n=== NATURALISTIC: BROADBAND SIGNAL ===\n');

    hg_elec_data = obj.elec_data;
    hg_bip_elec_data = obj.bip_elec_data;
    hg_sample_freq = obj.sample_freq;
    hg_stitch_index = obj.stitch_index;

    saved_IED_chans = obj.elec_ch_with_IED;
    saved_noise_chans = obj.elec_ch_with_noise;
    saved_user_chans = obj.elec_ch_user_deselect;
    saved_prelim_chans = obj.elec_ch_prelim_deselect;

    obj.first_step('doneVisualInspection', doneVisualInspection);
    obj.elec_ch_with_IED = saved_IED_chans;
    obj.elec_ch_with_noise = saved_noise_chans;
    obj.elec_ch_user_deselect = saved_user_chans;
    obj.elec_ch_prelim_deselect = saved_prelim_chans;
    obj.define_clean_channels();
    obj.for_preproc.isPlotVisible = false;

    obj.highpass_filter();
    obj.notch_filter();
    obj.elec_ch_with_noise = saved_noise_chans;
    obj.define_clean_channels();
    obj.reference_signal('doCAR', true, 'doBipolarReferencing', true);
    obj.downsample_signal('decimationFreq', decimationFreq);

    bb_elec_data = obj.elec_data;
    bb_bip_elec_data = obj.bip_elec_data;
    bb_sample_freq = obj.sample_freq;

    obj.elec_data = hg_elec_data;
    obj.bip_elec_data = hg_bip_elec_data;
    obj.sample_freq = hg_sample_freq;
    obj.stitch_index = hg_stitch_index;
    obj.elec_ch_with_IED = saved_IED_chans;
    obj.elec_ch_with_noise = saved_noise_chans;
    obj.elec_ch_user_deselect = saved_user_chans;
    obj.elec_ch_prelim_deselect = saved_prelim_chans;
    obj.define_clean_channels();
end

%% Adjust event samples if downsampling occurred
fprintf('\n=== NATURALISTIC: SEGMENT MARKERS ===\n');
naturalEvents = adjust_natural_events_for_decimation(obj, naturalEvents);

allMarkers = taskConfig.markers;
markerInfo = struct('name', {}, 'samples', {}, 'times', {}, 'duration', {});

fprintf('Searching for markers: %s\n', strjoin(allMarkers, ', '));
fprintf('All events in data: %s\n', strjoin(unique(naturalEvents.labels), ', '));

for m = 1:numel(allMarkers)
    markerName = allMarkers{m};
    markerIdx = find(strcmpi(naturalEvents.labels, markerName));

    if isfield(taskConfig.segmentDurations, markerName)
        configuredDuration = taskConfig.segmentDurations.(markerName);
    else
        configuredDuration = 300;
        warning('Duration not specified for marker "%s"; using %.1f s.', ...
            markerName, configuredDuration);
    end

    markerInfo(m).name = markerName;
    markerInfo(m).duration = configuredDuration;
    if ~isempty(markerIdx)
        markerInfo(m).samples = naturalEvents.samples(markerIdx);
        markerInfo(m).times = naturalEvents.times(markerIdx);
        fprintf('  %s: %d instance(s), configured duration %.1f s\n', ...
            markerName, numel(markerIdx), configuredDuration);
    else
        warning('Marker "%s" not found.', markerName);
        markerInfo(m).samples = [];
        markerInfo(m).times = [];
    end
end

%% Collect and sort all marker instances
allMarkerSamples = [];
allMarkerNames = {};
allMarkerIdx = [];
allMarkerDurations = [];

for m = 1:numel(markerInfo)
    if isempty(markerInfo(m).samples)
        continue;
    end
    for j = 1:numel(markerInfo(m).samples)
        allMarkerSamples(end+1) = markerInfo(m).samples(j); %#ok<AGROW>
        allMarkerNames{end+1} = markerInfo(m).name; %#ok<AGROW>
        allMarkerIdx(end+1) = j; %#ok<AGROW>
        allMarkerDurations(end+1) = markerInfo(m).duration; %#ok<AGROW>
    end
end

if isempty(allMarkerSamples)
    error('No markers found. Expected: %s', strjoin(allMarkers, ', '));
end

[allMarkerSamples, sortIdx] = sort(allMarkerSamples);
allMarkerNames = allMarkerNames(sortIdx);
allMarkerIdx = allMarkerIdx(sortIdx);
allMarkerDurations = allMarkerDurations(sortIdx);

nSegments = numel(allMarkerSamples);
preBufferSamples = round(taskConfig.preBuffer * obj.sample_freq);
postBufferSamples = round(taskConfig.postBuffer * obj.sample_freq);
totalSamples = size(obj.elec_data, 2);

segmentFiles = cell(nSegments, 1);
rawSegmentFiles = cell(nSegments, 1);
segmentDurations = zeros(nSegments, 1);

fprintf('\n=== NATURALISTIC: EXTRACT SEGMENTS (%d) ===\n', nSegments);
fprintf('Segmentation mode: %s\n', taskConfig.segmentMode);

for s = 1:nSegments
    segmentName = allMarkerNames{s};
    segmentInstance = allMarkerIdx(s);
    markerSample = allMarkerSamples(s);
    configuredDuration = allMarkerDurations(s);

    startSample = markerSample - preBufferSamples;
    if strcmp(taskConfig.segmentMode, 'duration')
        durationSamples = round(configuredDuration * obj.sample_freq);
        endSample = markerSample + durationSamples - 1 + postBufferSamples;
    elseif s < nSegments
        endSample = allMarkerSamples(s+1) - 1;
    else
        endSample = totalSamples;
    end

    startSample = max(1, startSample);
    endSample = min(totalSamples, endSample);

    actualDuration = (endSample - startSample + 1) / obj.sample_freq;
    segmentDurations(s) = actualDuration;

    fprintf('\nSegment %d/%d: %s (instance %d), %.1f s\n', ...
        s, nSegments, segmentName, segmentInstance, actualDuration);

    meta = struct();
    meta.subject = params.SubjectName;
    meta.marker = segmentName;
    meta.instance = segmentInstance;
    meta.startSample = startSample;
    meta.endSample = endSample;
    meta.markerSample = markerSample;
    meta.configuredDuration = configuredDuration;
    meta.actualDuration = actualDuration;
    meta.preBuffer = taskConfig.preBuffer;
    meta.postBuffer = taskConfig.postBuffer;
    meta.segmentMode = taskConfig.segmentMode;
    meta.sample_freq = obj.sample_freq;

    segment = meta;
    segment.signalType = 'high_gamma';
    segment.elec_data = obj.elec_data(:, startSample:endSample);
    segment.elec_ch_label = obj.elec_ch_label;
    segment.elec_ch_clean = obj.elec_ch_clean;
    segment.elec_ch_type = obj.elec_ch_type;
    if ~isempty(obj.bip_elec_data)
        segment.bip_elec_data = obj.bip_elec_data(:, startSample:endSample);
        segment.bip_ch_label = obj.bip_ch_label;
        segment.bip_ch = obj.bip_ch;
    else
        segment.bip_elec_data = [];
        segment.bip_ch_label = {};
        segment.bip_ch = [];
    end
    if isfield(obj, 'anatomy') && ~isempty(obj.anatomy)
        segment.anatomy = obj.anatomy;
    end
    segment.preprocessing = obj.for_preproc;

    hgFilename = sprintf('%s_%s_%02d_%.0fsec_HG.mat', ...
        params.SubjectName, segmentName, segmentInstance, actualDuration);
    hgPath = fullfile(segmentDir, hgFilename);
    save(hgPath, 'segment', '-v7.3');
    segmentFiles{s} = hgPath;
    fprintf('  [HG] %s\n', hgFilename);

    if saveBroadbandSegments
        bbStart = map_samples_between_rates(startSample, hg_sample_freq, bb_sample_freq);
        bbEnd = map_samples_between_rates(endSample, hg_sample_freq, bb_sample_freq);
        bbEnd = min(bbEnd, size(bb_elec_data, 2));
        if bbEnd < bbStart
            warning('Segment %d (%s): broadband range empty after rate mapping.', s, segmentName);
            bbStart = min(bbStart, size(bb_elec_data, 2));
            bbEnd = bbStart;
        end

        rawSegment = meta;
        rawSegment.signalType = 'broadband_raw';
        rawSegment.sample_freq = bb_sample_freq;
        rawSegment.startSample = bbStart;
        rawSegment.endSample = bbEnd;
        rawSegment.markerSample = map_samples_between_rates(markerSample, hg_sample_freq, bb_sample_freq);
        rawSegment.actualDuration = (bbEnd - bbStart + 1) / bb_sample_freq;
        rawSegment.elec_data = bb_elec_data(:, bbStart:bbEnd);
        rawSegment.elec_ch_label = obj.elec_ch_label;
        rawSegment.elec_ch_clean = obj.elec_ch_clean;
        rawSegment.elec_ch_type = obj.elec_ch_type;
        if ~isempty(bb_bip_elec_data)
            rawSegment.bip_elec_data = bb_bip_elec_data(:, bbStart:bbEnd);
            rawSegment.bip_ch_label = obj.bip_ch_label;
            rawSegment.bip_ch = obj.bip_ch;
        else
            rawSegment.bip_elec_data = [];
            rawSegment.bip_ch_label = {};
            rawSegment.bip_ch = [];
        end
        if isfield(obj, 'anatomy') && ~isempty(obj.anatomy)
            rawSegment.anatomy = obj.anatomy;
        end
        rawSegment.preprocessing = obj.for_preproc;

        rawFilename = sprintf('%s_%s_%02d_%.0fsec_broadband.mat', ...
            params.SubjectName, segmentName, segmentInstance, actualDuration);
        rawPath = fullfile(rawSegmentDir, rawFilename);
        save(rawPath, 'rawSegment', '-v7.3');
        rawSegmentFiles{s} = rawPath;
        fprintf('  [BB] %s\n', rawFilename);
    end
end

%% Summary
fprintf('\n=== NATURALISTIC: SUMMARY ===\n');
summary = struct();
summary.subject = params.SubjectName;
summary.taskType = params.taskType;
summary.processDate = datetime('now');
summary.segmentMode = taskConfig.segmentMode;
summary.nSegments = nSegments;
summary.markerNames = allMarkerNames;
summary.markerInstances = allMarkerIdx;
summary.markerSamples = allMarkerSamples;
summary.configuredDurations = allMarkerDurations;
summary.actualDurations = segmentDurations;
summary.preBuffer = taskConfig.preBuffer;
summary.postBuffer = taskConfig.postBuffer;
summary.segmentFiles = segmentFiles;
summary.rawSegmentFiles = rawSegmentFiles;
summary.rawSegmentDir = rawSegmentDir;
summary.nChannels_unipolar = numel(obj.elec_ch_clean);
summary.nChannels_bipolar = size(obj.bip_elec_data, 1);
summary.sampleRate = obj.sample_freq;
summary.preprocessing = obj.for_preproc;
summary.crunchedFile = crunchedFile;
summary.taskConfig = taskConfig;

summaryFile = fullfile(outputDir, sprintf('%s_%s_summary.mat', params.SubjectName, params.taskType));
save(summaryFile, 'summary', '-v7.3');

summaryTxtFile = fullfile(outputDir, sprintf('%s_%s_summary.txt', params.SubjectName, params.taskType));
write_naturalistic_summary_txt(summaryTxtFile, summary, taskConfig, segmentFiles, rawSegmentFiles);

save(crunchedFile, 'obj', 'naturalEvents', 'summary', '-v7.3');

fprintf('Summary saved: %s\n', summaryFile);
fprintf('Segment directory: %s\n', segmentDir);
if saveBroadbandSegments
    fprintf('Broadband segment directory: %s\n', rawSegmentDir);
end
end

function obj = ensure_naturalistic_hg_ready(obj, decimationFreq)
if abs(obj.sample_freq - decimationFreq) <= 0.1
    return;
end

rawFs = obj.sample_freq;
if isfield(obj.for_preproc, 'sample_freq_raw') && ~isempty(obj.for_preproc.sample_freq_raw)
    rawFs = obj.for_preproc.sample_freq_raw;
end

if abs(obj.sample_freq - rawFs) <= 0.1
    fprintf('Naturalistic: extracting high-gamma envelope...\n');
    obj.extract_high_gamma('doNapLabFilterExtraction', true);
end

if abs(obj.sample_freq - decimationFreq) > 0.1
    fprintf('Naturalistic: downsampling HG to %d Hz...\n', decimationFreq);
    obj.downsample_signal('decimationFreq', decimationFreq);
end
end

function sample = map_samples_between_rates(sample, srcFs, tgtFs)
if abs(srcFs - tgtFs) <= 0.1
    return;
end
sample = max(1, round((sample - 1) * tgtFs / srcFs) + 1);
end

function naturalEvents = adjust_natural_events_for_decimation(obj, naturalEvents)
originalFs = obj.for_preproc.sample_freq_raw;
currentFs = obj.sample_freq;

if abs(currentFs - originalFs) <= 0.1
    fprintf('No downsampling detected (%.0f Hz). Event samples unchanged.\n', currentFs);
    return;
end

fprintf('Downsampled %.0f Hz -> %.0f Hz; adjusting event samples...\n', originalFs, currentFs);
stitch_raw = obj.for_preproc.stitch_index_raw;
stitch_dec = obj.for_preproc.stitch_index_dec;

for i = 1:numel(naturalEvents.samples)
    rawSample = naturalEvents.samples(i);
    runIdx = find(rawSample >= stitch_raw, 1, 'last');
    rawOffset = stitch_raw(runIdx) - 1;
    decOffset = stitch_dec(runIdx) - 1;
    newSample = decOffset + round((rawSample - rawOffset - 1) * currentFs / originalFs) + 1;
    naturalEvents.samples(i) = newSample;
    naturalEvents.times(i) = newSample / currentFs;
end
end

function write_naturalistic_summary_txt(summaryTxtFile, summary, taskConfig, segmentFiles, rawSegmentFiles)
fid = fopen(summaryTxtFile, 'w');
if fid < 0
    warning('Could not write summary text file: %s', summaryTxtFile);
    return;
end

cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'NATURALISTIC PREPROCESSING SUMMARY\n');
fprintf(fid, '===================================\n\n');
fprintf(fid, 'Subject: %s\n', summary.subject);
fprintf(fid, 'Task Type: %s\n', summary.taskType);
fprintf(fid, 'Process Date: %s\n', char(summary.processDate));
fprintf(fid, 'Segmentation Mode: %s\n\n', summary.segmentMode);
fprintf(fid, 'Number of Segments: %d\n', summary.nSegments);
fprintf(fid, 'Unipolar Channels: %d\n', summary.nChannels_unipolar);
fprintf(fid, 'Bipolar Channels: %d\n', summary.nChannels_bipolar);
fprintf(fid, 'Sample Rate: %.1f Hz\n\n', summary.sampleRate);

fprintf(fid, 'CONFIGURED SEGMENT DURATIONS:\n');
fprintf(fid, '-----------------------------\n');
for m = 1:numel(taskConfig.markers)
    markerName = taskConfig.markers{m};
    duration = taskConfig.segmentDurations.(markerName);
    fprintf(fid, '  %s: %.1f seconds (%.1f minutes)\n', markerName, duration, duration/60);
end

fprintf(fid, '\nHIGH-GAMMA SEGMENT FILES:\n');
fprintf(fid, '-------------------------\n');
for s = 1:numel(segmentFiles)
    [~, fname, ext] = fileparts(segmentFiles{s});
    fprintf(fid, '%d. %s (instance %d)\n', s, summary.markerNames{s}, summary.markerInstances(s));
    fprintf(fid, '   Configured: %.1f sec | Actual: %.1f sec\n', ...
        summary.configuredDurations(s), summary.actualDurations(s));
    fprintf(fid, '   File: %s%s\n\n', fname, ext);
end

if ~isempty(rawSegmentFiles) && ~isempty(rawSegmentFiles{1})
    fprintf(fid, '\nRAW (BROADBAND) SEGMENT FILES:\n');
    fprintf(fid, '------------------------------\n');
    for s = 1:numel(rawSegmentFiles)
        [~, fname, ext] = fileparts(rawSegmentFiles{s});
        fprintf(fid, '%d. %s%s\n', s, fname, ext);
    end
end

fprintf(fid, '\nSource file: %s\n', summary.crunchedFile);
end

function val = get_naturalistic_opt(opts, name, defaultVal)
if isfield(opts, name) && ~isempty(opts.(name))
    val = opts.(name);
else
    val = defaultVal;
end
end
