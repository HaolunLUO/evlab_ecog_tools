function naturalEvents = extract_naturalistic_events(allDataFiles, taskConfig, obj)
% EXTRACT_NATURALISTIC_EVENTS  Scan Brainstorm exports for NA/NA2/NA3 markers.
%
%   naturalEvents = extract_naturalistic_events(allDataFiles, taskConfig, obj)
%
%   taskConfig.markers must list the event labels to collect (e.g. NA, NA2, NA3).
%   Event sample indices are returned in the coordinate system of the concatenated
%   raw recording at obj.sample_freq (before any downsampling).

if nargin < 3 || isempty(obj)
    error('obj with sample_freq is required.');
end
if ~isfield(taskConfig, 'markers') || isempty(taskConfig.markers)
    error('taskConfig.markers is required.');
end

markers = taskConfig.markers;
naturalEvents = struct('labels', {{}}, 'samples', [], 'times', []);

fprintf('Extracting naturalistic events from %d Brainstorm file(s)...\n', numel(allDataFiles));

cumulativeSamples = 0;
sfreq = obj.sample_freq;

for f = 1:numel(allDataFiles)
    fprintf('  Scanning file %d/%d: %s\n', f, numel(allDataFiles), allDataFiles{f});
    bstData = load(allDataFiles{f});

    if f == 1
        if isfield(bstData, 'F') && isstruct(bstData.F) && isfield(bstData.F, 'prop') ...
                && isfield(bstData.F.prop, 'sfreq')
            sfreq = bstData.F.prop.sfreq;
        elseif isfield(bstData, 'Fs')
            sfreq = bstData.Fs;
        end
    end

    if isfield(bstData, 'F') && isnumeric(bstData.F)
        nSamplesInFile = size(bstData.F, 2);
    elseif isfield(bstData, 'Time')
        nSamplesInFile = numel(bstData.Time);
    else
        nSamplesInFile = 0;
    end

    if isfield(bstData, 'Time') && ~isempty(bstData.Time)
        fileStartTime = bstData.Time(1);
    else
        fileStartTime = 0;
    end

    if isfield(bstData, 'Events') && ~isempty(bstData.Events)
        for e = 1:numel(bstData.Events)
            evt = bstData.Events(e);
            evtLabel = evt.label;
            if ~any(strcmpi(evtLabel, markers))
                continue;
            end

            fprintf('    Found marker: %s\n', evtLabel);

            if isfield(evt, 'times') && ~isempty(evt.times)
                evtTimes = evt.times;
            elseif isfield(evt, 'samples') && ~isempty(evt.samples)
                evtTimes = evt.samples / sfreq;
            else
                continue;
            end

            if size(evtTimes, 1) > 1
                evtTimes = evtTimes(1, :);
            end
            evtTimes = evtTimes(:)';

            for t = 1:numel(evtTimes)
                evtTime = evtTimes(t);
                sampleInFile = round((evtTime - fileStartTime) * sfreq) + 1;
                globalSample = cumulativeSamples + sampleInFile;
                globalTime = globalSample / obj.sample_freq;

                naturalEvents.labels{end+1} = evtLabel; %#ok<AGROW>
                naturalEvents.samples(end+1) = globalSample; %#ok<AGROW>
                naturalEvents.times(end+1) = globalTime; %#ok<AGROW>

                fprintf('      Event at %.2f sec -> global sample %d (%.2f sec)\n', ...
                    evtTime, globalSample, globalTime);
            end
        end
    else
        fprintf('    No Events field found.\n');
    end

    cumulativeSamples = cumulativeSamples + nSamplesInFile;
end

if isempty(naturalEvents.labels)
    error('No naturalistic markers (%s) found in data files.', strjoin(markers, ', '));
end

[naturalEvents.samples, sortIdx] = sort(naturalEvents.samples);
naturalEvents.labels = naturalEvents.labels(sortIdx);
naturalEvents.times = naturalEvents.times(sortIdx);

fprintf('Found %d naturalistic event(s).\n', numel(naturalEvents.labels));
end
