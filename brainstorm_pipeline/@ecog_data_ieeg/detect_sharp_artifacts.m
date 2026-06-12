function detect_sharp_artifacts(obj, ops)
% Detect sharp inter-ictal transients on significant bipolar *and* unipolar
% channels. When the dataset comes from the MIT Naturalistic Stories task
% the detector runs **only inside the story epochs**.
% 
% Modified for z-scored high gamma envelope data - removes sharp spiky artifacts
% by analyzing amplitude changes without frequency normalization.
% Uses OR logic: keeps peaks that pass EITHER amplitude OR slope criteria.

% ─────────────────────────── inputs ────────────────────────────
arguments
    obj ecog_data_ieeg
    ops.min_amplitude double = 15   % z-score units
    ops.min_slope     double = 10   % z-score units per sample
    
end
fs = obj.sample_freq;

% ──────────────────── decide which samples to keep ─────────────
isStories = isprop(obj,'experiment') && ...
            strcmp(obj.experiment,'MITNaturalisticStoriesTask');

if isStories
    % collect start/end indices of every story trial
    epochs = [];                         % [nEpoch × 2] (start, end)
    for t = 1:size(obj.trial_timing,1)
        tbl = obj.trial_timing{t,1};
        if startsWith(tbl.key,"story_")      % story_1 … story_4 …
            epochs = [epochs ; tbl.start tbl.end]; %#ok<AGROW>
        end
    end
    if isempty(epochs)
        warning('No story epochs found – analysing full recording instead');
        epochMask = true(1,size(obj.elec_data,2));
    else
        epochMask = false(1,size(obj.elec_data,2));
        for e = 1:size(epochs,1)
            s = max(1, floor(epochs(e,1)));
            f = min(length(epochMask), ceil(epochs(e,2)));
            epochMask(s:f) = true;
        end
    end
else
    epochMask = true(1,size(obj.elec_data,2));    % keep everything
end

% ───────────────────────── run detector ───────────────────────
if ~isempty(obj.bip_elec_data)
    obj.stats.artifact_stats_bipolar  = ...
        run_detector(obj.bip_elec_data(:,epochMask), ...
                     fs, ops);
end

obj.stats.artifact_stats_unipolar = ...
    run_detector(obj.elec_data(:,epochMask), ...
                 fs, ops);

% store mask for downstream plotting
obj.stats.epochMask = epochMask;

% ────────────────────────── figures ───────────────────────────
% plot_artifact_summary(obj);
% plot_artifact_timeseries(obj);       % now story-restricted
end

function stats = run_detector(data, fs, ops)
% Detect sharp artefacts on z-scored high gamma envelope data
% Uses OR logic: detects peaks that satisfy EITHER slope OR amplitude criteria
% ----------------------------------------------------------
% Input
%   data     : [nChan × nSamples] z-scored high gamma envelope traces
%   fs       : sampling frequency (Hz)
%   ops.*    : detection thresholds (already validated upstream)
%
% Output
%   stats(k) – struct with per-channel fields:
%       • channel_idx     – original index in the full montage
%       • artifact_count  – number of detected events
%       • max_amplitude   – largest peak amplitude (z-score) or NaN
%       • mean_slope      – mean absolute slope of accepted peaks (z-score/sample)
%       • artifact_times  – 1×M vector of event times (s)

% ---------------- Channel selection ------------------------
sigData  = data;
chan_ids = 1:size(data,1);
n        = numel(chan_ids);

% ---------------- Pre-allocate struct ----------------------
tmpl  = struct('channel_idx',[],'artifact_count',0,'max_amplitude',NaN,...
               'mean_slope',NaN,'artifact_times',[]);
stats = repmat(tmpl,n,1);

% ---------------- Main loop -------------------------------
parfor k = 1:n
    x  = sigData(k,:);
    dx = diff(x);                            % slope (z-score/sample) - NO fs division

    % --- METHOD 1: slope-based detection -------------------
    [slopePkAmp, slopePkLoc] = findpeaks(abs(dx), ...
                     'MinPeakHeight', ops.min_slope);

    % --- METHOD 2: amplitude-based detection ---------------
    [ampPkAmp, ampPkLoc] = findpeaks(abs(x), ...
                     'MinPeakHeight', ops.min_amplitude);

    % --- combine and deduplicate peaks (OR logic) ----------
    allLoc = [slopePkLoc, ampPkLoc];
    allAmp = [slopePkAmp, ampPkAmp];
    
    if ~isempty(allLoc)
        % Remove duplicates (peaks detected by both methods)
        [uniqueLoc, uniqueIdx] = unique(allLoc);
        uniqueAmp = allAmp(uniqueIdx);
        
        % Sort by location
        [keepLoc, sortIdx] = sort(uniqueLoc);
        keepAmp = uniqueAmp(sortIdx);
        
        % For slope values, we need to get the actual slope at each location
        keepSlope = abs(dx(keepLoc));
    else
        keepLoc = [];
        keepAmp = [];
        keepSlope = [];
    end

    % --- store per-channel results -------------------------
    stats(k).channel_idx     = chan_ids(k);
    stats(k).artifact_count  = numel(keepLoc);
    stats(k).artifact_times  = keepLoc/fs;

    % ensure scalar output
    if isempty(keepLoc)
        stats(k).max_amplitude = NaN;
        stats(k).mean_slope    = NaN;
    else
        stats(k).max_amplitude = max(abs(x(keepLoc)));
        stats(k).mean_slope    = mean(keepSlope,'omitnan');
    end
end
end

function plot_artifact_summary(obj)
if ~isempty(obj.bip_elec_data)
    b = obj.stats.artifact_stats_bipolar;
end
u = obj.stats.artifact_stats_unipolar;

figure('Name','Sharp-artifact overview (OR logic)','Color','w');
if ~isempty(obj.bip_elec_data)
    % --- Bipolar plots ---
    subplot(2,2,1)
    plot([b.channel_idx],[b.artifact_count],'o-'); grid on
    title('Bipolar – event count'); xlabel('channel'); ylabel('#artifacts');
    
    subplot(2,2,2)
    plot([u.channel_idx],[u.artifact_count],'o-'); grid on
    title('Unipolar – event count'); xlabel('channel'); ylabel('#artifacts');
    
    % --- Bipolar scatter (with NaN filtering) ---
    subplot(2,2,3)
    b_slope = [b.mean_slope];
    b_amp   = [b.max_amplitude];
    b_count = [b.artifact_count];
    
    % Remove NaN values and ensure matching lengths
    valid_b = ~isnan(b_slope) & ~isnan(b_amp);
    if any(valid_b)
        scatter(b_slope(valid_b), b_amp(valid_b), 40, b_count(valid_b), 'filled');
        colorbar; 
        if max(b_count) > 0
            clim([0 max(b_count)]);
        end
    end
    grid on
    title('Bipolar: slope vs amplitude'); 
    xlabel('mean slope (z-score/sample)'); ylabel('max amp (z-score)');
    
    % --- Unipolar scatter (with NaN filtering) ---
    subplot(2,2,4)
    u_slope = [u.mean_slope];
    u_amp   = [u.max_amplitude];
    u_count = [u.artifact_count];
    
    % Remove NaN values and ensure matching lengths
    valid_u = ~isnan(u_slope) & ~isnan(u_amp);
    if any(valid_u)
        scatter(u_slope(valid_u), u_amp(valid_u), 40, u_count(valid_u), 'filled');
        colorbar;
        if max(u_count) > 0
            clim([0 max(u_count)]);
        end
    end
    grid on
    title('Unipolar: slope vs amplitude'); 
    xlabel('mean slope (z-score/sample)'); ylabel('max amp (z-score)');
else
    subplot(1,2,1)
    plot([u.channel_idx],[u.artifact_count],'o-'); grid on
    title('Unipolar – event count'); xlabel('channel'); ylabel('#artifacts');
    
    % --- Unipolar scatter (with NaN filtering) ---
    subplot(1,2,2)
    u_slope = [u.mean_slope];
    u_amp   = [u.max_amplitude];
    u_count = [u.artifact_count];
    
    % Remove NaN values and ensure matching lengths
    valid_u = ~isnan(u_slope) & ~isnan(u_amp);
    if any(valid_u)
        scatter(u_slope(valid_u), u_amp(valid_u), 40, u_count(valid_u), 'filled');
        colorbar;
        if max(u_count) > 0
            clim([0 max(u_count)]);
        end
    end
    grid on
    title('Unipolar: slope vs amplitude'); 
    xlabel('mean slope (z-score/sample)'); ylabel('max amp (z-score)');
end
end

function plot_artifact_timeseries(obj)
u = obj.stats.artifact_stats_unipolar;

fs        = obj.sample_freq;
emask     = obj.stats.epochMask;
t         = find(emask)/fs;                    % x-axis (s)

% ─ bipolar ─
if ~isempty(obj.bip_elec_data)
    b = obj.stats.artifact_stats_bipolar;
    plot_group(obj.bip_elec_data, b, true, 'Bipolar – channels WITH artifacts');
    plot_group(obj.bip_elec_data, b, false, 'Bipolar – channels WITHOUT artifacts');
end

% ─ unipolar ─
plot_group(obj.elec_data, u, true, 'Unipolar – channels WITH artifacts');
plot_group(obj.elec_data, u, false, 'Unipolar – channels WITHOUT artifacts');

    function plot_group(data, stats, withArt, figName)
        sel = ([stats.artifact_count] > 0);
        if ~withArt,  sel = ~sel;  end
        ch  = [stats(sel).channel_idx];
        if isempty(ch), return, end

        figure('Name',figName,'Color','w');
        plot(t, data(ch,emask).'); grid on
        xlabel('Time (s)'); ylabel('z-score');
        title(sprintf('%s (n = %d)', figName, numel(ch)));
    end
end
