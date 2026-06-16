function plot_crosscheck_channels(sn_obj, overlapIdx, montage, maxPlots)
% PLOT_CROSSCHECK_CHANNELS  Plot HG z-score traces for overlap channels.
%
%   plot_crosscheck_channels(sn_obj, overlapIdx, montage)
%   plot_crosscheck_channels(sn_obj, overlapIdx, montage, maxPlots)
%
%   overlapIdx - channel indices (unipolar or bipolar, per montage)
%   montage    - 'unipolar' or 'bipolar'

if nargin < 4 || isempty(maxPlots)
    maxPlots = 6;
end

if isempty(overlapIdx)
    fprintf('No overlap channels to plot.\n');
    return;
end

montage = lower(montage);
nPlot = min(numel(overlapIdx), maxPlots);

if strcmp(montage, 'unipolar')
    data = sn_obj.elec_data;
    labels = sn_obj.elec_ch_label;
    artStats = sn_obj.stats.artifact_stats_unipolar;
else
    data = sn_obj.bip_elec_data;
    labels = sn_obj.bip_ch_label;
    artStats = sn_obj.stats.artifact_stats_bipolar;
end

fs = sn_obj.sample_freq;
t = (1:size(data, 2)) / fs;
artCount = artifact_counts_for_plot(artStats, size(data, 1));

for i = 1:nPlot
    ch = overlapIdx(i);
    figure('Name', sprintf('%s overlap ch %d', montage, ch), 'Color', 'w');
    plot(t, data(ch, :), 'Color', [0.2 0.2 0.8]);
    grid on;
    xlabel('Time (s)');
    ylabel('HG z-score');
    title(sprintf('%s | nArtifacts=%d', labels{ch}, artCount(ch)));

    statIdx = find([artStats.channel_idx] == ch, 1, 'first');
    if ~isempty(statIdx) && artStats(statIdx).artifact_count > 0 ...
            && ~isempty(artStats(statIdx).artifact_times)
        hold on;
        yl = ylim;
        evtT = artStats(statIdx).artifact_times;
        for j = 1:numel(evtT)
            plot([evtT(j) evtT(j)], yl, '--', 'Color', [0.85 0.2 0.2]);
        end
    end
end
end


function artCount = artifact_counts_for_plot(artStats, nChannels)
artCount = zeros(nChannels, 1);
for k = 1:numel(artStats)
    idx = artStats(k).channel_idx;
    if idx >= 1 && idx <= nChannels
        artCount(idx) = artStats(k).artifact_count;
    end
end
end
