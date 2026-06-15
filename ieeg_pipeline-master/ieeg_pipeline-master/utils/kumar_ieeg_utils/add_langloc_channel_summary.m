function add_langloc_channel_summary(rpt, obj)
% ADD_LANGLOC_CHANNEL_SUMMARY  Langloc significant-channel section for PDF reports.
%
%   Uses word-boundary results when present; otherwise falls back to the
%   sentence-vs-nonword contrast stored by high_gamma_plot_langloc.

import mlreportgen.report.*
import mlreportgen.dom.*

if isfield(obj.stats.time_series, 'pSigChan_wordboundaries_langloc')
    add(rpt, Heading2('Unipolar Channels with Significant Time Clusters (word boundaries)'));
    sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), ...
        obj.stats.time_series.pSigChan_wordboundaries_langloc));
    sigUnipolarChannels = find(sigChannels > 0);
elseif isfield(obj.stats.time_series, 'pSigChan_langloc_All_Trials')
    add(rpt, Heading2('Unipolar Channels with Significant Sentence vs Nonword Contrast'));
    pSig = obj.stats.time_series.pSigChan_langloc_All_Trials;
    sigUnipolarChannels = find(cellfun(@(x) isstruct(x) && any(x.h_sig_05), pSig));
else
    add(rpt, Paragraph('Langloc channel summary not available (word-boundary analysis was not run).'));
    sigUnipolarChannels = [];
end

if ~isempty(sigUnipolarChannels)
    unipolarList = cell(length(sigUnipolarChannels), 1);
    for i = 1:length(sigUnipolarChannels)
        unipolarList{i} = obj.elec_ch_label{sigUnipolarChannels(i)};
    end
    add(rpt, UnorderedList(unipolarList));
    add(rpt, PageBreak());
else
    add(rpt, Paragraph('No significant langloc unipolar channels found.'));
end

if isfield(obj.stats.time_series, 'pSigChan_bip_wordboundaries_langloc')
    add(rpt, Heading2('Bipolar Channels with Significant Time Clusters (word boundaries)'));
    sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), ...
        obj.stats.time_series.pSigChan_bip_wordboundaries_langloc));
    sigBipolarChannels = find(sigChannels > 0);
elseif isfield(obj.stats.time_series, 'pSigChan_bip_langloc_All_Trials')
    add(rpt, Heading2('Bipolar Channels with Significant Sentence vs Nonword Contrast'));
    pSigBip = obj.stats.time_series.pSigChan_bip_langloc_All_Trials;
    sigBipolarChannels = find(cellfun(@(x) isstruct(x) && any(x.h_sig_05), pSigBip));
else
    sigBipolarChannels = [];
end

if ~isempty(sigBipolarChannels)
    bipolarList = cell(length(sigBipolarChannels), 1);
    for i = 1:length(sigBipolarChannels)
        bipolarList{i} = obj.bip_ch_label{sigBipolarChannels(i)};
    end
    add(rpt, UnorderedList(bipolarList));
    add(rpt, PageBreak());
elseif isfield(obj.stats.time_series, 'pSigChan_bip') || ...
        isfield(obj.stats.time_series, 'pSigChan_bip_langloc_All_Trials') || ...
        isfield(obj.stats.time_series, 'pSigChan_bip_wordboundaries_langloc')
    add(rpt, Paragraph('No significant langloc bipolar channels found.'));
end
end
