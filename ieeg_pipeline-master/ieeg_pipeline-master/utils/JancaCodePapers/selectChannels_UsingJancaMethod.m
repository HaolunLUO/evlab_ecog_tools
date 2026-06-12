
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JB Eichenlaub - Last modification 01.12.2017
% Main pre-processing step:
% 	- combine results from Janca approach across segments uing discharges.MV and d_decim
%   - compute IEDs rate
%   - and aplly threshold
% Output: tableChanSelection containing number of IEDs per chan, chans below threshold, ...

function [tableChanSelection] = selectChannels_UsingJancaMethod...
    (IEDsDetected, threshold, nameChanIni, settings, patientID, segmentsID,TimeRangeBefore,TimeRangeAfter)

% nameChanIni=nameChans.Ini
% IEDsDetected=detectionIEDs.segments
% threshold=detectionIEDs.threshold
% segmentsID= length(segmentsID);
numChanIni = size(nameChanIni,1);
  
% Combine segments to assess overall IEDs across all available segments
% dischages.MV: [n x c] matrix where n number of IEDs detected and c channels
% d_decim:      [s x c] matrix where s number of samples (after dowsnsampling) and c channels
% 0=no detected, 1=detected. When k2 ON, ambiguous IEDs are scored 0.5. Default k2 OFF.
currIEDs.fs             = 200; % default downsampling during automatic detection (-dec 200)
currIEDs.discharges.MV  = [];
currIEDs.numSamples     = 0;
for seg = 1:numel(segmentsID)
    if ~isempty(IEDsDetected{seg})
        currIEDs.discharges.MV = [currIEDs.discharges.MV; IEDsDetected{seg}.discharges.MV];
        currIEDs.numSamples    = currIEDs.numSamples + size(IEDsDetected{seg}.d_decim, 1);
    end
end

% Sanity check
 if numChanIni ~= size(currIEDs.discharges.MV, 2); error(['Either: number of channel labels vs channels in the chan x time sample differ!', ...
     'OR NO detected IIDs, data could be flat! I am kind of a silly and confusing warning message!']); end

% Compute number of detected spike per channel: [c x 1] where c channels
numSpikes     = []; numSpikes     = sum(currIEDs.discharges.MV==1, 1);
totalDuration = []; totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60;
numSpikes_min = []; numSpikes_min = numSpikes / totalDuration;
numSpikes     = transpose(numSpikes);
numSpikes_min = transpose(numSpikes_min);

% Select channels with IEDs / minute below threshold - [c x 1] where c channels
indChanSelected = [];
indChanSelected = find(numSpikes_min < threshold);
tableChanSelection.numSpikesAll           = numSpikes_min;
tableChanSelection.indChansSelected       = indChanSelected;
tableChanSelection.nameChansSelected      = nameChanIni(indChanSelected,:);
tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);

% plot figure with threshold
figure; set(gcf, 'Position', [200, 200, 1700, 600]);
numSpikes = []; numSpikes = tableChanSelection.numSpikesAll;
subplot(2,2,1)
if isempty(TimeRangeBefore)==0
%     plot(TimeRangeBefore')
%         title('voltages before task')
end
box off
subplot(2,2,2)
if isempty(TimeRangeAfter)==0
%     plot(TimeRangeAfter')
%     title('voltages after task')
end
box off
subplot(2,2,3:4)
plot(numSpikes); hold on; plot([0 numChanIni+1], [threshold threshold] , 'black');
set(gca, 'XTick', 1:numChanIni, 'XTickLabel', nameChanIni(1:numChanIni,:), 'FontSize', 10, 'XTickLabelRotation', 45);
xlim([0 numChanIni+1]); ylim([0 20]);
ylabel(['Num IEDs detected /minute']);
box off
title(['Channel Selection - Patient: ' patientID ' - setting: ' settings ' - Threshold = ' num2str(threshold)]);
