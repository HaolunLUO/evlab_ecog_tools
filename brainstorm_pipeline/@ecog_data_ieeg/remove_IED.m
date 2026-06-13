function remove_IED(obj)
    % This function marks electrodes with significant Interictal 
    % Epileptiform Discharges (IEDs) using a protocol adopted from 
    % Radek Janca (see directory below for sample papers).
    % Once marked, these electrodes will be removed from subsequent 
    % analyses.
    %
    % This is the only function in the pipeline where the parameters
    % are hardcoded into the function itself instead of defined in the 
    % define_parameters function. That is because this script was 
    % designed by another lab and should only be modified with extreme
    % caution. Similarly, this is the only part of the pipeline that 
    % calls on functions not contained within this ecog_data.m file. 
    % This is, again, because the procedure was designed by another 
    % lab and it was determined that their script should be self-
    % contained. 
    %
    % NOTE - the following directory must be added to your MATLAB path
    % /mindhive/evlab/u/Shared/merged_ecog_pipeline/utils/JancaCodePapers 

    signal = double(obj.elec_data');

    fprintf(1, '\n> Finding electrodes with significant Interictal Epileptiform Discharges (IEDs) ... \n');
    
    detectionIEDs          = []; % output from Janca et al. script - 1 cell array per segment
    detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1'; % if you change "-dec 200" here, do not forget to change in selectChannels_Using ... below
    detectionIEDs.segments = [];
        
    % NOT DEFINED IN THIS CLASS
    detectionIEDs = automaticSpikeDetection_UsingJancaMethod(signal,obj.sample_freq,detectionIEDs.settings);
    
    if obj.for_preproc.isPlotVisible
        obj.plot_channels(detectionIEDs.envelope,...
                        obj.elec_ch_label,...
                        obj.elec_ch_clean,...
                        obj.elec_ch_valid,...
                        't_len',100,...
                        'sample_freq',200,...
                        'plotIEDs',true,...
                        'chanIEDs',detectionIEDs.out.chan,...
                        'posIEDs',detectionIEDs.out.pos...
        ); 
    end
        
    % From this automatic assessment, selection of the final pool of channels
    detectionIEDs.tableChanSelection = [];  % info regarding chan selection
    detectionIEDs.threshold = 6.5; % channels with IEDs higher than threshold are removed  
        
    currIEDs.fs = 200; % default downsampling during automatic detection (-dec 200)
    currIEDs.discharges.MV = [];
    currIEDs.numSamples = 0;
        
    currIEDs.discharges.MV = detectionIEDs.discharges.MV;
    currIEDs.numSamples= currIEDs.numSamples + size(detectionIEDs.d_decim, 1);
        
    % Compute number of detected spike per channel: [c x 1] where c channels
    numSpikes     = []; numSpikes     = sum(currIEDs.discharges.MV==1, 1);
    totalDuration = []; totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60; % in minutes
    numSpikes_min = []; numSpikes_min = numSpikes / totalDuration;
    numSpikes     = transpose(numSpikes);
    numSpikes_min = transpose(numSpikes_min);
        
    % Select channels with IEDs / minute below threshold - [c x 1] where c channels
    indChanSelected = [];
    indChanSelected = find(numSpikes_min < detectionIEDs.threshold);
    tableChanSelection.numSpikesAll           = numSpikes_min;
    tableChanSelection.indChansSelected       = indChanSelected;
    tableChanSelection.indChansDeselected     = setdiff(obj.elec_ch,indChanSelected);
    tableChanSelection.nameChansSelected      = transpose(obj.elec_ch_label(indChanSelected));
    tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);
        
        
    obj.for_preproc.IEDRemoval_results=tableChanSelection;
    obj.for_preproc.IEDRemoval_results.threshold=detectionIEDs.threshold;

    

    % Display total number of electrodes with IEDs.
    % NOTE: obj.elec_ch is a column vector ([nChannels x 1], created as
    % (1:nChannels)' in the conversion), so the total channel count must be
    % taken with numel(). Using size(obj.elec_ch,2) returned 1 regardless of
    % how many channels there were, which made the skip threshold ceil(1/3)=1
    % and incorrectly aborted IED removal (e.g. "16 / 1 (1600.0%)").
    nIED    = length(tableChanSelection.indChansDeselected);
    nTotal  = numel(obj.elec_ch);
    pctIED  = 100 * nIED / nTotal;
    thresh  = ceil(nTotal / 3);

    fprintf('  Electrodes with IEDs:  %d / %d (%.1f%%) — threshold: >%d\n', ...
        nIED, nTotal, pctIED, thresh);

    % Single decision point: skip the step if too many electrodes were flagged,
    % otherwise mark the IED electrodes and recompute the clean-channel set.
    if nIED > thresh
        fprintf('  *** Too many electrodes with significant IEDs, SKIPPING STEP ***\n');

        new_order_mask = cell2mat(cellfun(@(x) strcmp(x,'IEDRemoval'),obj.for_preproc.order,'UniformOutput',false));
        obj.for_preproc.order = obj.for_preproc.order(~new_order_mask);

    else
        fprintf('  IED count within threshold — proceeding with IED removal\n');

        obj.elec_ch_with_IED = tableChanSelection.indChansDeselected;
        obj.elec_ch_with_IED = intersect(obj.elec_ch_clean,obj.elec_ch_with_IED); % don't mark already noisy electrodes
        obj.define_clean_channels();

        fprintf(1,'Electrodes with significant IEDs: ');
        fprintf(1,'%d ', obj.elec_ch_with_IED(:)); fprintf('\n');

    end
         
end