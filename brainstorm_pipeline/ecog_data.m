classdef ecog_data < dynamicprops
% ECOG_DATA  Monolithic ECoG/SEEG data container and preprocessing pipeline.
%
% This class is the Brainstorm-pipeline counterpart of ecog_data_v2 in
% MGH_utils/@ecog_data_v2/. The two classes share an identical constructor
% signature and property set, so crunched .mat files produced by either
% class are interchangeable for downstream analysis.
%
% Key improvements over ecog_data_v2:
%   - extract_shanks() operates only on clean channels (avoids label-parse
%     failures on excluded channels)
%   - reference_signal() implements full bipolar referencing inline
%   - combine_data_files() uses explicit sample offsets per stitch segment
%   - Additional preprocessing order: 'defaultSEEGorBOTH' with CAR before
%     bipolar referencing (matches SEEG clinical convention)
%
% Usage: crunched .mat files written by brainstorm_to_mit_crunched_new.m
% contain a variable named 'obj' of this class type.

properties
    %% ---- DATA ----
    elec_data
    bip_elec_data
    stitch_index
    sample_freq
    for_preproc             % preproc
    trial_data              % trial

    %% ---- INFO ----
    subject
    experiment
    trial_timing
    events_table
    condition
    session
    crunched_file_name      % output
    crunched_file_path
    raw_file_name           % input
    raw_file_path

    %% ---- LABELS ----
    elec_ch                 % unipolar
    elec_ch_label 
    elec_ch_prelim_deselect
    elec_ch_with_IED
    elec_ch_with_noise
    elec_ch_user_deselect
    elec_ch_clean
    elec_ch_valid
    elec_ch_type
    bip_ch                  % bipolar
    bip_ch_label
    bip_ch_valid
    bip_ch_grp             
    bip_ch_label_grp

    %% ---- ANATOMY ----
    anatomy

    %% ---- STATS ----
    % Populated by extract_normalization_metrics(), extract_time_significance(),
    % and extract_significant_channel().
    stats
    
end


methods
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ecog_data(...
            for_preproc,...             % DATA
            subject,...                 % INFO
            experiment,...
            crunched_file_name,...
            crunched_file_path,...
            raw_file_name,...
            raw_file_path,...
            elec_ch_label,...           % LABELS
            elec_ch,...
            elec_ch_prelim_deselect,...
            elec_ch_type)

        obj.for_preproc=for_preproc;

        obj.subject=subject;
        obj.experiment=experiment;
        obj.crunched_file_name=crunched_file_name;
        obj.crunched_file_path=crunched_file_path;
        obj.raw_file_name=raw_file_name;
        obj.raw_file_path=raw_file_path;

        obj.elec_ch=elec_ch;
        obj.elec_ch_label=elec_ch_label;
        obj.elec_ch_prelim_deselect=elec_ch_prelim_deselect;
        obj.elec_ch_type=elec_ch_type;
        
    end


    %% PREPROCESSING PIPELINE


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESS SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function preprocess_signal(obj,varargin)
        p = inputParser();
        addParameter(p,'order','defaultBOTH');
        addParameter(p,'isPlotVisible',true);
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        if strcmp(ops.order,'defaultECOG')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',... 
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'defaultSEEGorBOTH')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',... 
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...       
            };

        elseif strcmp(ops.order,'SEEGorBOTHbyShank')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',...
                     'ShankCSR',... 
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...       
            };

        elseif strcmp(ops.order,'defaultMCJandBJH') 
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'BipolarReferencing'...
                     'GaussianFilterExtraction',...
                     'removeOutliers',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'preEnvelopeExtractionECOG')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'CAR',...
                     'downsample'...
            };

        elseif strcmp(ops.order,'preEnvelopeExtractionSEEGorBOTH')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'GlobalMeanRemoval',...
                     'BipolarReferencing'...
                     'downsample'...
            };
    
        elseif strcmp(ops.order,'preEnvelopeExtractionMCJandBJH') 
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',...
                     'IEDRemoval',...
                     'visualInspection',...
                     'BipolarReferencing'...      
                     'downsample'...
            };

        elseif strcmp(ops.order,'test')
            order = {'downsample'};

        else
            error('Preprocessing order not recognized: %s', ops.order);
        end

        obj.for_preproc.order = order;
        obj.for_preproc.isPlotVisible = ops.isPlotVisible; 

        names = {'highpassFilter',...
                 'notchFilter',...
                 'IEDRemoval',...
                 'visualInspection',...
                 'GlobalMeanRemoval',...
                 'CAR',...
                 'ShankCSR',...
                 'BipolarReferencing',...
                 'GaussianFilterExtraction',...
                 'BandpassExtraction',...
                 'zscore',...
                 'downsample',...
                 'removeOutliers'...
        };
        functions = {'highpass_filter',...
                     'notch_filter',...
                     'remove_IED',...
                     'visual_inspection',...
                     'reference_signal',...
                     'reference_signal',...
                     'reference_signal',...
                     'reference_signal',...
                     'extract_high_gamma',...
                     'extract_high_gamma',...
                     'zscore_signal',...
                     'downsample_signal',...
                     'remove_outliers'...
        };

        name_to_function = containers.Map(names,functions);
        function_order = cellfun(@(x) name_to_function(x),order,'UniformOutput',false);

        fprintf(1,'\nSTARTING TO PREPROCESS SIGNAL\n');

        obj.first_step('doneVisualInspection',ops.doneVisualInspection);

        prev_step = '';
        for i=1:length(function_order)
            step = function_order{i};

            if strcmp(step,prev_step)
                continue
            end

            if strcmp(step,'reference_signal')
                j = i; flags = ""; 
                while  1
                    flags = strcat(flags,"'do",order{j},"',true,");
                    j = j+1;
                    try
                        next_step = function_order{j};
                    catch
                        break;
                    end
                    if ~strcmp(step,next_step)
                        break;
                    end
                end
                flags = char(flags); flags = flags(1:end-1);
                eval(strcat("obj.",step,"(",flags,");"));

            elseif strcmp(step,'extract_high_gamma')
                flags = strcat("'do",order{i},"',true");
                eval(strcat("obj.",step,"(",flags,");"));

            elseif strcmp(step,'visual_inspection')
                eval(strcat("obj.",step,"('doneVisualInspection',",num2str(ops.doneVisualInspection),");"));
                
            else
                eval(strcat("obj.",step,"();"));
            end

            prev_step = step;

        end

        fprintf(1,'\nDONE PREPROCESSING SIGNAL \n');

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HIGHPASS FILTER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function highpass_filter(obj)
        signal = obj.elec_data';
        highpass = obj.for_preproc.highpass;

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Highpass filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);

            for idx_channel=1:size(signal_,2)
                warning('off', 'signal:filtfilt:ParseSOS');
                signal_(:,idx_channel) = filtfilt(highpass.sos,highpass.g,double(signal_(:,idx_channel)));
                fprintf(1,'.');
            end

            signal(obj.stitch_index(k):stop,:) = signal_;
            fprintf(1,'] done\n');
        end

        obj.elec_data = signal';
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTCH FILTER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function notch_filter(obj)
        signal = obj.elec_data';
        notch = obj.for_preproc.notch;

        signal_noise_before = obj.measure_line_noise(signal);

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Notch filtering signal from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);

            for idx_channel=1:size(signal_,2)
                signal_preliminary = double(signal_(:,idx_channel));
                for idx = 1:length(obj.for_preproc.filter_params.notch.fcenter) %#ok<PFBNS>
                    signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
                end 
                signal_(:,idx_channel) = signal_preliminary;
                fprintf(1,'.');
            end

            signal(obj.stitch_index(k):stop,:) = signal_;
            fprintf(1,'] done\n');
        end

        obj.elec_data = signal';

        signal_noise_after = obj.measure_line_noise(signal);

        obj.elec_ch_with_noise = obj.elec_ch(signal_noise_after(:,2) > (mean(signal_noise_after(:,2))+5*std(signal_noise_after(:,2))));
        obj.elec_ch_with_noise = intersect(obj.elec_ch_clean,obj.elec_ch_with_noise);
        obj.define_clean_channels();

        obj.for_preproc.notchFilter_results.signal_noise_before_notch = signal_noise_before;
        obj.for_preproc.notchFilter_results.mean_signal_noise_before_notch = mean(signal_noise_before(obj.elec_ch_clean,2));
        obj.for_preproc.notchFilter_results.signal_noise_after_notch = signal_noise_after;
        obj.for_preproc.notchFilter_results.mean_signal_noise_after_notch = mean(signal_noise_after(obj.elec_ch_clean,2));

        fprintf(1,'\nReduced 50 Hz noise from %.2f to %.2f uV\n',mean(signal_noise_before(obj.elec_ch_clean,2)),mean(signal_noise_after(obj.elec_ch_clean,2)));
        fprintf(1,'Electrodes with significant line noise: ');
        fprintf(1,'%d ', obj.elec_ch_with_noise(:)); fprintf('\n');

        if obj.for_preproc.isPlotVisible
            f = obj.plot_line_noise(signal_noise_before,signal_noise_after);
            prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove due to significant line noise? (format: [1,2]) - ';
            more_line_noise = input(prompt1)';
            obj.elec_ch_with_noise = union(obj.elec_ch_with_noise,more_line_noise);
            obj.define_clean_channels();
            close(f);
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IED REMOVAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function remove_IED(obj)
        signal = double(obj.elec_data');

        fprintf(1, '\n> Finding electrodes with significant IEDs ... \n');
        
        detectionIEDs = [];
        detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1';
        detectionIEDs.segments = [];
            
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
            
        detectionIEDs.tableChanSelection = [];  
        detectionIEDs.threshold = 6.5;
            
        currIEDs.fs = 200;
        currIEDs.discharges.MV = [];
        currIEDs.numSamples = 0;
        currIEDs.discharges.MV = detectionIEDs.discharges.MV;
        currIEDs.numSamples= currIEDs.numSamples + size(detectionIEDs.d_decim, 1);
            
        numSpikes     = sum(currIEDs.discharges.MV==1, 1);
        totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60;
        numSpikes_min = numSpikes / totalDuration;
        numSpikes     = transpose(numSpikes);
        numSpikes_min = transpose(numSpikes_min);
            
        indChanSelected = find(numSpikes_min < detectionIEDs.threshold);
        tableChanSelection.numSpikesAll           = numSpikes_min;
        tableChanSelection.indChansSelected       = indChanSelected;
        tableChanSelection.indChansDeselected     = setdiff(obj.elec_ch,indChanSelected);
        tableChanSelection.nameChansSelected      = transpose(obj.elec_ch_label(indChanSelected));
        tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);
            
        obj.for_preproc.IEDRemoval_results=tableChanSelection;
        obj.for_preproc.IEDRemoval_results.threshold=detectionIEDs.threshold;
  
        if length(tableChanSelection.indChansDeselected) > ceil(size(obj.elec_ch,1)/3)
            fprintf(1,'Too many electrodes with significant IEDs, SKIPPING STEP\n')
            new_order_mask = cell2mat(cellfun(@(x) strcmp(x,'IEDRemoval'),obj.for_preproc.order,'UniformOutput',false));
            obj.for_preproc.order = obj.for_preproc.order(~new_order_mask);
        else 
            obj.elec_ch_with_IED = tableChanSelection.indChansDeselected;
            obj.elec_ch_with_IED = intersect(obj.elec_ch_clean,obj.elec_ch_with_IED);
            obj.define_clean_channels();
            fprintf(1,'Electrodes with significant IEDs: ');
            fprintf(1,'%d ', obj.elec_ch_with_IED(:)); fprintf('\n');
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUAL INSPECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function visual_inspection(obj,varargin)
        p = inputParser();
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        signal = obj.elec_data';

        fprintf(1,'\n> Visually inspecting signal ...\n');

        if obj.for_preproc.isPlotVisible
            obj.plot_channels(signal,...
                            obj.elec_ch_label,...
                            obj.elec_ch_clean,...
                            obj.elec_ch_valid,...
                            'stitch_index',obj.stitch_index,...
                            'sample_freq',obj.sample_freq,...
                            'downsample',true,...
                            'decimation_freq',obj.for_preproc.decimation_freq...
            );
        end

        if ~ops.doneVisualInspection
            prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove from visual inspection? (format: [1,2]) - ';
            prompt2 = 'Your name please :) - ';
            obj.elec_ch_user_deselect = input(prompt1)';
            vi_ops.inspected = 1;
            vi_ops.inspected_by = input(prompt2,'s');
            vi_ops.inspection_date = datestr(now, 'yyyy/mm/dd-HH:MM');
            obj.for_preproc.visualInspection_results = vi_ops;
        end
        
        obj.define_clean_channels();
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REFERENCE SIGNAL (CAR / GlobalMean / ShankCSR / Bipolar)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reference_signal(obj,varargin)
        p = inputParser();
        addParameter(p,'doGlobalMeanRemoval',false)
        addParameter(p,'doCAR',false);
        addParameter(p,'doShankCSR',false);
        addParameter(p,'doBipolarReferencing',false);
        parse(p, varargin{:});
        ops = p.Results;

        signal = obj.elec_data';
        signal_for_bip = obj.elec_data';

        ecog_chans = find(strcmp(obj.elec_ch_type, 'ecog_grid') | strcmp(obj.elec_ch_type, 'ecog_strip'));
        seeg_chans = find(strcmp(obj.elec_ch_type, 'seeg'));

        signal_bipolar = zeros(size(signal));

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Referencing signal from file %d of %d ... \n',k,length(obj.stitch_index));
            
            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            % GLOBAL MEAN REMOVAL
            if ops.doGlobalMeanRemoval
                fprintf(1, '\n>> Removing global mean of signal \n');
                signal_ = signal(obj.stitch_index(k):stop,:);
                overall_mean = mean(signal_(:,obj.elec_ch_clean),2);
                signal_ = signal_ - repmat(overall_mean,1,size(signal_,2));
                signal(obj.stitch_index(k):stop,:) = signal_;
            end

            % COMMON AVERAGE REFERENCING
            if ops.doCAR && ~isempty(ecog_chans)
                fprintf(1, '\n>> Common average filtering signal \n');
                fprintf(1,'[');
                signal_ = signal(obj.stitch_index(k):stop,:);
                num_amps = ceil(size(signal_,2) / obj.for_preproc.elecs_per_amp);
                eligible_channels_ecog = intersect(ecog_chans,obj.elec_ch_clean);
                for idx_amp = 1:num_amps
                    idx_low  = (idx_amp-1)*obj.for_preproc.elecs_per_amp+1;
                    idx_high = min((idx_amp-0)*obj.for_preproc.elecs_per_amp+0, max(eligible_channels_ecog));
                    list_channels = intersect(eligible_channels_ecog,idx_low:idx_high);
                    if ~isempty(list_channels) && length(list_channels)>1
                        signal_mean = mean(signal_(:,list_channels),2);
                        for idx_ch = list_channels
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    end
                end
                fprintf(1,'] done\n');
                signal(obj.stitch_index(k):stop,:) = signal_;
            end

            % SHANK COMMON SOURCE REMOVAL
            if ops.doShankCSR && ~isempty(seeg_chans)
                fprintf(1, '\n>> Shank common source removal \n');
                fprintf(1,'[');
                signal_ = signal(obj.stitch_index(k):stop,:);
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);
                [shank_locs,~,~] = obj.extract_shanks();
                for idx_shk = 1:size(shank_locs,1)
                    same_shank = intersect(eligible_channels_seeg, shank_locs{idx_shk});
                    if ~isempty(same_shank) && length(same_shank)>1
                        signal_mean = mean(signal_(:,same_shank),1);
                        for idx_ch = same_shank
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    end
                end
                signal(obj.stitch_index(k):stop,:) = signal_;
                fprintf(1,'] done\n');
            elseif ops.doShankCSR 
                error('No SEEG channels to perform shank CSR on')
            end

            % BIPOLAR REFERENCING
            if ops.doBipolarReferencing && ~isempty(seeg_chans)
                fprintf(1, '\n>> Bipolar referencing signal \n');
                fprintf(1,'[');
                signal_ = signal_for_bip(obj.stitch_index(k):stop,:);
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);

                [shank_locs,~,~] = obj.extract_shanks();

                chan_idx_for_bip = cellfun(@(x) intersect(x,eligible_channels_seeg,'stable'),shank_locs,'uni',false);
                
                chan_num_for_bip = cell(length(chan_idx_for_bip), 1);
                for iShank = 1:length(chan_idx_for_bip)
                    original_indices = chan_idx_for_bip{iShank};
                    electrode_numbers = zeros(length(original_indices), 1);
                    for j = 1:length(original_indices)
                        label = obj.elec_ch_label{original_indices(j)};
                        num_str = regexp(label, '\d+', 'match');
                        if ~isempty(num_str)
                            electrode_numbers(j) = str2double(num_str{1});
                        else
                            electrode_numbers(j) = NaN;
                        end
                    end
                    chan_num_for_bip{iShank} = electrode_numbers;
                end
                
                bipolar_diffs_idx_grp = cellfun(@(x,y) [y(find(diff(x)==1 & diff(y)==1))+1, y(find(diff(x)==1 & diff(y)==1))], chan_num_for_bip, chan_idx_for_bip, 'uni', false)';
                bipolar_diffs_name_grp = cellfun(@(x) obj.elec_ch_label(x), bipolar_diffs_idx_grp, 'uni', false);
                
                bipolar_diffs_idx = cell2mat(bipolar_diffs_idx_grp');
                bipolar_diffs_name = obj.elec_ch_label(bipolar_diffs_idx);
                
                bipolar_idxs = 1:size(bipolar_diffs_idx,1);
                bipolar_valid = ones(size(bipolar_diffs_idx,1),1);
            
                signal_bipolar_= double([]); 
                for bipolar_id = 1:size(bipolar_diffs_idx,1)
                    bipol_ch_1 = bipolar_diffs_idx(bipolar_id,1);
                    bipol_ch_2 = bipolar_diffs_idx(bipolar_id,2);
                    bip_ch_name1 = bipolar_diffs_name{bipolar_id,1};
                    bip_ch_name2 = bipolar_diffs_name{bipolar_id,2};
              
                    B = cell2mat(extract(bip_ch_name2,lettersPattern));
                    A = cell2mat(extract(bip_ch_name1,lettersPattern));
                    assert(all(A==B),"Some channels are not on the same shank");
                    B = str2num(cell2mat(extract(bip_ch_name2,digitsPattern))); %#ok<ST2NM>
                    A = str2num(cell2mat(extract(bip_ch_name1,digitsPattern))); %#ok<ST2NM>
                    assert(A-B==1, "Some channels are more than one apart");
                    
                    signal_bipolar_(:,bipolar_id) = signal_(:,bipol_ch_1)-signal_(:,bipol_ch_2);
                    fprintf(1,'.');
                end
            
                bipolar_diffs_name = arrayfun(@(x) {[bipolar_diffs_name{x,1} '-' bipolar_diffs_name{x,2}]}, ...
                    [1:size(bipolar_diffs_name,1)])';
            
                signal_bipolar(obj.stitch_index(k):stop,1:size(bipolar_diffs_idx,1)) = signal_bipolar_;
            
                fprintf(1,'] done\n');
            
            elseif ops.doBipolarReferencing
                error('No SEEG channels to perform bipolar referencing on')
            end
        end

        if ops.doGlobalMeanRemoval || ops.doCAR || ops.doShankCSR
            obj.elec_data = signal';
        end

        if ops.doBipolarReferencing 
            signal_bipolar = signal_bipolar(:,1:size(bipolar_diffs_idx,1));
            obj.bip_elec_data    = signal_bipolar';
            obj.bip_ch           = bipolar_idxs';
            obj.bip_ch_label     = bipolar_diffs_name; 
            obj.bip_ch_valid     = bipolar_valid;
            obj.bip_ch_grp       = bipolar_diffs_idx_grp';
            obj.bip_ch_label_grp = bipolar_diffs_name_grp';
              
            if obj.for_preproc.isPlotVisible
                obj.plot_channels(signal_bipolar,...
                                obj.bip_ch_label,...
                                obj.bip_ch,...
                                obj.bip_ch_valid,...
                                'stitch_index',obj.stitch_index,...
                                'sample_freq',obj.sample_freq,...
                                'downsample',true,...
                                'decimation_freq',obj.for_preproc.decimation_freq...
                );
            end
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT HIGH GAMMA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_high_gamma(obj,varargin)
        p = inputParser();
        addParameter(p,'doGaussianFilterExtraction',false);
        addParameter(p,'doBandpassExtraction',false);
        % NAPLAB/Columbia filterbank — alternative to Chang-lab Gaussian method.
        % Calls CUprocessingHilbertTransform_filterbankGUI (embedded below).
        addParameter(p,'doNapLabFilterExtraction',false);
        parse(p, varargin{:});
        ops = p.Results;

        nMethods = ops.doGaussianFilterExtraction + ops.doBandpassExtraction + ops.doNapLabFilterExtraction;
        if nMethods ~= 1
            error('Must specify exactly one extraction method (doGaussianFilterExtraction, doBandpassExtraction, or doNapLabFilterExtraction).')
        end

        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        if ~isfield(obj.for_preproc,'gaussian') || ~isfield(obj.for_preproc,'bandpass')
            obj.define_parameters()
        end

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Extracting high gamma signal from file %d of %d ... \n',k,length(obj.stitch_index));

            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);

            if ops.doGaussianFilterExtraction
                cfs = obj.for_preproc.gaussian.cfs;
                sds = obj.for_preproc.gaussian.sds;

                filter_bank={};
                for s=1:length(cfs)
                    filter_bank{s} = obj.gaussian_filter(transpose(signal_(:,1)),obj.sample_freq,cfs(s),sds(s));
                end
                obj.for_preproc.gaussian.filter_banks = filter_bank;

                fprintf(1, '\n>> Extracting unipolar high gamma envelope based on gaussian filtering \n');
                fprintf(1,'[');
                signal_hilbert = nan*signal_;
                for kk=1:size(signal_,2)
                    signal_hilbert_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                    signal_hilbert(:,kk) = transpose(mean(signal_hilbert_all,1));
                    fprintf(1,'.');
                end 
                fprintf(1,'] done\n');
                signal(obj.stitch_index(k):stop,:) = signal_hilbert;

                if ~isempty(obj.bip_elec_data)
                    fprintf(1, '\n>> Extracting bipolar high gamma envelope based on gaussian filtering \n');
                    fprintf(1,'[');
                    signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                    signal_hilbert_bipolar = nan*signal_bipolar_;
                    for kk=1:size(signal_bipolar_,2)
                        signal_hilbert_bipolar_all = cell2mat(cellfun(@abs,obj.hilbert_transform(double(transpose(signal_bipolar_(:,kk))),obj.sample_freq,filter_bank),'UniformOutput',false));
                        signal_hilbert_bipolar(:,kk) = transpose(mean(signal_hilbert_bipolar_all,1));
                        fprintf(1,'.');
                    end
                    fprintf(1,'] done\n')
                    signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                end

            elseif ops.doBandpassExtraction
                B = obj.for_preproc.bandpass.B;
                A = obj.for_preproc.bandpass.A;
                signal_hilbert = filtfilt(B,A,double(signal_));
                signal_hilbert = abs(hilbert(signal_hilbert));
                signal_hilbert(signal_hilbert < 0) = 0;

                if ~isempty(obj.bip_elec_data)
                    signal_hilbert_bipolar = filtfilt(B,A,double(signal_bipolar_));
                    signal_hilbert_bipolar = abs(hilbert(signal_hilbert_bipolar));
                    signal_hilbert_bipolar(signal_hilbert_bipolar < 0) = 0;
                    signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                end

            % ----------------------------------------------------------
            % NAPLAB / COLUMBIA UNIVERSITY FILTERBANK
            % Log-spaced Gaussian filterbank from the Neural Acoustic
            % Processing Lab, Columbia University (naplab.ee.columbia.edu).
            % Identical in design philosophy to the Chang-lab method but
            % with octave spacing of 1/7 and slightly different sigma.
            % ----------------------------------------------------------
            elseif ops.doNapLabFilterExtraction
                freqRange = [obj.for_preproc.filter_params.bandpass.f_gamma_low, ...
                             obj.for_preproc.filter_params.bandpass.f_gamma_high];

                fprintf(1, '\n>> Extracting unipolar high gamma (NAPLAB filterbank) \n');
                [dh,~,~] = ecog_data.naplab_filterbank(signal_', obj.sample_freq, freqRange);
                signal_hilbert = mean(abs(dh), 3)';
                signal(obj.stitch_index(k):stop,:) = signal_hilbert;

                if ~isempty(obj.bip_elec_data)
                    fprintf(1, '\n>> Extracting bipolar high gamma (NAPLAB filterbank) \n');
                    signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                    [dh_bip,~,~] = ecog_data.naplab_filterbank(signal_bipolar_', obj.sample_freq, freqRange);
                    signal_hilbert_bipolar = mean(abs(dh_bip), 3)';
                    signal_bipolar(obj.stitch_index(k):stop,:) = signal_hilbert_bipolar;
                end
            end
        end

        obj.elec_data = signal';
        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ZSCORE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function zscore_signal(obj)
        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        for k=1:length(obj.stitch_index)
            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end
            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_ = zscore(signal_);
            signal(obj.stitch_index(k):stop,:) = signal_;
        
            if ~isempty(obj.bip_elec_data)
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_bipolar_ = zscore(signal_bipolar_);
                signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
            end
        end
        
        obj.elec_data = signal';
        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DOWNSAMPLE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function downsample_signal(obj,varargin)
        p = inputParser();
        addParameter(p,'decimationFreq',obj.for_preproc.decimation_freq)
        parse(p, varargin{:});
        ops = p.Results;
        
        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        signal_dec = [];
        signal_bipolar_dec = [];
        curr_stitch = 1;
        stitch_index = [];

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Resampling signal from file %d of %d ... \n',k,length(obj.stitch_index));

            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_ = resample(double(signal_),ops.decimationFreq,obj.sample_freq);
            signal_dec = [signal_dec; signal_];
        
            if ~isempty(obj.bip_elec_data)
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_bipolar_ = resample(double(signal_bipolar_),ops.decimationFreq,obj.sample_freq);
                signal_bipolar_dec = [signal_bipolar_dec; signal_bipolar_];
            end

            stitch_index = [stitch_index; curr_stitch];
            curr_stitch = curr_stitch + size(signal_,1);
        end

        obj.elec_data = signal_dec';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar_dec';
        end

        if ~isempty(obj.trial_timing) && (ops.decimationFreq ~= obj.for_preproc.decimation_freq)
            decimation_factor = obj.sample_freq / ops.decimationFreq;
            trial_timing = cell(size(obj.trial_timing));
            for i=1:size(trial_timing,1)
                tmp_table = obj.trial_timing{i,1};
                tmp_table.start = round(tmp_table.start / decimation_factor);
                tmp_table.end = round(tmp_table.end / decimation_factor);
                trial_timing{i,1} = tmp_table;
            end
            obj.trial_timing = trial_timing;
        end

        if ~isempty(obj.trial_timing) && (ops.decimationFreq == obj.for_preproc.decimation_freq)
            obj.trial_timing = obj.for_preproc.trial_timing_dec;
        end

        obj.stitch_index = stitch_index;
        obj.sample_freq = round(ops.decimationFreq);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE OUTLIERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function remove_outliers(obj)
        ops.trimmed      = obj.for_preproc.outlier.trimmed;
        ops.threshold    = obj.for_preproc.outlier.threshold;
        ops.percentile   = obj.for_preproc.outlier.percentile;
        ops.buffer       = obj.for_preproc.outlier.buffer;
        ops.interpMethod = obj.for_preproc.outlier.interpMethod;

        signal = obj.elec_data';
        signal = obj.zero_out_signal(signal);

        outlierRemoval_results.idxs    = [];
        outlierRemoval_results.prcnts  = [];
        outlierRemoval_results.ignored = [];

        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
            signal_bipolar = obj.zero_out_signal(signal_bipolar);
            outlierRemoval_results.idxs_bipolar    = [];
            outlierRemoval_results.prcnts_bipolar  = [];
            outlierRemoval_results.ignored_bipolar = [];
        end

        for k=1:length(obj.stitch_index)
            fprintf(1, '\n> Removing outliers from file %d of %d ... \n',k,length(obj.stitch_index));
            fprintf(1,'[');

            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end

            signal_ = signal(obj.stitch_index(k):stop,:);
            [signal_,idxs,prcnts,ignored] = obj.envelope_outliers(signal_,ops);
            signal(obj.stitch_index(k):stop,:) = signal_;
            fprintf(1,'] done\n');

            outlierRemoval_results.idxs = [outlierRemoval_results.idxs, idxs];
            outlierRemoval_results.prcnts = [outlierRemoval_results.prcnts, prcnts];
            outlierRemoval_results.ignored = [outlierRemoval_results.ignored; ismember(obj.elec_ch,ignored)'];
            
            ignored_not_noisy = intersect(obj.elec_ch_clean,ignored);
            fprintf(1,'Ignored unipolar channels: ');
            fprintf(1,'%d ',ignored_not_noisy(:)); fprintf('\n');

            if ~isempty(obj.bip_elec_data)
                fprintf(1,'[');
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                [signal_bipolar_,idxs_bipolar,prcnts_bipolar,ignored_bipolar] = obj.envelope_outliers(signal_bipolar_,ops);
                signal_bipolar(obj.stitch_index(k):stop,:) = signal_bipolar_;
                fprintf(1,'] done\n');

                outlierRemoval_results.idxs_bipolar = [outlierRemoval_results.idxs_bipolar, idxs_bipolar];
                outlierRemoval_results.prcnts_bipolar = [outlierRemoval_results.prcnts_bipolar, prcnts_bipolar];
                outlierRemoval_results.ignored_bipolar = [outlierRemoval_results.ignored_bipolar, ismember(obj.bip_ch,ignored_bipolar)'];

                ignored_bipolar_not_noisy = intersect(obj.bip_ch,ignored_bipolar);
                fprintf(1,'Ignored bipolar channels: ');
                fprintf(1,'%d ',ignored_bipolar_not_noisy(:)); fprintf('\n');
            end
        end

        obj.elec_data = signal';
        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar';
        end

        obj.for_preproc.outlierRemoval_results = outlierRemoval_results;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADD ANATOMY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function add_anatomy(obj,anatomy_path,varargin)
        p = inputParser();
        addParameter(p,'veraFolder',true);
        addParameter(p,'subdirName','VERA_');
        addParameter(p,'templateName','cvs_avg35_inMNI152');
        parse(p, varargin{:});
        ops = p.Results;

        if ops.veraFolder
            folder = [ops.subdirName obj.subject filesep];
            template_folder = [ops.templateName filesep];
        else
            folder = '';
            template_folder = '';
        end

        filename = [anatomy_path folder obj.subject '_brain.mat'];
        subject_space = load(filename);
        obj.anatomy.subject_space = subject_space;

        filename = [anatomy_path folder obj.subject '_MNI_brain.mat'];
        mni_space = load(filename);
        obj.anatomy.mni_space = mni_space;

        filename = [anatomy_path template_folder ops.templateName '.mat'];
        template_brain = load(filename);
        obj.anatomy.template_brain = template_brain;

        [mapping,labels] = obj.channel_mapping_anatomical(subject_space);
        obj.anatomy.mapping = mapping;
        obj.anatomy.labels = labels;
        
        hemisphere = cell(size(subject_space.tala.electrodes,1),1);
        right_idxs = find(subject_space.tala.electrodes(:,1) > 0);
        hemisphere(right_idxs,1) = {'right'};
        left_idxs = find(subject_space.tala.electrodes(:,1) < 0);
        hemisphere(left_idxs,1) = {'left'};
        obj.anatomy.hemisphere = hemisphere;
    end


    %% HELPER METHODS


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIRST STEP (reset to raw state)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = first_step(obj,varargin)
        p = inputParser();
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        obj.stitch_index = obj.for_preproc.stitch_index_raw;
        obj.sample_freq  = round(obj.for_preproc.sample_freq_raw);
        obj.define_parameters();

        obj.elec_ch_with_IED = [];
        obj.elec_ch_with_noise = [];
        if ~ops.doneVisualInspection
            obj.elec_ch_user_deselect = [];
        end

        obj.define_clean_channels()

        obj.for_preproc.notchFilter_results = [];
        obj.for_preproc.IEDRemoval_results = [];
        obj.for_preproc.visualInspection_results = [];
        obj.for_preproc.outlierRemoval_results = [];

        obj.bip_elec_data    = [];
        obj.bip_ch           = [];
        obj.bip_ch_label     = [];
        obj.bip_ch_valid     = [];
        obj.bip_ch_grp       = [];
        obj.bip_ch_label_grp = [];

        if isfield(obj.for_preproc,'trial_timing_raw')
            obj.trial_timing = obj.for_preproc.trial_timing_raw;
        end

        obj.elec_data = obj.for_preproc.elec_data_raw;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_parameters(obj)
        param = struct;
            
        param.highpass.Wp = 0.50;
        param.highpass.Ws = 0.05;
        param.highpass.Rp = 3;
        param.highpass.Rs = 30;
            
        param.peak.fcenter = [45,50,55];
        param.peak.bw      = ones(1,length(param.peak.fcenter)).*0.001;
            
        param.notch.fcenter = [50:50:250];
        param.notch.bw = ones(1,length(param.notch.fcenter)).*0.001;

        param.gaussian.f_gamma_low = 70;
        param.gaussian.f_gamma_high = 150;

        param.bandpass.f_gamma_low = 70;
        param.bandpass.f_gamma_high = 150;
        param.bandpass.filter_order = 6;

        highpass.Wp = param.highpass.Wp/(obj.sample_freq/2); 
        highpass.Ws = param.highpass.Ws/(obj.sample_freq/2);
        highpass.Rp = param.highpass.Rp; 
        highpass.Rs = param.highpass.Rs;
        [highpass.n,highpass.Wn] = buttord(highpass.Wp,highpass.Ws,highpass.Rp,highpass.Rs);
        highpass.n = highpass.n + rem(highpass.n,2);
        [highpass.z,highpass.p,highpass.k] = butter(highpass.n,highpass.Wn,'high');
        [highpass.sos,highpass.g] = zp2sos(highpass.z,highpass.p,highpass.k);
        highpass.h = dfilt.df2sos(highpass.sos,highpass.g);

        for idx = 1:length(param.peak.fcenter)
            peak{idx}.wo = param.peak.fcenter(idx)/(obj.sample_freq/2);
            peak{idx}.bw = param.peak.bw(idx);  
            [peak{idx}.b,peak{idx}.a] = iirpeak(peak{idx}.wo,peak{idx}.bw);
        end

        for idx = 1:length(param.notch.fcenter)
            notch{idx}.wo = param.notch.fcenter(idx)/(obj.sample_freq/2);  
            notch{idx}.bw = param.notch.bw(idx);
            [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
        end 

        [gaussian.cfs,gaussian.sds] = obj.get_filter_param_chang_lab(param.gaussian.f_gamma_low,param.gaussian.f_gamma_high);
        gaussian.cfs(1) = 73.0;
        gaussian.f_gamma_low = param.gaussian.f_gamma_low;
        gaussian.f_gamma_high = param.gaussian.f_gamma_high;

        bandpass.h = fdesign.bandpass('N,F3dB1,F3dB2',param.bandpass.filter_order,param.bandpass.f_gamma_low,param.bandpass.f_gamma_high,obj.sample_freq);
        bandpass.Hd = design(bandpass.h,'butter');
        [bandpass.B, bandpass.A] = sos2tf(bandpass.Hd.sosMatrix,bandpass.Hd.scaleValues);
        bandpass.filter_order = param.bandpass.filter_order;
        bandpass.f_gamma_low = param.bandpass.f_gamma_low;
        bandpass.f_gamma_high = param.bandpass.f_gamma_high;

        zero = 1;
        outlier.trimmed = 1;
        outlier.threshold = 5;
        outlier.percentile = 0.9;
        outlier.buffer = 20;
        outlier.interpMethod = 'linear';

        obj.for_preproc.filter_params   = param;
        obj.for_preproc.highpass        = highpass;
        obj.for_preproc.peak            = peak;
        obj.for_preproc.notch           = notch;
        obj.for_preproc.gaussian        = gaussian;
        obj.for_preproc.bandpass        = bandpass;
        obj.for_preproc.zero            = zero;
        obj.for_preproc.outlier         = outlier;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ZERO OUT SIGNAL EDGES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function signal=zero_out_signal(obj,signal)
        samples_to_remove = obj.sample_freq*obj.for_preproc.zero;
        for k=1:length(obj.stitch_index)
            if k == length(obj.stitch_index)
                stop = size(signal,1);
            else
                stop = obj.stitch_index(k+1)-1;
            end
            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_(1:samples_to_remove,:) = zeros(samples_to_remove,size(signal_,2));
            signal_((end-samples_to_remove+1):end,:) = zeros(samples_to_remove,size(signal_,2));
            signal(obj.stitch_index(k):stop,:) = signal_;
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASURE LINE NOISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function signal_noise=measure_line_noise(obj,signal)
        fprintf(1, '\n> Measuring 50Hz noise power ...\n');
        fprintf(1,'[');
        peak = obj.for_preproc.peak;
        signal_noise = zeros(size(signal,2),length(peak));
        for idx_channel=1:size(signal,2)
            for idx_filter=1:length(peak)
                signal_noise(idx_channel,idx_filter) = mean(abs(filter(peak{idx_filter}.b,peak{idx_filter}.a,signal(:,idx_channel))));
            end
            fprintf(1,'.');
        end
        fprintf(1,'] done\n');
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT LINE NOISE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f=plot_line_noise(obj,noise_before,noise_after)
        fprintf(1, '\n> Plotting 50Hz noise power ...\n');
        close all
        f = figure;
        set(gcf,'position',[30,30,2300,900]);
        c= [0.4660 0.6740 0.1880];
        x = find(~obj.elec_ch_valid);
        idxs = ~obj.elec_ch_valid;

        currsub = subplot(2,2,1);       
        stem(noise_before,'filled'); 
        axis tight; hold on;
        stem(x,noise_before(idxs,:),'filled','Color','k')
        legend({'45Hz noise','50Hz noise','55Hz noise','MARKED NOISY'},'Location','best','FontSize',16,'Box','off');
        ylabel('Noise (uV)','FontSize',18);
        title('BEFORE NOTCH FILTERING','FontSize',22);
        obj.update_position(currsub);

        currsub = subplot(2,2,3);
        stem(noise_before(:,2)./mean(noise_before(:,[1,3]),2),'filled','Color',c); 
        axis tight; hold on;
        stem(x,noise_before(idxs,2)./mean(noise_before(idxs,[1,3]),2),'filled','Color','k')
        xlabel('Channel #','FontSize',18);
        ylabel('50Hz noise / mean 45Hz+55Hz noise','FontSize',18)
        obj.update_position(currsub);

        currsub = subplot(2,2,2);
        stem(noise_after,'filled'); 
        axis tight; hold on;
        stem(x,noise_after(idxs,:),'filled','Color','k')
        title('AFTER NOTCH FILTERING','FontSize',22)
        obj.update_position(currsub);

        currsub = subplot(2,2,4);
        stem(noise_after(:,2)./mean(noise_after(:,[1,3]),2),'filled','Color',c); 
        axis tight; hold on;
        stem(x,noise_after(idxs,2)./mean(noise_after(idxs,[1,3]),2),'filled','Color','k')
        xlabel('Channel #','FontSize',18)
        obj.update_position(currsub);

        filename = split(obj.crunched_file_name,'/');
        filename = split(filename{end},'.');
        filename = [filename{1} '_line_noise.png'];
        saveas(gcf,strcat(filename));
        set(0, 'CurrentFigure', f);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE SUBPLOT POSITION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function update_position(obj,currsub)
        pos = get(currsub, 'Position');
        new_pos = pos + [-0.05 -0.05 0.07 0.07];
        set(currsub, 'Position', new_pos)
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE CLEAN CHANNELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_clean_channels(obj)
        bad_elecs = union(obj.elec_ch_prelim_deselect,obj.elec_ch_with_IED);
        bad_elecs = union(bad_elecs,obj.elec_ch_with_noise);
        bad_elecs = union(bad_elecs,obj.elec_ch_user_deselect);
        obj.elec_ch_clean = setdiff(obj.elec_ch,bad_elecs,'stable');
        obj.elec_ch_valid = ismember(obj.elec_ch,obj.elec_ch_clean);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT SHANKS
    % Operates only on clean channels to avoid label-parse failures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [locs,ch_tags,ch_num]=extract_shanks(obj)
        clean_labels = obj.elec_ch_label(obj.elec_ch_clean);
        
        ch_num_cell = extract(clean_labels, digitsPattern);
        ch_tags_cell = extract(clean_labels, lettersPattern);
        
        ch_num = zeros(length(ch_num_cell), 1);
        for i = 1:length(ch_num_cell)
            if ~isempty(ch_num_cell{i})
                ch_num(i) = str2double(ch_num_cell{i});
            else
                ch_num(i) = NaN;
            end
        end
        
        ch_tags = cell(length(ch_tags_cell), 1);
        for i = 1:length(ch_tags_cell)
            if ~isempty(ch_tags_cell{i})
                ch_tags{i} = ch_tags_cell{i};
            else
                ch_tags{i} = '';
            end
        end
        
        valid = ~isnan(ch_num) & ~cellfun(@isempty, ch_tags);
        ch_num = ch_num(valid);
        ch_tags = ch_tags(valid);
        clean_labels = clean_labels(valid);
        
        [tags, ~, ~] = unique(ch_tags, 'stable');
        
        locs = cell(length(tags), 1);
        for i = 1:length(tags)
            tag_mask = strcmp(ch_tags, tags{i});
            channels_with_tag = clean_labels(tag_mask);
            locs{i} = [];
            for j = 1:length(channels_with_tag)
                idx = find(strcmp(obj.elec_ch_label, channels_with_tag{j}));
                if ~isempty(idx)
                    locs{i} = [locs{i}; idx];
                end
            end
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHANG LAB FILTER PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [cfs,sds]=get_filter_param_chang_lab(obj,f_low,f_high)
        cfs_round_factor = 1;
        cfs_round_val = 10^(cfs_round_factor);
        sds_round_factor = 2;
        sds_round_val = 10^(sds_round_factor);
        fq_min = 4.0749286538265;
        fq_max = 200.;
        scale = 7.;
        cfs = 2.^((log2(fq_min) * scale:1: log2(fq_max) * scale) / scale);
        sds = 10.^( log10(.39) + .5 * (log10(cfs)));
        sds = (sds) * sqrt(2.);
        sds=round(sds*sds_round_val)/sds_round_val;
        cfs=round(cfs*cfs_round_val)/cfs_round_val;
        if nargin<2
            return
        end
        index=(cfs<f_high) & (cfs>f_low);
        cfs=cfs(index);
        sds=sds(index);
    end
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GAUSSIAN FILTER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function k=gaussian_filter(obj,X,rate,center,sd)
        N = size(X,2);
        d = 1./rate;
        freq = [0:ceil(N/2-1), ceil(-(N)/2):1:-1]/(d*N);
        k = exp((-(abs(freq) - center).^2)/(2 * (sd^2)));
        k = k/norm(k);
    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HILBERT TRANSFORM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Xh_out,X_fft_h_out]=hilbert_transform(obj,X,rate,varargin)
        switch nargin 
        case 4
            filters=varargin{1};
            phase='None';
            X_fft_h='None';
        case 5
            filters=varargin{1};
            phase=varargin{2};
            X_fft_h='None';
        case 6
            filters=varargin{1};
            phase=varargin{2};
            X_fft_h=varargin{3};
        otherwise 
            filters='None';
            phase='None';
            X_fft_h='None';        
        end 
            
        if ~iscell(filters)
            filters = {filters};
        end 

        time = size(X,2);
        d=1./rate;
        freq=[0:ceil(time/2-1), ceil(-(time)/2):1:-1]/(d*time);
        Xh=cell(length(filters),1);

        if strcmp(X_fft_h, 'None')
            h = zeros(1,length(freq));
            h(freq > 0) = 2.;
            h(1) = 1.;
            X_fft_h = fft(X) .* h;
            if ~strcmp(phase,'None')
                X_fft_h=fft(X).*phase;
            end
        end
            
        for ii=1:size(filters,2)
            f=filters{ii};
            if strcmp(f,'None')
                Xh{ii} = ifft(X_fft_h);
            else
                f = f / max(f);
                Xh{ii} = ifft(X_fft_h .* f);
            end 
        end

        if size(Xh,1) == 1
            Xh_out= Xh{1};
            X_fft_h_out=X_fft_h;
        else 
            Xh_out= Xh;
            X_fft_h_out=X_fft_h;
        end 
    end
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENVELOPE OUTLIER REMOVAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [envelopes_rjct,outliers_idx,outlier_prcnt,ignored_channels]=envelope_outliers(obj,envelopes,ops,varargin)
        n_channels = size(envelopes,2);

        if isvector(ops.threshold)
            ops.threshold = repmat(ops.threshold(:), 1, n_channels);
        else
            assert(size(ops.threshold) == n_channels);
        end

        s = quantile(envelopes, ops.percentile);
        Z = envelopes ./ repmat(s, size(envelopes,1),1);
        clear s;

        outliers = Z > ops.threshold;
        envelopes_rjct = nan*envelopes;
        outlier_prcnt = nan*zeros(size(envelopes,2),1);
        outliers_idx = {};
        ignored_channels = [];

        for q=1:n_channels
            outlie = outliers(:,q);
            outlier_prcnt(q) = sum(outlie)./numel(outlie)*1e2;
            
            envl = envelopes(:,q);
            envl_rejct = envl;
            
            rise = find(diff(outlie)==1).';
            fall = find(diff(outlie)==-1).';

            try 
                assert(numel(rise)==numel(fall));
                process = true;
            catch err 
                ignored_channels = [ignored_channels, q];
                process = false;
            end 

            outlie_idx ={ };
            if numel(rise)>0 && process
                for r=rise
                    f = fall(find(fall>r,1,'first'));
                    interp_win = r:f;
                    interp_val = envl(interp_win).';
                    lookup_window = max([r-ops.buffer,1]):min([(f+ops.buffer),numel(outlie)]);
                    main_sig = envl(lookup_window);
                    [C,ia,ib] = intersect(lookup_window,interp_win);
                    aux = main_sig;
                    fixed_sig = main_sig;
                    aux(ia) = [];
                    intrp_sig = interp1(setdiff(lookup_window,interp_win,'stable'),aux,interp_win,ops.interpMethod);
                    fixed_sig(ia) = intrp_sig;
                    envl_rejct(lookup_window) = fixed_sig;
                    outlie_idx = [outlie_idx,[interp_win;interp_val]];
                end 
            end

            envelopes_rjct(:,q) = envl_rejct;
            outliers_idx = [outliers_idx;{outlie_idx}];
            fprintf(1,'.');
        end 
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMBINE DATA FILES (multi-block sessions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function combine_data_files(obj)
        if isa(obj.trial_timing{1},'table')
            fprintf(1,'No need to combine data files!\n')
            return
        end

        trial_timing_raw = [];
        trial_timing_dec = [];
        condition        = [];
        session          = [];
        events_table     = [];

        fprintf(1, '\n> Combining trial info from data files \n');
        fprintf(1,'[');

        for i=1:length(obj.stitch_index)
            samples_to_add_raw = obj.for_preproc.stitch_index_raw(i)-1;
            samples_to_add_dec = obj.for_preproc.stitch_index_dec(i)-1;
            trial_timing_raw_  = obj.for_preproc.trial_timing_raw{i};
            trial_timing_dec_  = obj.for_preproc.trial_timing_dec{i};

            for j=1:size(trial_timing_raw_,1)
                trial_timing_raw__ = trial_timing_raw_{j};
                trial_timing_dec__ = trial_timing_dec_{j};

                trial_timing_raw__(:,'start') = table(trial_timing_raw__.start + samples_to_add_raw); 
                trial_timing_raw__(:,'end')   = table(trial_timing_raw__.end + samples_to_add_raw); 
                trial_timing_dec__(:,'start') = table(trial_timing_dec__.start + samples_to_add_dec); 
                trial_timing_dec__(:,'end')   = table(trial_timing_dec__.end + samples_to_add_dec); 

                trial_timing_raw_{j} = trial_timing_raw__;
                trial_timing_dec_{j} = trial_timing_dec__;
                fprintf(1,'.');
            end

            trial_timing_raw = [trial_timing_raw; trial_timing_raw_];
            trial_timing_dec = [trial_timing_dec; trial_timing_dec_];

            condition    = [condition; obj.condition{i}];
            session      = [session; obj.session{i}];
            events_table = [events_table; obj.events_table{i}];
        end

        fprintf(1,'] done\n');
            
        obj.for_preproc.trial_timing_raw = trial_timing_raw;
        obj.for_preproc.trial_timing_dec = trial_timing_dec;

        if obj.sample_freq == obj.for_preproc.decimation_freq
            obj.trial_timing = trial_timing_dec;
        else
            obj.trial_timing = trial_timing_raw;
        end

        obj.condition    = condition;
        obj.session      = session;
        obj.events_table = events_table;
    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHANNEL MAPPING TO ANATOMY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [chan_mapping,chan_labels]=channel_mapping_anatomical(obj,vera_mat)
        fprintf(1, '\n> Adding anatomical files to object\n');

        chan_labels = obj.elec_ch_label;
        chan_types = obj.elec_ch_type;

        if strcmp(obj.subject,'BJH006')
            reduced_localized_names = cellfun(@(x) split(x,'-'),vera_mat.electrodeNames,'UniformOutput',false);
            localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{end},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
            localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        elseif strcmp(obj.subject,'SLCH002')
            reduced_localized_names = cellfun(@(x) split(x,'^'),vera_mat.electrodeNames,'UniformOutput',false);
            localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{2},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
            localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        elseif contains(obj.subject,'BJH') 
            first_thing = cellfun(@(x) x(1:3),vera_mat.electrodeNames,'UniformOutput',false);
            second_thing = cellfun(@(x) extract(x,digitsPattern),vera_mat.electrodeNames,'UniformOutput',false); 
            second_thing = cellfun(@(x) x{end},second_thing,'UniformOutput',false);
            localized_names = cellfun(@(x,y) strcat(x([1,3]),'_',y),first_thing,second_thing,'UniformOutput',false);
            if strcmp(obj.subject,'BJH008')
                for i=1:length(localized_names)
                    if contains(localized_names{i},'ER')
                        localized_names(i) = {['O' localized_names{i}(3:end)]};
                    elseif contains(localized_names{i},'FR')
                        localized_names(i) = {['P' localized_names{i}(3:end)]};
                    else
                        localized_names(i) = {localized_names{i}([1,3:end])};
                    end
                end
            end
            chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
            chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

        else
            num_chans = sum(contains(chan_types,'ecog') | strcmp('seeg',chan_types));
            chan_mapping = num2cell(1:num_chans)';
            chan_labels = chan_labels(1:num_chans);
            empty_to_add = cell(size(obj.elec_data,1)-num_chans,1);
            chan_mapping = [chan_mapping; empty_to_add];
            chan_labels = [chan_labels; empty_to_add];
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKE TRIALS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function make_trials(obj)
        if ~isempty(obj.bip_elec_data)
            trial_keys = {'key','string','elec_data','bip_elec_data'};
        else
            trial_keys = {'key','string','elec_data'};
        end

        fprintf(1, '\n> Cutting signal into trial data ... \n');
        fprintf(1,'[');

        for k = 1:size(obj.trial_timing,1)
            trial_time_tbl = obj.trial_timing{k};
            trial_elec_data = arrayfun(@(x) obj.elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
        
            if ~isempty(obj.bip_elec_data)
                trial_bip_elec_data = arrayfun(@(x) obj.bip_elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
                obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,trial_bip_elec_data,'VariableNames',trial_keys);
            else
                obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,'VariableNames',trial_keys);
            end
            fprintf(1,'.');
        end
        fprintf(1,'] done\n');
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET CONDITION TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cond_data=get_cond_resp(obj,condition,varargin)
        p = inputParser();
        addParameter(p,'keep_trials',[]);
        parse(p, varargin{:});
        ops = p.Results;

        if isempty(obj.trial_data)
            obj.make_trials();
        end
        
        if ~isempty(ops.keep_trials)
            assert(length(obj.condition)==length(ops.keep_trials));
            cond_id = find(cell2mat(arrayfun(@(x) (strcmp(obj.condition{x},condition) && ops.keep_trials(x)),1:length(obj.condition),'UniformOutput',false)));
        else
            cond_id = find(cell2mat(arrayfun(@(x) strcmp(obj.condition{x},condition),1:length(obj.condition),'UniformOutput',false)));
        end

        cond_data = obj.trial_data(cond_id);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function input_d_ave=get_average(obj,input_d)
        input_d_ave = input_d;
        for k=1:size(input_d,1)
            B = input_d{k};
            [keys,strings,values] = obj.get_columns(B);
            values_ave = varfun( @(x) cellfun(@(y) nanmean(y,2), x,'uni',false), values,'OutputFormat','table');
            values_ave.Properties.VariableNames = values.Properties.VariableNames;
            input_d_ave{k} = [keys,strings,values_ave];
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET EVENTS BY KEY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_d=get_value(obj,input_d,varargin)
        p = inputParser();
        addParameter(p, 'key', 'word');
        addParameter(p, 'type', 'match');
        parse(p, varargin{:});
        ops = p.Results;

        if strcmp(ops.type,'match')
            func = @(x,y) ismember(x,y);
        else
            func = @(x,y) contains(x,y);
        end

        output_d = input_d;
        for k = 1:size(input_d,1)
            B = input_d{k};
            output_d{k} = B(func(B.key,ops.key),:);
        end
    end
 

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET TABLE COLUMNS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [keys,strings,values]=get_columns(obj,B)
        keys = B(:,ismember(B.Properties.VariableNames,'key'));
        strings = B(:,ismember(B.Properties.VariableNames,'string'));
        values = B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMBINE TRIAL CONDITIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_d=combine_trial_cond(obj,input_d)
        all_keys = cellfun(@(x) x.key,input_d,'uni',false);
        [X,Y] = ndgrid(1:numel(all_keys));
        Z = tril(true(numel(all_keys)),-1);
        assert(all(arrayfun(@(x,y) isequal(all_keys{x},all_keys{y}),X(Z),Y(Z))));
        func = @(x) cell2mat(reshape(vertcat(x),1,[]));
        output_d = table();
        for k=1:numel(all_keys{1})
            each_key = all_keys{1}{k};
            temp = (cellfun(@(x) x(ismember(x.key,each_key),:),input_d,'uni',false));
            cond_tbl = vertcat(temp{:});
            [~,strings,values] = obj.get_columns(cond_tbl);
            values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            temp_table = cell2table(horzcat(each_key,{strings.string},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
            output_d = [output_d; temp_table];
        end 
    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE SIGNAL FOR CONDITION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output_tbl,cond_table]=get_ave_cond_trial(obj,varargin)
        p = inputParser();
        addParameter(p,'words',1:8);
        addParameter(p,'condition',[]);
        addParameter(p,'keep_trials',[])
        parse(p, varargin{:});
        ops = p.Results;
            
        func = @(x) cell2mat(permute(x,[3,2,1]));
            
        if ops.condition
            condition_flag = ops.condition;
            cond_data = obj.get_cond_resp(condition_flag,'keep_trials',ops.keep_trials);
        else 
            if isempty(obj.trial_data)
                obj.make_trials();
            end
            condition_flag = 'all';
            cond_data = obj.trial_data;
        end
        
        cond_data_ave = obj.get_average(cond_data);
        word_data = obj.get_value(cond_data_ave,'key','word','type','contain');
            
        B = obj.combine_trial_cond(word_data);
        [keys,strings,values] = obj.get_columns(B);
        values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
        cond_table = cell2table(horzcat(condition_flag,{strings.string},values_comb),'VariableNames',B.Properties.VariableNames);

        func_1 = @(x) x{1}(:,:,ops.words);
        func_2 = @(x) nanmean(x,3);
        
        B = cond_table;
        [keys,strings,values] = obj.get_columns(B);
        condition_ave = cellfun(@(X) func_2(func_1(values.(X))),values.Properties.VariableNames,'uni',false);
        
        output_tbl = cell2table(horzcat(condition_flag,strings.string,condition_ave),'VariableNames',B.Properties.VariableNames);
    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET CONDITION IDs (logical index vector)
    % Returns a logical row vector, one entry per trial.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cond_id = get_cond_id(obj, condition, varargin)
        p = inputParser();
        addParameter(p,'keep_trials',[]);
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.keep_trials)
            assert(length(obj.condition)==length(ops.keep_trials));
            cond_id = cell2mat(arrayfun(@(x) (strcmp(obj.condition{x},condition) && ops.keep_trials(x)),...
                               1:length(obj.condition),'UniformOutput',false));
        else
            cond_id = cell2mat(arrayfun(@(x) strcmp(obj.condition{x},condition),...
                               1:length(obj.condition),'UniformOutput',false));
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT TRIAL EPOCHS
    % Returns a 3-D array [nChans x nTrials x nSamples].
    % Unlike make_trials() (which uses a cell-table format), this method
    % is designed for downstream normalization and significance testing.
    %
    %   epoch_tw        - [start stop] in seconds relative to event onset
    %   key             - event key to align to (default 'fix')
    %   selectChannels  - channel indices to extract (default: all)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [trial_data_epoch, trial_bip_data_epoch] = extract_trial_epochs(obj, varargin)
        p = inputParser();
        addParameter(p,'epoch_tw',  [-0.5 3]);
        addParameter(p,'key',       'fix');
        addParameter(p,'selectChannels', 1:size(obj.elec_data,1));
        parse(p, varargin{:});
        ops = p.Results;

        epoch_tw_samples = round(ops.epoch_tw .* obj.sample_freq);

        fprintf(1, '\n> Cutting signal into trial epochs ... \n');
        fprintf(1,'[');
        trial_bip_data_epoch = [];
        trial_data_epoch     = [];

        for k = 1:size(obj.trial_timing,1)
            trial_time_tbl = obj.trial_timing{k};
            probe_key = find(ismember(trial_time_tbl.key, ops.key));
            assert(length(probe_key)==1, 'Key ''%s'' not found or not unique in trial %d', ops.key, k);

            t_start = trial_time_tbl(probe_key,:).start + epoch_tw_samples(1);
            t_stop  = trial_time_tbl(probe_key,:).start + epoch_tw_samples(2);

            trial_data_epoch(:,k,:) = obj.elec_data(ops.selectChannels, t_start:t_stop);

            if ~isempty(obj.bip_elec_data)
                trial_bip_data_epoch(:,k,:) = obj.bip_elec_data(:, t_start:t_stop);
            end

            fprintf(1,'.');
        end
        fprintf(1,'] done\n');
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT NORMALIZATION METRICS
    % Computes per-channel [mean, std] from fixation/baseline epochs and
    % stores the result in obj.stats.normMetrics.
    %
    % Requires kumar_ieeg_utils/ (remove_bad_trials, extractCommonTrials)
    % to be on the MATLAB path.
    %
    %   baseTimeRange   - [start stop] seconds for baseline window (default [-0.5 0])
    %   timepad         - extra seconds added either side before epoch extraction
    %   key             - event key that marks baseline onset (default 'fix')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_normalization_metrics(obj, varargin)
        p = inputParser();
        addParameter(p,'baseTimeRange', [-0.5 0]);
        addParameter(p,'timepad',       0.5);
        addParameter(p,'key',           'fix');
        parse(p, varargin{:});
        ops = p.Results;

        baseTimeExtract = [ops.baseTimeRange(1)-ops.timepad, ...
                           ops.baseTimeRange(2)+ops.timepad];

        [baseData, baseData_bip] = obj.extract_trial_epochs(...
            'epoch_tw', baseTimeExtract, 'key', ops.key);

        [~, goodtrials]  = remove_bad_trials(baseData);
        goodTrialsCommon = extractCommonTrials(goodtrials);

        fprintf(1, '\n>> Extracting normalization metrics (unipolar) \n');
        fprintf(1,'[');
        normFactor = zeros(size(baseData,1), 2);
        for iChan = 1:size(baseData,1)
            slice = squeeze(baseData(iChan, goodTrialsCommon, :));
            normFactor(iChan,:) = [mean(slice(:)), std(slice(:))];
            fprintf(1,'.');
        end
        fprintf(1,'] done\n');
        normMetrics.normFactor = normFactor;

        if ~isempty(baseData_bip)
            fprintf(1, '\n>> Extracting normalization metrics (bipolar) \n');
            fprintf(1,'[');
            [~, goodtrials_bip]  = remove_bad_trials(baseData_bip);
            goodTrialsCommon_bip = extractCommonTrials(goodtrials_bip);
            normFactor_bip = zeros(size(baseData_bip,1), 2);
            for iChan = 1:size(baseData_bip,1)
                slice = squeeze(baseData_bip(iChan, goodTrialsCommon_bip, :));
                normFactor_bip(iChan,:) = [mean(slice(:)), std(slice(:))];
                fprintf(1,'.');
            end
            fprintf(1,'] done\n');
            normMetrics.normFactor_bip = normFactor_bip;
        end

        obj.stats.normMetrics = normMetrics;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NORMALIZE SIGNAL
    % Applies a normalization method to obj.elec_data (and bip_elec_data)
    % using metrics computed by extract_normalization_metrics().
    %
    %   normtype - one of:
    %     'z-score'    (x - mean) / std             (default)
    %     'mean-sub'   (x - mean)
    %     'perc-change'(x - mean) / mean
    %     'ratio'       x / mean
    %     'log-ratio'  10*log10(x / mean)
    %     'norm'       (x - mean) / (x + mean)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function normalize_signal(obj, varargin)
        p = inputParser();
        addParameter(p,'normtype','z-score');
        parse(p, varargin{:});
        ops = p.Results;

        if ~isfield(obj.stats,'normMetrics')
            obj.extract_normalization_metrics();
        end
        nM = obj.stats.normMetrics;

        obj.elec_data = apply_norm(obj.elec_data, nM.normFactor, ops.normtype);

        if ~isempty(obj.bip_elec_data) && isfield(nM,'normFactor_bip')
            obj.bip_elec_data = apply_norm(obj.bip_elec_data, nM.normFactor_bip, ops.normtype);
        end

        function data = apply_norm(data, fac, ntype)
            mu  = fac(:,1);
            sig = fac(:,2);
            switch ntype
                case 'z-score'
                    data = (data - mu) ./ sig;
                case 'mean-sub'
                    data = data - mu;
                case 'perc-change'
                    data = (data - mu) ./ mu;
                case 'ratio'
                    data = data ./ mu;
                case 'log-ratio'
                    data = 10 .* log10(data ./ mu);
                case 'norm'
                    data = (data - mu) ./ (data + mu);
                otherwise
                    error('Unknown normtype: %s', ntype);
            end
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT TIME SIGNIFICANCE
    % Runs permutation cluster test at every time point for every channel.
    % Results stored in obj.stats.time_series.pSigChan (and _bip).
    %
    % Requires kumar_ieeg_utils/ on the MATLAB path:
    %   timePermCluster, extendTimeEpoch
    %
    %   baseTime  - [start stop] seconds for baseline window
    %   epochTime - [start stop] seconds for analysis window
    %   p_val     - significance threshold (default 0.05)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_time_significance(obj, varargin)
        p = inputParser();
        addParameter(p,'baseTime',  [-0.5 0]);
        addParameter(p,'epochTime', [-0.5 6]);
        addParameter(p,'numPerm',   10000);
        addParameter(p,'p_val',     0.05);
        parse(p, varargin{:});
        ops = p.Results;

        [baseData,  baseData_bip]  = obj.extract_trial_epochs('epoch_tw', ops.baseTime);
        [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', ops.epochTime);

        baseDataExtend     = extendTimeEpoch(baseData,     size(epochData,3));
        baseData_bip_Ext   = extendTimeEpoch(baseData_bip, size(epochData,3));

        fprintf(1,'\n>> Permutation cluster test — unipolar (%d channels) ...\n', size(baseDataExtend,1));
        pSigChan = cell(1, size(baseDataExtend,1));
        parfor iChan = 1:size(baseDataExtend,1)
            aSig = squeeze(epochData(iChan,:,:));
            bSig = squeeze(baseDataExtend(iChan,:,:));
            pSigChan{iChan} = timePermCluster(aSig, bSig, 'pThresh', ops.p_val);
        end

        pSigChan_bip = {};
        if ~isempty(epochData_bip)
            fprintf(1,'\n>> Permutation cluster test — bipolar (%d channels) ...\n', size(baseData_bip_Ext,1));
            pSigChan_bip = cell(1, size(baseData_bip_Ext,1));
            parfor iChan = 1:size(baseData_bip_Ext,1)
                aSig = squeeze(epochData_bip(iChan,:,:));
                bSig = squeeze(baseData_bip_Ext(iChan,:,:));
                pSigChan_bip{iChan} = timePermCluster(aSig, bSig, 'pThresh', ops.p_val);
            end
        end

        obj.stats.time_series.pSigChan     = pSigChan;
        obj.stats.time_series.pSigChan_bip = pSigChan_bip;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT SIGNIFICANT CHANNELS
    % One-sided permutation test (epoch power > baseline power) per
    % channel, then FDR correction.
    % Results stored in obj.stats.sig_hg_channel.
    %
    % Requires on MATLAB path:
    %   remove_bad_trials  (kumar_ieeg_utils/)
    %   fdr_bh             (fdr_bh/)
    %
    %   baseTime  - [start stop] seconds baseline window
    %   epochTime - [start stop] seconds analysis window
    %   numPerm   - permutation count (default 10000)
    %   p_val     - threshold (default 0.05)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_significant_channel(obj, varargin)
        p = inputParser();
        addParameter(p,'baseTime',  [-0.5 0]);
        addParameter(p,'epochTime', [0 0.5]);
        addParameter(p,'numPerm',   10000);
        addParameter(p,'p_val',     0.05);
        parse(p, varargin{:});
        ops = p.Results;

        [baseData,  baseData_bip]  = obj.extract_trial_epochs('epoch_tw', ops.baseTime);
        [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', ops.epochTime);

        basePower  = mean(baseData.^2,  3);
        epochPower = mean(epochData.^2, 3);

        fprintf(1, '\n>> Extracting significant channels (unipolar) ...\n');
        fprintf(1,'[');
        [~, goodtrials] = remove_bad_trials(epochData);
        pSig.pChan = zeros(1, size(basePower,1));
        for iChan = 1:size(basePower,1)
            pSig.pChan(iChan) = ecog_data.local_permtest(...
                epochPower(iChan, goodtrials(iChan,:)), ...
                basePower(iChan,  goodtrials(iChan,:)), ...
                ops.numPerm);
            fprintf(1,'.');
        end
        fprintf(1,'] done\n');
        pSig.h_fdr_05 = fdr_bh(pSig.pChan, 0.05);
        pSig.h_fdr_01 = fdr_bh(pSig.pChan, 0.01);

        if ~isempty(epochData_bip)
            basePower_bip  = mean(baseData_bip.^2,  3);
            epochPower_bip = mean(epochData_bip.^2, 3);
            fprintf(1, '\n>> Extracting significant channels (bipolar) ...\n');
            fprintf(1,'[');
            [~, goodtrials_bip] = remove_bad_trials(epochData_bip);
            pSig.pChan_bip = zeros(1, size(basePower_bip,1));
            for iChan = 1:size(basePower_bip,1)
                pSig.pChan_bip(iChan) = ecog_data.local_permtest(...
                    epochPower_bip(iChan, goodtrials_bip(iChan,:)), ...
                    basePower_bip(iChan,  goodtrials_bip(iChan,:)), ...
                    ops.numPerm);
                fprintf(1,'.');
            end
            fprintf(1,'] done\n');
            pSig.h_bip_fdr_05 = fdr_bh(pSig.pChan_bip, 0.05);
            pSig.h_bip_fdr_01 = fdr_bh(pSig.pChan_bip, 0.01);
        end

        obj.stats.sig_hg_channel = pSig;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT CHANNELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_channels(obj,varargin)
        p = inputParser();
        addRequired(p,'signal');
        addRequired(p,'channel_labels');
        addRequired(p,'clean_channels');
        addRequired(p,'valid_channels');
        addParameter(p,'stitch_index',[1]);
        addParameter(p,'t_len',60);
        addParameter(p,'sample_freq',1200);
        addParameter(p,'downsample',false);
        addParameter(p,'decimation_freq',300)
        addParameter(p,'plotIEDs',false);
        addParameter(p,'chanIEDs',[]);
        addParameter(p,'posIEDs',[]);
        addParameter(p,'save',false); 
        parse(p, varargin{:});
        ops = p.Results;
        
        curr_sample_freq = ops.sample_freq;
        D = ops.signal;
        
        x_norm_cell = [];
        for k=1:length(ops.stitch_index)
            if k == length(ops.stitch_index)
                stop = size(D,1);
            else
                stop = ops.stitch_index(k+1)-1;
            end 
            D_ = D(ops.stitch_index(k):stop,:);
            x_cell = mat2cell(D_',ones(1,size(D_,2)));

            if ops.downsample
                decimation_factor = ops.sample_freq/ops.decimation_freq;
                x_cell = cellfun(@(x) downsample(x,decimation_factor),x_cell,'uni',false);
                curr_sample_freq = ops.decimation_freq;
            end

            min_max = mean(cell2mat(cellfun(@(y) prctile(y,[3 97]),x_cell,'UniformOutput',false)));
            x_norm_cell_ = cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
            x_norm_cell_ = arrayfun(@(x) x_norm_cell_{x}*ops.valid_channels(x),1:size(x_norm_cell_,1),'uni',false)';
            
            if length(ops.stitch_index) > 1
                x_norm_cell = [x_norm_cell, x_norm_cell_];
            else 
                x_norm_cell = x_norm_cell_;
            end
        end

        if length(ops.stitch_index) > 1
            for k=1:length(obj.stitch_index)-1
                x_norm_cell(:,k+1) = arrayfun(@(x) {[x_norm_cell{x,k}, x_norm_cell{x,k+1}]},[1:size(x_norm_cell,1)])';
            end
            x_norm_cell = x_norm_cell(:,length(ops.stitch_index));
        end
        assert(size(x_norm_cell,2)==1,'x_norm_cell not in the correct format');

        figure(1); clf;
        t_length = ops.t_len*curr_sample_freq;
        col_inf=inferno(floor(.8*size(x_norm_cell,1)));
        col_vir=viridis(floor(.8*size(x_norm_cell,1)));
        colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
        close all
        figure(1);
        clf;
        set(gcf,'position',[31,1,1700,900]);
        ax = axes('position',[.05,.1,.93,.88]);
        hold on
        time_stamps = [1:size(x_norm_cell{1},2)]/curr_sample_freq;
        hold on
        H = arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,ops.channel_labels{x})),[1:size(x_norm_cell,1)]);

        if ops.plotIEDs
            spike_chan = ops.chanIEDs;
            [spike_chan_sort,sort_idx] = sort(spike_chan);
            spike_times_sort = ops.posIEDs(sort_idx);
            yval = cell2mat(arrayfun(@(x) x_norm_cell{spike_chan_sort(x)}(floor(spike_times_sort(x)*200))+spike_chan_sort(x),1:length(spike_chan_sort),'uni',false));
            H1 = scatter(spike_times_sort,yval,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 .7 .7],'LineWidth',2,'MarkerFaceAlpha',.5);
        end
            
        set(gcf,'doublebuffer','on');
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','');
        arrayfun(@(x) text(0,x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),ops.clean_channels);
        set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
        ax.XAxis.TickLength = [0.005,0.01];
        ax.YAxis.TickLength = [0.005,0.01];
        set(ax,'xlim',[0 ops.t_len]);
        pos = get(ax,'position');
        Newpos = [pos(1) pos(2)-0.1 pos(3) 0.05];
        xmax=max(time_stamps);
        S = ['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(ops.t_len) '])'];
        h = uicontrol('style','slider','units','normalized','position',Newpos,'callback',S,'min',0,'max',xmax-ops.t_len);
        datacursormode on
    end

    
end % end methods block


methods (Static)

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NAPLAB FILTERBANK  (static — can also be called standalone)
    %
    %   [filteredData, cfs, sigma_fs] = ecog_data.naplab_filterbank(d, Fs, freqRange)
    %
    %   d          - [nChans x nSamples]  input signal
    %   Fs         - sampling rate (Hz)
    %   freqRange  - [fLow fHigh] (Hz), default [70 150]
    %
    %   filteredData - [nChans x nSamples x nBands]  complex analytic signal
    %   cfs          - [1 x nBands]  center frequencies (Hz)
    %   sigma_fs     - [1 x nBands]  Gaussian sigma (Hz)
    %
    % Source: Neural Acoustic Processing Lab, Columbia University
    %         naplab.ee.columbia.edu
    %         Embedded here to avoid external dependency.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONE-SIDED PERMUTATION TEST  (static helper for extract_significant_channel)
    %
    %   p = ecog_data.local_permtest(sample1, sample2, numperm)
    %
    %   Tests H0: mean(sample1) <= mean(sample2).
    %   Returns proportion of shuffled differences exceeding observed difference.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = local_permtest(sample1, sample2, numperm)
        samples    = [sample1, sample2];
        samplediff = mean(sample1) - mean(sample2);
        n1         = length(sample1);
        diffshuff  = zeros(1, numperm);
        for n = 1:numperm
            s = samples(randperm(length(samples)));
            diffshuff(n) = mean(s(1:n1)) - mean(s(n1+1:end));
        end
        p = sum(diffshuff > samplediff) / numperm;
    end


    function [filteredData, cfs, sigma_fs] = naplab_filterbank(d, Fs, freqRange)
        if nargin < 3 || isempty(freqRange)
            freqRange = [70 150];
        end

        a = [log10(.39); .5];
        f0       = 0.018;
        octspace = 1/7;
        minf = freqRange(1);
        maxf = freqRange(2);
        maxfo = log2(maxf/f0);

        cfs = f0;
        sigma_f = 10^(a(1) + a(2)*log10(cfs(end)));

        while log2(cfs(end)/f0) < maxfo
            cfo = log2(cfs(end)/f0) + octspace;
            if cfs(end) < 4
                cfs = [cfs, cfs(end) + sigma_f];
            else
                cfs = [cfs, f0*(2^cfo)];
            end
            sigma_f = 10^(a(1) + a(2)*log10(cfs(end)));
        end

        cfs = cfs(cfs >= minf & cfs <= maxf);
        npbs = length(cfs);
        sigma_fs = (10.^([ones(length(cfs),1), log10(cfs')]*a))';

        % Remove problematic frequency bands near power-line harmonics
        badfs    = [find(cfs>340 & cfs<480), find(cfs>720 & cfs<890)];
        sigma_fs = sigma_fs(setdiff(1:npbs, badfs));
        cfs      = cfs(setdiff(1:npbs, badfs));
        npbs     = length(cfs);
        sds      = sigma_fs .* sqrt(2);

        T      = size(d, 2);
        freqs  = (0:floor(T/2)) .* (Fs/T);
        nfreqs = length(freqs);

        % One-sided spectrum (Hilbert convention)
        h = zeros(1, T);
        if mod(T,2) == 0
            h([1, T/2+1]) = 1;
            h(2:T/2) = 2;
        else
            h(1) = 1;
            h(2:(T+1)/2) = 2;
        end

        filteredData = zeros(size(d,1), T, npbs);
        for c = 1:size(d,1)
            adat = fft(d(c,:), T);
            for f = 1:npbs
                H = zeros(1, T);
                k_vec = freqs - cfs(f);
                H(1:nfreqs) = exp((-0.5) .* ((k_vec ./ sds(f)).^2));
                H(nfreqs+1:end) = fliplr(H(2:ceil(T/2)));
                H(1) = 0;
                filteredData(c,:,f) = ifft(adat .* (H .* h), T);
            end
        end
        filteredData = abs(filteredData);
    end

end % end static methods block


end % end classdef
