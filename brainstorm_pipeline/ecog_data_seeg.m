classdef ecog_data_seeg < ecog_data_ieeg
% ECOG_DATA_SEEG  SEEG-specific subclass of the advanced ieeg_pipeline engine.
%
% This class re-bases the Brainstorm/MIT SEEG pipeline onto the more advanced
% EvLab ieeg_pipeline engine (vendored here as @ecog_data_ieeg) by INHERITING
% everything from that engine and overriding ONLY the methods that must differ
% for stereo-EEG (SEEG) data.
%
% Why subclass the vendored @ecog_data_ieeg instead of the upstream
% ieeg_pipeline-master/@ecog_data directly:
%   - The upstream class is named `ecog_data`, which already exists twice more
%     in this repo (./ecog_data.m and MGH_utils/ecog_data.m). Subclassing the
%     bare name `ecog_data` would resolve by MATLAB path order and could silently
%     pick a legacy class. The uniquely-named @ecog_data_ieeg removes that
%     collision so path ORDER does not matter.
%   - @ecog_data_ieeg is kept verbatim from upstream (only the class name and the
%     `arguments obj` type validators are renamed), so it can be re-synced from
%     ieeg_pipeline-master with a simple diff.
%
% Inherited advanced engine methods (NOT overridden here) now come from the
% ieeg_pipeline engine: extract_high_gamma, normalize_signal (with envelope
% smoothing), downsample_signal, make_trials, measure_line_noise, remove_IED,
% visual_inspection, extract_significant_channel, extract_time_significance,
% output_data_structures, output_xarray(_minimal), plus (concatenation), etc.
%
% Additional SEEG-tuned helpers kept on this subclass (thin variants of, or
% extras alongside, the engine equivalents):
%   - detect_sharp_artifacts : sharp-transient flagging on the z-scored HG
%                          envelope (inputParser-based; avoids parfor).
%   - smooth_high_gamma  : explicit Gaussian smoothing of the HG envelope (the
%                          engine also smooths inside normalize_signal).
%   - extract_bandpass_signal : segment-wise bandpass; errors clearly if eegfilt
%                          (EEGLAB) is not on the path.
%   - saveUpdatedObject  : convenience save of the processed object.
%
% SEEG-specific overrides (everything else is inherited from ecog_data_ieeg):
%   - define_parameters  : 50 Hz line-noise standard (peak [45 50 55] Hz, notch
%                          50/100/150/... Hz) instead of the 60 Hz US standard,
%                          and excludes Gaussian high-gamma bands that fall on
%                          50 Hz harmonics.
%   - notch_filter       : reports the actual line-noise frequency (from
%                          for_preproc.filter_params.line_noise_hz) and runs the
%                          interactive noisy-channel review.
%   - plot_line_noise    : 50 Hz axis labels; figure name from crunched_file_name.
%   - extract_shanks     : operates only on clean channels, so excluded channels
%                          with unparseable labels cannot break shank parsing.
%   - reference_signal   : derives bipolar pairs directly from channel labels
%                          (consistent with the clean-only extract_shanks) and
%                          adds along-shank LAPLACIAN referencing
%                          (doLaplacianReferencing).
%   - preprocess_signal  : SEEG preprocessing orders that do NOT apply CAR before
%                          bipolar referencing, plus Laplacian-based orders
%                          ('defaultSEEGLaplacian', 'preEnvelopeExtractionSEEGLaplacian').
%                          Unrecognized orders are delegated to the engine.
%   - extract_trial_epochs: string-key epoching ('key','word_1', ...) with the
%                          window rounded to whole samples (the engine uses a
%                          numeric probe_key instead).
%   - extract_normalization_metrics : retains the `key` argument so the baseline
%                          can be anchored to a word onset (the engine dropped it).
%   - get_cond_id / get_cond_resp / get_value / get_ave_cond_trial : kept for the
%                          string condition flags and the exact table layout the
%                          S-vs-N analysis layer (ecog_sn_data_seeg) consumes.
%
% Crunched .mat files produced by brainstorm_to_mit_crunched_new.m contain a
% variable named 'obj' of this class type.

methods
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ecog_data_seeg(varargin)
        % Same constructor signature as ecog_data_ieeg; just forwards arguments.
        obj@ecog_data_ieeg(varargin{:});
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESS SIGNAL  (SEEG orders, no CAR before bipolar referencing)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function preprocess_signal(obj,varargin)
        p = inputParser();
        addParameter(p,'order','defaultSEEG');
        addParameter(p,'isPlotVisible',true);
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        % SEEG preprocessing orders. CAR is intentionally NOT included before
        % bipolar referencing for SEEG data (see class header).
        switch ops.order
            case 'defaultSEEG'
                order = {'highpassFilter',...
                         'notchFilter',...
                         'IEDRemoval',...
                         'visualInspection',...
                         'BipolarReferencing',...
                         'GaussianFilterExtraction',...
                         'removeOutliers',...
                         'downsample'...
                };

            case 'defaultSEEGbyShank'
                order = {'highpassFilter',...
                         'notchFilter',...
                         'IEDRemoval',...
                         'visualInspection',...
                         'ShankCSR',...
                         'BipolarReferencing',...
                         'GaussianFilterExtraction',...
                         'removeOutliers',...
                         'downsample'...
                };

            case 'preEnvelopeExtractionSEEG'
                order = {'highpassFilter',...
                         'notchFilter',...
                         'IEDRemoval',...
                         'visualInspection',...
                         'BipolarReferencing',...
                         'downsample'...
                };

            case 'defaultSEEGLaplacian'
                % Local Laplacian (along-shank) referencing instead of bipolar.
                % Ported from the ieeg_pipeline advanced preprocessing.
                order = {'highpassFilter',...
                         'notchFilter',...
                         'IEDRemoval',...
                         'visualInspection',...
                         'LaplacianReferencing',...
                         'GaussianFilterExtraction',...
                         'removeOutliers',...
                         'downsample'...
                };

            case 'preEnvelopeExtractionSEEGLaplacian'
                order = {'highpassFilter',...
                         'notchFilter',...
                         'IEDRemoval',...
                         'visualInspection',...
                         'LaplacianReferencing',...
                         'downsample'...
                };

            case 'test'
                order = {'downsample'};

            otherwise
                % Defer ECoG / engine orders (defaultECOG, defaultSEEGorBOTH,
                % broadband variants, etc.) to the ieeg_pipeline engine.
                preprocess_signal@ecog_data_ieeg(obj,varargin{:});
                return
        end

        obj.for_preproc.order = order;
        obj.for_preproc.isPlotVisible = ops.isPlotVisible;

        % mapping from preprocessing name to its method
        names = {'highpassFilter',...
                 'notchFilter',...
                 'IEDRemoval',...
                 'visualInspection',...
                 'GlobalMeanRemoval',...
                 'CAR',...
                 'ShankCSR',...
                 'LaplacianReferencing',...
                 'BipolarReferencing',...
                 'GaussianFilterExtraction',...
                 'BandpassExtraction',...
                 'NapLabFilterExtraction',...
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
                     'reference_signal',...
                     'extract_high_gamma',...
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
                while 1
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
    % NOTCH FILTER  (reports actual line-noise frequency)
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

        ln_hz = obj.for_preproc.filter_params.line_noise_hz;
        fprintf(1,'\nReduced %d Hz noise from %.2f to %.2f uV\n', ln_hz, mean(signal_noise_before(obj.elec_ch_clean,2)),mean(signal_noise_after(obj.elec_ch_clean,2)));
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
    % PLOT LINE NOISE
    % Labels use the actual line-noise frequency (line_noise_hz) and the
    % saved-figure name is derived from crunched_file_name. The ecog_data_v2
    % version references for_preproc.log_file_name, which the brainstorm
    % converter never sets.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f=plot_line_noise(obj,noise_before,noise_after)
        ln_hz  = obj.for_preproc.filter_params.line_noise_hz;
        lo_hz  = ln_hz - 5;   % peak filter centre - 5 Hz
        hi_hz  = ln_hz + 5;   % peak filter centre + 5 Hz

        fprintf(1, '\n> Plotting %d Hz noise power ...\n', ln_hz);
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
        legend({sprintf('%dHz noise',lo_hz), sprintf('%dHz noise',ln_hz), ...
                sprintf('%dHz noise',hi_hz), 'MARKED NOISY'}, ...
               'Location','best','FontSize',16,'Box','off');
        ylabel('Noise (uV)','FontSize',18);
        title('BEFORE NOTCH FILTERING','FontSize',22);
        obj.update_position(currsub);

        currsub = subplot(2,2,3);
        stem(noise_before(:,2)./mean(noise_before(:,[1,3]),2),'filled','Color',c);
        axis tight; hold on;
        stem(x,noise_before(idxs,2)./mean(noise_before(idxs,[1,3]),2),'filled','Color','k')
        xlabel('Channel #','FontSize',18);
        ylabel(sprintf('%dHz noise / mean %dHz+%dHz noise', ln_hz, lo_hz, hi_hz),'FontSize',18)
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
    % REFERENCE SIGNAL  (label-driven bipolar pairs; consistent with the
    % clean-only extract_shanks override below)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reference_signal(obj,varargin)
        p = inputParser();
        addParameter(p,'doGlobalMeanRemoval',false)
        addParameter(p,'doCAR',false);
        addParameter(p,'doShankCSR',false);
        addParameter(p,'doLaplacianReferencing',false);
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

            % COMMON AVERAGE REFERENCING (ECoG only)
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

            % LAPLACIAN REFERENCING (SEEG only)
            % Ported from the ieeg_pipeline. Re-references each clean contact to
            % its neighbour(s) along the same shank: endpoints use the single
            % adjacent contact, interior contacts subtract the mean of both
            % neighbours. Operates only on clean channels (consistent with the
            % SEEG extract_shanks override above).
            if ops.doLaplacianReferencing && ~isempty(seeg_chans)
                fprintf(1, '\n>> Performing Laplacian referencing for SEEG \n');
                fprintf(1,'[');

                signal_ = signal(obj.stitch_index(k):stop,:);
                signal_laplace = signal_;   % channels on skipped shanks keep original values

                [shank_locs,~,~] = obj.extract_shanks();
                for idx_shk = 1:length(shank_locs)
                    % sorted ascending by channel index (≈ contact order along shank)
                    shank_channels = intersect(shank_locs{idx_shk}, obj.elec_ch_clean);
                    if length(shank_channels) >= 2
                        for i = 1:length(shank_channels)
                            curr_channel = shank_channels(i);
                            if i == 1
                                next_channel = shank_channels(i+1);
                                laplacian = signal_(:,curr_channel) - signal_(:,next_channel);
                            elseif i == length(shank_channels)
                                prev_channel = shank_channels(i-1);
                                laplacian = signal_(:,curr_channel) - signal_(:,prev_channel);
                            else
                                prev_channel = shank_channels(i-1);
                                next_channel = shank_channels(i+1);
                                laplacian = signal_(:,curr_channel) - ...
                                            0.5 * (signal_(:,prev_channel) + signal_(:,next_channel));
                            end
                            signal_laplace(:,curr_channel) = laplacian;
                            fprintf(1,'.');
                        end
                    else
                        fprintf(1, '\nWarning: shank %d has fewer than 2 clean channels. Skipping Laplacian referencing for this shank.\n', idx_shk);
                    end
                end

                signal(obj.stitch_index(k):stop,:) = signal_laplace;
                fprintf(1,'] done\n');
            elseif ops.doLaplacianReferencing
                error('No SEEG channels to perform Laplacian referencing on')
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

        if ops.doGlobalMeanRemoval || ops.doCAR || ops.doShankCSR || ops.doLaplacianReferencing
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
    % EXTRACT SHANKS  (clean channels only, to avoid label-parse failures
    % on excluded channels)
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
    % DEFINE PARAMETERS  (50 Hz line-noise standard)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_parameters(obj)
        % ---------------------------------------------------------------
        % LINE NOISE STANDARD: 50 Hz  (European / Chinese / Asian)
        %   peak filter   : [45, 50, 55] Hz  - measures noise power
        %   notch filter  : 50, 100, 150, 200, 250 Hz - removes harmonics
        %   Gaussian bank : centres within +/-5 Hz of any 50 Hz harmonic
        %                   that falls inside the gamma band are excluded
        %
        % To switch to 60 Hz (US standard) replace:
        %   param.peak.fcenter  = [55, 60, 65];
        %   param.notch.fcenter = [60, 120, 180, 240];
        %   line_noise_hz       = 60;
        % ---------------------------------------------------------------
        param = struct;
        param.line_noise_hz = 50;           % <-- single source of truth

        % --- highpass filter ---
        param.highpass.Wp = 0.50;           % Hz
        param.highpass.Ws = 0.05;           % Hz
        param.highpass.Rp = 3;              % dB
        param.highpass.Rs = 30;             % dB

        % --- IIR peak filter (for line-noise QC measurement) ---
        % Measures power at {line_noise - 5, line_noise, line_noise + 5} Hz
        param.peak.fcenter = param.line_noise_hz + [-5, 0, 5];   % [45 50 55]
        param.peak.bw      = ones(1,3) .* 0.001;

        % --- notch filter: fundamental + harmonics up to Nyquist/2 ---
        max_notch = min(obj.sample_freq/2 - 50, 250);
        param.notch.fcenter = param.line_noise_hz : param.line_noise_hz : max_notch;
        param.notch.bw      = ones(1, length(param.notch.fcenter)) .* 0.001;

        % --- high-gamma band for Gaussian and bandpass extraction ---
        param.gaussian.f_gamma_low  = 70;
        param.gaussian.f_gamma_high = 150;

        param.bandpass.f_gamma_low  = 70;
        param.bandpass.f_gamma_high = 150;
        param.bandpass.filter_order = 6;


        % ----- BUILD FILTERS -----

        % --- highpass ---
        highpass.Wp = param.highpass.Wp/(obj.sample_freq/2);
        highpass.Ws = param.highpass.Ws/(obj.sample_freq/2);
        highpass.Rp = param.highpass.Rp;
        highpass.Rs = param.highpass.Rs;
        [highpass.n,highpass.Wn] = buttord(highpass.Wp,highpass.Ws,highpass.Rp,highpass.Rs);
        highpass.n = highpass.n + rem(highpass.n,2);
        [highpass.z,highpass.p,highpass.k] = butter(highpass.n,highpass.Wn,'high');
        [highpass.sos,highpass.g] = zp2sos(highpass.z,highpass.p,highpass.k);
        highpass.h = dfilt.df2sos(highpass.sos,highpass.g);

        % --- IIR peak (noise measurement) ---
        for idx = 1:length(param.peak.fcenter)
            peak{idx}.wo = param.peak.fcenter(idx)/(obj.sample_freq/2);
            peak{idx}.bw = param.peak.bw(idx);
            [peak{idx}.b,peak{idx}.a] = iirpeak(peak{idx}.wo,peak{idx}.bw);
        end

        % --- notch (line-noise removal) ---
        for idx = 1:length(param.notch.fcenter)
            notch{idx}.wo = param.notch.fcenter(idx)/(obj.sample_freq/2);
            notch{idx}.bw = param.notch.bw(idx);
            [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);
        end

        % --- Chang-lab Gaussian filterbank ---
        [gaussian.cfs, gaussian.sds] = obj.get_filter_param_chang_lab(...
            param.gaussian.f_gamma_low, param.gaussian.f_gamma_high);

        % Exclude Gaussian centres within +/-5 Hz of any 50 Hz harmonic
        % that falls inside the gamma band (e.g. 100 Hz for 50 Hz datasets).
        % This prevents residual line-noise from leaking into the HG envelope.
        harmonics_in_band = param.line_noise_hz * ...
            (ceil(param.gaussian.f_gamma_low  / param.line_noise_hz) : ...
             floor(param.gaussian.f_gamma_high / param.line_noise_hz));
        bad_bands = false(size(gaussian.cfs));
        for h = harmonics_in_band
            bad_bands = bad_bands | (abs(gaussian.cfs - h) < 5);
        end
        gaussian.cfs = gaussian.cfs(~bad_bands);
        gaussian.sds = gaussian.sds(~bad_bands);
        fprintf(1,'Gaussian filterbank: %d bands (excluded %d near 50 Hz harmonics)\n', ...
            numel(gaussian.cfs), sum(bad_bands));

        gaussian.f_gamma_low  = param.gaussian.f_gamma_low;
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
    % TRIAL / CONDITION HELPERS
    %
    % These are kept as SEEG overrides (rather than inheriting the engine
    % versions) so the S-vs-N analysis layer (ecog_sn_data_seeg) sees exactly
    % the table layout it expects, with string condition flags (e.g.
    % 'SENTENCES') and the plural 'words' name-value. They form a self-consistent
    % set (get_ave_cond_trial calls get_cond_resp and get_value below).
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


    function cond_id = get_cond_id(obj, condition, varargin)
        % Returns a logical row vector, one entry per trial.
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


    function [output_tbl,cond_table]=get_ave_cond_trial(obj,varargin)
        % Averages signal for a given condition across selected words.
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
        [keys,strings,values] = obj.get_columns(B); %#ok<ASGLU>
        values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
        cond_table = cell2table(horzcat(condition_flag,{strings.string},values_comb),'VariableNames',B.Properties.VariableNames);

        func_1 = @(x) x{1}(:,:,ops.words);
        func_2 = @(x) nanmean(x,3);

        B = cond_table;
        [keys,strings,values] = obj.get_columns(B); %#ok<ASGLU>
        condition_ave = cellfun(@(X) func_2(func_1(values.(X))),values.Properties.VariableNames,'uni',false);

        output_tbl = cell2table(horzcat(condition_flag,strings.string,condition_ave),'VariableNames',B.Properties.VariableNames);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT TRIAL EPOCHS
    % Returns a 3-D array [nChans x nTrials x nSamples]. Rounds the epoch
    % window to whole samples so non-integer windows cannot produce
    % non-integer indices.
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
    % EXTRACT NORMALIZATION METRICS  (retains the `key` baseline anchor)
    % The ieeg_pipeline engine dropped the `key` argument and samples the
    % baseline relative to a fixed probe index. The brainstorm trial_timing
    % tables have no 'fix' marker, so the baseline must be anchored to a named
    % key (e.g. 'word_1') with the window taken just before it. This mirrors the
    % MGH_utils/@ecog_data_v2 behaviour and feeds obj.stats.normMetrics, which
    % the inherited normalize_signal consumes.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function extract_normalization_metrics(obj,varargin)
        p = inputParser();
        addParameter(p,'baseTimeRange',[-0.5 0]);
        addParameter(p,'timepad',0.5);
        addParameter(p,'key','fix');
        parse(p, varargin{:});
        epoch_args = p.Results;

        timePad         = epoch_args.timepad;
        baseTimeRange   = epoch_args.baseTimeRange;
        baseTimeExtract = [baseTimeRange(1)-timePad baseTimeRange(2)+timePad];

        [baseData,baseData_bip] = obj.extract_trial_epochs('epoch_tw',baseTimeExtract,'key',epoch_args.key);
        [~,goodtrials] = remove_bad_trials(baseData);
        goodTrialsCommon = extractCommonTrials(goodtrials);

        fprintf(1, '\n>> Extracting normalization metrics for unipolar high gamma envelope \n');
        normFactor = zeros(size(baseData, 1), 2);
        fprintf(1,'[');
        for iChan = 1:size(baseData, 1)
            normFactor(iChan, :) = [mean2(squeeze(baseData(iChan, goodTrialsCommon, :))), std2(squeeze(baseData(iChan, goodTrialsCommon, :)))];
            fprintf(1,'.');
        end
        fprintf(1,'] done\n')
        normMetrics.normFactor = normFactor;

        if(~isempty(baseData_bip))
            fprintf(1, '\n>> Extracting normalization metrics for bipolar high gamma envelope \n');
            normFactor = zeros(size(baseData_bip, 1), 2);
            [~,goodtrials] = remove_bad_trials(baseData_bip);
            goodTrialsCommon = extractCommonTrials(goodtrials);
            fprintf(1,'[');
            for iChan = 1:size(baseData_bip, 1)
                normFactor(iChan, :) = [mean2(squeeze(baseData_bip(iChan, goodTrialsCommon, :))), std2(squeeze(baseData_bip(iChan, goodTrialsCommon, :)))];
                fprintf(1,'.');
            end
            fprintf(1,'] done\n')
            normMetrics.normFactor_bip = normFactor;
        end

        obj.stats.normMetrics = normMetrics;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DETECT SHARP ARTIFACTS  (SEEG variant of the ieeg_pipeline engine method)
    % Flags sharp inter-ictal transients on the z-scored high-gamma envelope
    % of bipolar (if present) and unipolar channels, using OR logic between an
    % amplitude criterion and a slope criterion. Results are stored in
    % obj.stats.artifact_stats_unipolar / .artifact_stats_bipolar.
    %
    % NOTE: run this AFTER normalize_signal (z-score), since the thresholds are
    % expressed in z-score units. For the MIT Naturalistic Stories task the
    % detector restricts analysis to story epochs.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function detect_sharp_artifacts(obj, varargin)
        p = inputParser();
        addParameter(p,'min_amplitude',15);   % z-score units
        addParameter(p,'min_slope',10);       % z-score units per sample
        parse(p, varargin{:});
        ops = p.Results;

        fs = obj.sample_freq;

        % decide which samples to analyse
        isStories = ~isempty(obj.experiment) && ...
                    strcmp(obj.experiment,'MITNaturalisticStoriesTask');
        if isStories
            epochs = [];
            for t = 1:size(obj.trial_timing,1)
                tbl = obj.trial_timing{t,1};
                if any(startsWith(string(tbl.key),"story_"))
                    storyRows = startsWith(string(tbl.key),"story_");
                    epochs = [epochs ; tbl.start(storyRows) tbl.end(storyRows)]; %#ok<AGROW>
                end
            end
            if isempty(epochs)
                warning('No story epochs found - analysing full recording instead');
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
            epochMask = true(1,size(obj.elec_data,2));
        end

        if ~isempty(obj.bip_elec_data)
            obj.stats.artifact_stats_bipolar = ...
                run_sharp_artifact_detector(obj.bip_elec_data(:,epochMask), fs, ops);
        end
        obj.stats.artifact_stats_unipolar = ...
            run_sharp_artifact_detector(obj.elec_data(:,epochMask), fs, ops);
        obj.stats.epochMask = epochMask;

        n_uni = sum([obj.stats.artifact_stats_unipolar.artifact_count] > 0);
        fprintf(1,'\nSharp-artifact detection: %d/%d unipolar channels flagged\n', ...
            n_uni, numel(obj.stats.artifact_stats_unipolar));
        if ~isempty(obj.bip_elec_data)
            n_bip = sum([obj.stats.artifact_stats_bipolar.artifact_count] > 0);
            fprintf(1,'Sharp-artifact detection: %d/%d bipolar channels flagged\n', ...
                n_bip, numel(obj.stats.artifact_stats_bipolar));
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SMOOTH HIGH GAMMA  (ported from ieeg_pipeline normalize step)
    % Gaussian smoothing (default 100 ms window) of the high-gamma envelope.
    % The ieeg_pipeline applies this inside normalize_signal; here it is kept
    % as an explicit, opt-in step so the inherited normalize_signal behaviour
    % is preserved by default.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function smooth_high_gamma(obj,varargin)
        p = inputParser();
        addParameter(p,'window_s',0.1);   % smoothing window in seconds
        parse(p, varargin{:});
        ops = p.Results;

        windowSize = max(1, round(ops.window_s * obj.sample_freq));
        fprintf(1,'\n> Gaussian-smoothing high gamma (%.0f ms window) ...\n', ops.window_s*1000);

        if ~isempty(obj.elec_data)
            obj.elec_data = smoothdata(obj.elec_data', 'gaussian', windowSize)';
        end
        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = smoothdata(obj.bip_elec_data', 'gaussian', windowSize)';
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT BANDPASS SIGNAL  (SEEG variant of the ieeg_pipeline engine method)
    % Segment-wise (stitch-aware) bandpass filtering of the raw/continuous
    % signal. Overwrites obj.elec_data (and obj.bip_elec_data, if present) with
    % the filtered signal and also returns it.
    %
    % REQUIRES eegfilt (EEGLAB) on the MATLAB path.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [filtered_signal, filtered_signal_bipolar] = extract_bandpass_signal(obj, low_cutoff, high_cutoff)
        if exist('eegfilt','file') ~= 2
            error(['extract_bandpass_signal requires eegfilt (EEGLAB) on the ' ...
                   'MATLAB path. Add EEGLAB or use GaussianFilterExtraction instead.']);
        end

        signal = double(obj.elec_data');
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        else
            signal_bipolar = [];
        end

        sample_freq  = obj.sample_freq;
        stitch_index = obj.stitch_index(:);

        filtered_signal = nan(size(signal));
        if ~isempty(signal_bipolar)
            filtered_signal_bipolar = nan(size(signal_bipolar));
        else
            filtered_signal_bipolar = [];
        end

        for k = 1:length(stitch_index)
            fprintf(1, '\n> Bandpass filtering segment %d of %d ... \n', k, length(stitch_index));
            if k == length(stitch_index)
                stop = size(signal,1);
            else
                stop = stitch_index(k+1)-1;
            end
            seg_idx = stitch_index(k):stop;

            segment = signal(seg_idx, :);
            filtered_segment = eegfilt(segment', sample_freq, low_cutoff, high_cutoff, 0, 200)';
            filtered_signal(seg_idx, :) = filtered_segment;

            if ~isempty(signal_bipolar)
                segment_bip = signal_bipolar(seg_idx, :);
                filtered_segment_bip = eegfilt(segment_bip', sample_freq, low_cutoff, high_cutoff, 0, 200)';
                filtered_signal_bipolar(seg_idx, :) = filtered_segment_bip;
            end
        end

        obj.elec_data = filtered_signal';
        if ~isempty(filtered_signal_bipolar)
            obj.bip_elec_data = filtered_signal_bipolar';
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE UPDATED OBJECT  (mirrors the ieeg_pipeline engine method)
    % Convenience save of the processed object to
    % <crunched_file_path>/<subject>_<experiment>_crunched_HG_ZScore.mat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function saveUpdatedObject(obj)
        crunchedFolder = fullfile(obj.crunched_file_path);
        if ~isempty(crunchedFolder) && ~exist(crunchedFolder, 'dir')
            mkdir(crunchedFolder);
        end
        filename = fullfile(crunchedFolder, [obj.subject '_' obj.experiment '_crunched_HG_ZScore.mat']);
        save(filename, 'obj', '-v7.3');
        fprintf('Updated object saved as: %s\n', filename);
    end

end % methods

end % classdef


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS  (visible to the class methods above)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats = run_sharp_artifact_detector(data, fs, ops)
% Detect sharp artefacts on z-scored high-gamma envelope data using OR logic:
% an event is kept if it satisfies EITHER the slope OR the amplitude criterion.
%
% Inputs
%   data : [nChan x nSamples] z-scored high-gamma envelope traces
%   fs   : sampling frequency (Hz)
%   ops  : struct with fields min_amplitude, min_slope
%
% Output (per-channel struct array)
%   channel_idx, artifact_count, max_amplitude, mean_slope, artifact_times

    chan_ids = 1:size(data,1);
    n        = numel(chan_ids);

    tmpl  = struct('channel_idx',[],'artifact_count',0,'max_amplitude',NaN,...
                   'mean_slope',NaN,'artifact_times',[]);
    stats = repmat(tmpl,n,1);

    for k = 1:n
        x  = data(k,:);
        dx = diff(x);   % slope (z-score/sample)

        [slopePkAmp, slopePkLoc] = findpeaks(abs(dx), 'MinPeakHeight', ops.min_slope); %#ok<ASGLU>
        [ampPkAmp,   ampPkLoc]   = findpeaks(abs(x),  'MinPeakHeight', ops.min_amplitude); %#ok<ASGLU>

        allLoc = [slopePkLoc, ampPkLoc];
        if ~isempty(allLoc)
            uniqueLoc = unique(allLoc);
            keepLoc   = sort(uniqueLoc);
            keepSlope = abs(dx(min(keepLoc, numel(dx))));
        else
            keepLoc   = [];
            keepSlope = [];
        end

        stats(k).channel_idx    = chan_ids(k);
        stats(k).artifact_count = numel(keepLoc);
        stats(k).artifact_times = keepLoc/fs;
        if isempty(keepLoc)
            stats(k).max_amplitude = NaN;
            stats(k).mean_slope    = NaN;
        else
            stats(k).max_amplitude = max(abs(x(keepLoc)));
            stats(k).mean_slope    = mean(keepSlope,'omitnan');
        end
    end
end
