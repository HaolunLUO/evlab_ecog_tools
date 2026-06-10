classdef ecog_data_seeg < ecog_data_v2
% ECOG_DATA_SEEG  SEEG-specific subclass of ecog_data_v2.
%
% This class keeps the Brainstorm/MIT SEEG pipeline aligned with the canonical
% MGH_utils/@ecog_data_v2 class by INHERITING everything from ecog_data_v2 and
% overriding ONLY the methods that must differ for stereo-EEG (SEEG) data.
%
% Why a subclass (instead of editing ecog_data_v2 directly or keeping a full
% standalone copy):
%   - ecog_data_v2 is shared by the ECoG (grid/strip) pipelines; editing it in
%     place would risk those analyses.
%   - A uniquely-named subclass removes the previous class-name collision (there
%     used to be three different `ecog_data` classes on the path, so the
%     brainstorm copy was silently shadowed and never actually ran).
%   - Inheriting from ecog_data_v2 means non-SEEG behavior stays automatically in
%     sync with the canonical class (no copy-drift).
%
% SEEG-specific overrides (everything else is inherited from ecog_data_v2):
%   - define_parameters  : 50 Hz line-noise standard (peak [45 50 55] Hz, notch
%                          50/100/150/... Hz) instead of the 60 Hz US standard,
%                          and excludes Gaussian high-gamma bands that fall on
%                          50 Hz harmonics.
%   - notch_filter       : reports the actual line-noise frequency (from
%                          for_preproc.filter_params.line_noise_hz).
%   - extract_shanks     : operates only on clean channels, so excluded channels
%                          with unparseable labels cannot break shank parsing.
%   - reference_signal   : derives bipolar pairs directly from channel labels,
%                          keeping it consistent with the clean-only
%                          extract_shanks above.
%   - preprocess_signal  : SEEG preprocessing orders that do NOT apply CAR
%                          before bipolar referencing (bipolar referencing
%                          already removes shared/common signal between adjacent
%                          contacts, so a prior common-average step is both
%                          unnecessary and can distort the local estimate).
%
% Crunched .mat files produced by brainstorm_to_mit_crunched_new.m contain a
% variable named 'obj' of this class type.

methods
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ecog_data_seeg(varargin)
        % Same constructor signature as ecog_data_v2; just forwards arguments.
        obj@ecog_data_v2(varargin{:});
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

            case 'test'
                order = {'downsample'};

            otherwise
                % Defer ECoG / legacy orders (including any that intentionally
                % use CAR) to the canonical ecog_data_v2 implementation.
                preprocess_signal@ecog_data_v2(obj,varargin{:});
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
    % REFERENCE SIGNAL  (label-driven bipolar pairs; consistent with the
    % clean-only extract_shanks override below)
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
    % These override the ecog_data_v2 versions, which are validated against
    % the legacy `ecog_data` class (so they reject an ecog_data_v2 subclass)
    % and/or assume numeric condition codes. The versions below accept string
    % condition flags (e.g. 'SENTENCES') and the plural 'words' name-value
    % used by the S-vs-N analysis layer.
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

end % methods

end % classdef
