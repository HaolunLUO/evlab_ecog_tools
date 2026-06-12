classdef ecog_data_ieeg < dynamicprops
% ECOG_DATA_IEEG  Advanced iEEG/SEEG preprocessing engine.
%
% This is a uniquely-named, vendored copy of the EvLab ieeg_pipeline core class
% (originally ieeg_pipeline-master/@ecog_data, Duraivel/Casto). It is renamed to
% ecog_data_ieeg ONLY to avoid the class-name collision with the legacy
% ecog_data classes elsewhere in this repo (./ecog_data.m, MGH_utils/ecog_data.m);
% the implementation is otherwise kept verbatim so it can be re-synced from
% upstream with a simple diff.
%
% The SEEG pipeline class ecog_data_seeg subclasses this engine and overrides
% only the methods that must differ for stereo-EEG (see ecog_data_seeg.m).
% TODO - SUMMARY

properties
    %% ---- DATA ----
    elec_data
    bip_elec_data
    stitch_index
    sample_freq
    for_preproc             % preproc
    trial_data              % trial
    stats
    recording_duration
   
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
    processing_history
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
    
end


methods
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = ecog_data_ieeg(...
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
        
        % Construct an instance of this class

        %% ---- DATA ----
        obj.for_preproc=for_preproc;

        %% ---- INFO ----
        obj.subject=subject;
        obj.experiment=experiment;
        obj.crunched_file_name=crunched_file_name;
        obj.crunched_file_path=crunched_file_path;
        obj.raw_file_name=raw_file_name;
        obj.raw_file_path=raw_file_path;

        %% ---- LABELS ----
        obj.elec_ch=elec_ch;
        obj.elec_ch_label=elec_ch_label;
        obj.elec_ch_prelim_deselect=elec_ch_prelim_deselect;
        obj.elec_ch_type=elec_ch_type;
        
    end


    %% FUNCTION TO RUN PIPELINE


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESS SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function preprocess_signal(obj,varargin)
        % TODO - description

        p = inputParser();
        addParameter(p,'order','defaultBOTH');
        addParameter(p,'isPlotVisible',true);
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;


        % ---------------------------
        % DEFINE PREPROCESSING ORDER
        % ---------------------------
        % this pipeline was constructed for the 'default' ordering
        % all new preprocessing orders should be checked for problems
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
                     'GlobalMeanRemoval',... % could z-score here
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
                     'GlobalMeanRemoval',... % could z-score here
                     'IEDRemoval',...
                     'visualInspection',...
                     'BipolarReferencing'...      
                     'downsample'...
            };

        elseif strcmp(ops.order,'visuallyInspectFromRaw')
            order = {'highpassFilter',... 
                     'notchFilter',...
                     'GlobalMeanRemoval',...
                     'visualInspection',...
            };

        elseif strcmp(ops.order,'oldPipeline') % as close as we can get
            order = {'visualInspection',...
                     'highpassFilter',... 
                     'CAR',...
                     'notchFilter',...
                     'GaussianFilterExtraction',...
                     'downsample'...
            };
        
        elseif strcmp(ops.order,'defaultECOGBroadBand')
            order = {'highpassFilter',... 
                        'notchFilter',...
                        'IEDRemoval',...
                        'visualInspection',...
                        'CAR',...                        
            };
        elseif strcmp(ops.order,'defaultSEEGorBOTHBroadBandNoIED')
            order = {'highpassFilter',... 
                        'notchFilter',...
                        'CAR',... 
                        'LaplacianReferencing',...
                        'BipolarReferencing'...                              
            };

        elseif strcmp(ops.order,'defaultSEEGorBOTHBroadBand')
            order = {'highpassFilter',... 
                        'notchFilter',...
                        'IEDRemoval',...
                        'CAR',... 
                        'LaplacianReferencing',...
                        'BipolarReferencing'...                              
            };

        elseif strcmp(ops.order,'defaultSEEGorBOTHBroadBandNoReference')
            order = {'highpassFilter',... 
                        'notchFilter',...
                        'IEDRemoval',... 
                        'GlobalMeanRemoval',...
            };

        elseif strcmp(ops.order,'test')
            order = {'downsample'};

        else
            error('Preprocessing order not specified');
        end

        obj.for_preproc.order = order;
        obj.for_preproc.isPlotVisible = ops.isPlotVisible; 

        % mappping from preprocessing name to its function in ecog_data
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
                 'NapLabFilterExtraction',...
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

        
        % -----------------------------
        % CALL PREPROCESSING FUNCTIONS
        % -----------------------------
        fprintf(1,'\nSTARTING TO PREPROCESS SIGNAL\n');

        obj.first_step('doneVisualInspection',ops.doneVisualInspection);

        prev_step = '';
        for i=1:length(function_order)
            step = function_order{i};

            % skip current step if same as last step 
            if strcmp(step,prev_step)
                continue
            end

            % add flags to reference_signal function
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
                flags = char(flags); flags = flags(1:end-1); % remove comma at end
                eval(strcat("obj.",step,"(",flags,");"));

            % add flags to extract_high_gamma function
            elseif strcmp(step,'extract_high_gamma')
                flags = strcat("'do",order{i},"',true");
                eval(strcat("obj.",step,"(",flags,");"));

            % visual inspection
            elseif strcmp(step,'visual_inspection')
                eval(strcat("obj.",step,"('doneVisualInspection',",num2str(ops.doneVisualInspection),");"));
                
            % functions without input flags/arguments
            else
                eval(strcat("obj.",step,"();"));
            end

            prev_step = step;

        end

        fprintf(1,'\nDONE PREPROCESSING SIGNAL \n');

    end


    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIRST STEP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj] = first_step(obj,varargin)
        % Reset preprocessed parameters to raw values and re-define 
        % parameters (e.g., filter params).
        % 
        % Returns : obj

        p = inputParser();
        addParameter(p,'doneVisualInspection',false);
        parse(p, varargin{:});
        ops = p.Results;

        obj.stitch_index = obj.for_preproc.stitch_index_raw;
        obj.sample_freq  = obj.for_preproc.sample_freq_raw;

        % re-define parameters
        obj.define_parameters();

        % % check if visual inspection needs to be done again
        % if isfield(obj.for_preproc,'visual_inspection') % if it has been done before
        %     if any(cellfun(@(x) strcmp(x,'visualInspection'),obj.for_preproc.order,'UniformOutput',false))
        %         obj.for_preproc.visual_inspection = []; % only clear if visual inspection will be done again
        %         obj.elec_ch_user_deselect = [];
        %     else % use old visual inspection 
        %         fprintf(1,'\nChannels to be removed because of previously completed visual inspection: ');
        %         fprintf(1,'%d ',obj.elec_ch_user_deselect(:)); fprintf(1,'\n');
        %     end
        % end

        % reset noisy electrodes
        obj.elec_ch_with_IED = [];
        obj.elec_ch_with_noise = [];
        if ~ops.doneVisualInspection
            obj.elec_ch_user_deselect = [];
        end

        % define clean electrodes
        obj.define_clean_channels()

        % clear results from preprocessing steps
        obj.for_preproc.notchFilter_results = [];
        obj.for_preproc.IEDRemoval_results = [];
        obj.for_preproc.visualInspection_results = [];
        obj.for_preproc.outlierRemoval_results = [];

        % reset all bipolar info
        obj.bip_elec_data    = [];
        obj.bip_ch           = [];
        obj.bip_ch_label     = [];
        obj.bip_ch_valid     = [];
        obj.bip_ch_grp       = [];
        obj.bip_ch_label_grp = [];

        % update trial timing table to correspond to the original signal 
        if isfield(obj.for_preproc,'trial_timing_raw') % might not exist if first time preproc
            obj.trial_timing = obj.for_preproc.trial_timing_raw;
        end

        obj.elec_data = obj.for_preproc.elec_data_raw;

    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function define_parameters(obj)
        % Define parameters to be used in preprocessing and store them 
        % in obj.for_preproc. This is the ONLY function where parameters 
        % should be modified by the user.
        % 
        % Returns : obj


        % SET FILTER PARAMETERS

        param = struct;
            
        % --- highpass filter ---
        param.highpass.Wp = 0.50; % Hz
        param.highpass.Ws = 0.05; % Hz
        param.highpass.Rp = 3;    % dB
        param.highpass.Rs = 30;   % dB
            
        % --- IIR peak filter ---
        param.peak.fcenter = [55,60,65];
        param.peak.bw      = ones(1,length(param.peak.fcenter)).*0.001;
            
        % --- notch filter ---
        param.notch.fcenter = [60,120,180,240];
        param.notch.bw = ones(1,length(param.notch.fcenter)).*0.001;

        % --- gaussian filter ---
        param.gaussian.f_gamma_low = 70;
        param.gaussian.f_gamma_high = 150;

        % --- bandpass filter ---
        param.bandpass.f_gamma_low = 70;
        param.bandpass.f_gamma_high = 150;
        param.bandpass.filter_order = 6;


        % CONSTRUCT FILTERS (shouldn't need to edit)

        % --- highpass filter ---
        highpass.Wp = param.highpass.Wp/(obj.sample_freq/2); 
        highpass.Ws = param.highpass.Ws/(obj.sample_freq/2);
        highpass.Rp = param.highpass.Rp; 
        highpass.Rs = param.highpass.Rs;

        [highpass.n,highpass.Wn] = buttord(highpass.Wp,highpass.Ws,highpass.Rp,highpass.Rs);
        highpass.n = highpass.n + rem(highpass.n,2);

        % caclulate the filter coefficients in Zero-Pole-Gain design
        [highpass.z,highpass.p,highpass.k] = butter(highpass.n,highpass.Wn,'high');
        [highpass.sos,highpass.g] = zp2sos(highpass.z,highpass.p,highpass.k);
        highpass.h = dfilt.df2sos(highpass.sos,highpass.g);

        % --- IIR peak filter ---
        % calculate the IIR-peak filter coefficients in a,b format
        for idx = 1:length(param.peak.fcenter)
            peak{idx}.wo = param.peak.fcenter(idx)/(obj.sample_freq/2);
            peak{idx}.bw = param.peak.bw(idx);  
            [peak{idx}.b,peak{idx}.a] = iirpeak(peak{idx}.wo,peak{idx}.bw);
        end

        % --- notch filter ---
        % calculate the IIR-peak filter coefficients in a,b format 
        for idx = 1:length(param.notch.fcenter)
            notch{idx}.wo = param.notch.fcenter(idx)/(obj.sample_freq/2);  
            notch{idx}.bw = param.notch.bw(idx);
            [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
        end 

        % --- gaussian filter ---
        %
        % See:
        %   Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
        %   The Control of Vocal Pitch in Human Laryngeal Motor Cortex. Cell 174 (1).  
        [gaussian.cfs,gaussian.sds] = obj.get_filter_param_chang_lab(param.gaussian.f_gamma_low,param.gaussian.f_gamma_high);
        gaussian.cfs(1) = 73.0; % correction to stay away from 60 hz
        gaussian.f_gamma_low = param.gaussian.f_gamma_low;
        gaussian.f_gamma_high = param.gaussian.f_gamma_high;

        % --- bandpass filter ---
        % Construct an FDESIGN object and call its BUTTER method.
        bandpass.h = fdesign.bandpass('N,F3dB1,F3dB2',param.bandpass.filter_order,param.bandpass.f_gamma_low,param.bandpass.f_gamma_high,obj.sample_freq);
        bandpass.Hd = design(bandpass.h,'butter');
        [bandpass.B, bandpass.A] = sos2tf(bandpass.Hd.sosMatrix,bandpass.Hd.scaleValues);
        bandpass.filter_order = param.bandpass.filter_order;
        bandpass.f_gamma_low = param.bandpass.f_gamma_low;
        bandpass.f_gamma_high = param.bandpass.f_gamma_high;


        % OTHER PARAMETERS

        % zero out signal
        zero = 1; % seconds

        % outlier removal
        outlier.trimmed = 1; % second
        outlier.threshold = 5;
        outlier.percentile = 0.9;
        outlier.buffer = 20;
        outlier.interpMethod = 'linear';


         % add fields to obj.info
        % WILL OVERWRITE OLD FILTER PARAMETERS
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
    % OUTPUT DATA STRUCTURES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function output_data_structures(obj)
        % Outputs a MATLAB structure and a python Xarray of object.
        %
        %

    end 


    
    

    
end


end

