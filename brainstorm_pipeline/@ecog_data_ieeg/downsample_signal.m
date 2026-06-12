function downsample_signal(obj,ops)
        % Downsamples signal and all other relevant parts of object
        % 
        arguments
            obj ecog_data_ieeg
            ops.decimationFreq = obj.for_preproc.decimation_freq; 
        end
        % p = inputParser();
        % addParameter(p,'decimationFreq',obj.for_preproc.decimation_freq)
        % parse(p, varargin{:});
        % ops = p.Results;
        
        signal = obj.elec_data';
        if ~isempty(obj.bip_elec_data)
            signal_bipolar = obj.bip_elec_data';
        end

        % assert(obj.sample_freq > obj.for_preproc.decimation_freq,'signal has already been downsampled');

        signal_dec = [];
        signal_bipolar_dec = [];

        curr_stitch = 1;
        stitch_index = []; % first sample in each data file

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Resampling signal from file %d of %d ... \n',k,length(obj.stitch_index));

            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end

            % --- UNIPOLAR ---
            signal_ = signal(obj.stitch_index(k):stop,:);
            signal_ = resample(double(signal_),ops.decimationFreq,obj.sample_freq);
            signal_dec = [signal_dec; signal_];
        
            % --- BIPOLAR ---
            if ~isempty(obj.bip_elec_data)
                signal_bipolar_ = signal_bipolar(obj.stitch_index(k):stop,:);
                signal_bipolar_ = resample(double(signal_bipolar_),ops.decimationFreq,obj.sample_freq);
                signal_bipolar_dec = [signal_bipolar_dec; signal_bipolar_];
            end

            % update stitch index
            stitch_index = [stitch_index; curr_stitch];
            curr_stitch = curr_stitch + size(signal_,1);

        end

        obj.elec_data = signal_dec';

        if ~isempty(obj.bip_elec_data)
            obj.bip_elec_data = signal_bipolar_dec';
        end

        % construct new trial timing table if downsample rate is not the preset downsample rate
        % (NOT ADVISED)
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

        % set trial timing table to the downsampled version if exists
        if ~isempty(obj.trial_timing) && (ops.decimationFreq == obj.for_preproc.decimation_freq)
            obj.trial_timing = obj.for_preproc.trial_timing_dec;
        end

        obj.stitch_index = stitch_index;
        obj.sample_freq = ops.decimationFreq;

    end
