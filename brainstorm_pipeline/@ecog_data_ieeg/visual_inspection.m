 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VISUALLY INSPECT SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function visual_inspection(obj,ops)
        % Visually inspect electrodes without sig line noise
        % User will be asked to provide additional electrodes to remove
        % Should be run after highpass and notch filtering signal

        arguments
            obj ecog_data_ieeg
            ops.doneVisualInspection = false; 
        end

        % p = inputParser();
        % addParameter(p,'doneVisualInspection',false);
        % parse(p, varargin{:});
        % ops = p.Results;

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

        % only do visual inspection if hasn't been done already
        if ~ops.doneVisualInspection
            prompt1 = '\nUSER INPUT REQUIRED: \nAdditional channels to remove from visual inspection? (format: [1,2]) - ';
            prompt2 = 'Your name please :) - ';

            obj.elec_ch_user_deselect = input(prompt1)';
        
            % output structure for object
            vi_ops.inspected = 1;
            vi_ops.inspected_by = input(prompt2,'s');
            vi_ops.inspection_date = datestr(now, 'yyyy/mm/dd-HH:MM');

            obj.for_preproc.visualInspection_results = vi_ops;
        end
        
        % update clean/valid electrodes
        obj.define_clean_channels();

    end