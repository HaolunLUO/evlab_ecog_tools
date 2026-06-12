%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE COMMMON NOISE AND REFERENCE SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reference_signal(obj,ops)
        % Should be run after visually_inspect
        % TODO - decide what to do with the EEG channels

        arguments
            obj ecog_data;
            ops.doGlobalMeanRemoval = false;
            ops.doCAR = true;
            ops.doShankCSR = false;
            ops.doBipolarReferencing = false;
            ops.doLaplacianReferencing = false;
        end
        
        % p = inputParser();
        % addParameter(p,'doGlobalMeanRemoval',false)
        % addParameter(p,'doCAR',true);
        % addParameter(p,'doShankCSR',false);
        % addParameter(p,'doBipolarReferencing',true);
        % parse(p, varargin{:});
        % ops = p.Results;

        signal = obj.elec_data';
        signal_for_bip = obj.elec_data';

        % split ECoG (grids+srips) and SEEG channels
        ecog_chans = find(strcmp(obj.elec_ch_type, 'ecog_grid') | strcmp(obj.elec_ch_type, 'ecog_strip'));
        seeg_chans = find(strcmp(obj.elec_ch_type, 'seeg'));

        % initialize bipolar signal
        signal_bipolar = zeros(size(signal)); % might not be used

        for k=1:length(obj.stitch_index) % number of separate data files with signal
            fprintf(1, '\n> Referencing signal from file %d of %d ... \n',k,length(obj.stitch_index));
            
            if k == length(obj.stitch_index) % signal for file stops at end of matrix
                stop = size(signal,1);
            else % signal for file stops before stitch index of next file
                stop = obj.stitch_index(k+1)-1;
            end


            % --------------------
            % GLOBAL MEAN REMOVAL
            % --------------------
            if ops.doGlobalMeanRemoval
                fprintf(1, '\n>> Removing global mean of signal \n');

                signal_ = signal(obj.stitch_index(k):stop,:);

                overall_mean = mean(signal_(:,obj.elec_ch_clean),2);
                signal_ = signal_ - repmat(overall_mean,1,size(signal_,2));

                signal(obj.stitch_index(k):stop,:) = signal_;

            end


            % ---------------------------------------
            % COMMON AVERAGE REFERENCING (ECoG only) <-- NOT ANYMORE
            % ---------------------------------------
            if ops.doCAR && ~isempty(ecog_chans)
                fprintf(1, '\n>> Common average filtering signal \n');
                fprintf(1,'[');

                signal_ = signal(obj.stitch_index(k):stop,:);

                % determine the number of amps
                num_amps = ceil(size(signal_,2) / obj.for_preproc.elecs_per_amp);

                % exclude the channels that were noisy or not ECoG
                eligible_channels_ecog = intersect(ecog_chans,obj.elec_ch_clean);
                % eligible_channels = obj.elec_ch_clean;

                % for each of these amps
                for idx_amp = 1:num_amps
                    idx_low  = (idx_amp-1)*obj.for_preproc.elecs_per_amp+1;
                    idx_high = min((idx_amp-0)*obj.for_preproc.elecs_per_amp+0, max(eligible_channels_ecog));
        
                    % exclude the channels that were not on this amp
                    list_channels = intersect(eligible_channels_ecog,idx_low:idx_high);
        
                    % check if any channels are left
                    if ~isempty(list_channels) && length(list_channels)>1
                        % calculate the common average reference signal
                        signal_mean = mean(signal_(:,list_channels),2);
                        % subtract the common average signal from each channel of amp
                        for idx_ch = list_channels
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    else
                        % don't change the signal at all
                        for idx_ch = list_channels
                            signal_(:,idx_ch) = signal_(:,idx_ch);
                            fprintf(1,'.');
                        end
                    end
                end

                fprintf(1,'] done\n');

                signal(obj.stitch_index(k):stop,:) = signal_;

            % elseif ops.doCAR
            %     error('No ECoG channels to perform CAR on')
            end


            % ----------------------------------------
            % SHANK COMMON SOURCE REMOVAL (SEEG only)
            % ----------------------------------------
            if ops.doShankCSR && ~isempty(seeg_chans)
                fprintf(1, '\n>> Shank common source removal \n');
                fprintf(1,'[');

                signal_ = signal(obj.stitch_index(k):stop,:);

                % exclude the channels that were noisy or not SEEG
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);

                if ops.doCAR & ~isempty(ecog_chans)
                    assert(isempty(intersect(eligible_channels_ecog,eligible_channels_seeg)),"Some channels are labelled as both ECoG and SEEG.")
                end

                % split labels into shanks
                [shank_locs,~,~] = obj.extract_shanks();
                
                % for each shank
                for idx_shk = 1:size(shank_locs,1);

                    % exlude channels that were not on this shank
                    same_shank = shank_locs{idx_shk};
                    same_shank = intersect(eligible_channels_seeg,same_shank);

                    % check if any channels are left
                    if ~isempty(same_shank) && length(same_shank)>1
                        % calculate the common average reference signal
                        signal_mean = mean(signal_(:,same_shank),1);
                        % subtract the common average signal from each channel of shank
                        for idx_ch = same_shank
                            signal_(:,idx_ch) = signal_(:,idx_ch) - signal_mean;
                            fprintf(1,'.');
                        end
                    else
                        % don't change the signal at all 
                        % (either none clean, only one clean, or none SEEG)
                        for idx_ch = same_shank
                            signal_(:,idx_ch) = signal_(:,idx_ch);
                            fprintf(1,'.');
                        end
                    end
                end

                signal(obj.stitch_index(k):stop,:) = signal_;

                fprintf(1,'] done\n');

            elseif ops.doShankCSR 
                error('No SEEG channels to perform shank CSR on')
            end

            % --------------------------------
            % LAPLACIAN REFERENCING (SEEG only)
            % --------------------------------
            if ops.doLaplacianReferencing && ~isempty(seeg_chans)
                fprintf(1, '\n>> Performing Laplacian referencing for SEEG \n');
                fprintf(1,'[');
    
                signal_ = signal(obj.stitch_index(k):stop,:);
    
                % Get SEEG shanks
                [shank_locs, shank_tags, ~] = obj.extract_shanks();
                signal_laplace = signal_;  % channels on skipped shanks keep original values
                % Perform Laplacian referencing for each shank
                for shank_idx = 1:length(shank_locs)
                    shank_channels = shank_locs{shank_idx};
                    shank_channels = intersect(shank_channels, obj.elec_ch_clean);  % Only use clean channels
                        
                    if length(shank_channels) >= 2  % Need at least 2 channels for Laplacian
                        for i = 1:length(shank_channels)
                            curr_channel = shank_channels(i);
                            fprintf(1, ['\n>> Processing channel ' obj.elec_ch_label{curr_channel}]);
                            if i == 1  % First channel
                                next_channel = shank_channels(i+1);
                                laplacian = signal_(:,curr_channel) - signal_(:,next_channel);
                                fprintf(1, ['\n>> Using next channel ' obj.elec_ch_label{next_channel}]);
                            elseif i == length(shank_channels)  % Last channel
                                prev_channel = shank_channels(i-1);
                                laplacian = signal_(:,curr_channel) - signal_(:,prev_channel);
                                fprintf(1, ['\n>> Using previous channel ' obj.elec_ch_label{prev_channel}]);
                            else  % Middle channels
                                prev_channel = shank_channels(i-1);
                                next_channel = shank_channels(i+1);
                                fprintf(1, ['\n>> Using previous channel ' obj.elec_ch_label{prev_channel}]);
                                fprintf(1, ['\n>> Using next channel ' obj.elec_ch_label{next_channel}]);
                                laplacian = signal_(:,curr_channel) - ...
                                            0.5 * (signal_(:,prev_channel) + signal_(:,next_channel));
                            end
                            
                            signal_laplace(:,curr_channel) = laplacian;
                            fprintf(1,'.');
                        end
                    else
                        fprintf(1, '\nWarning: Shank %d has fewer than 2 clean channels. Skipping Laplacian referencing for this shank.\n', shank_idx);
                    end
                end
    
                signal(obj.stitch_index(k):stop,:) = signal_laplace;
    
                fprintf(1,'] done\n');
            elseif ops.doLaplacianReferencing
                error('No SEEG channels to perform Laplacian referencing on')
            end


            % --------------------------------
            % BIPOLAR REFERENCING (SEEG only)
            % --------------------------------
            if ops.doBipolarReferencing && ~isempty(seeg_chans)
                fprintf(1, '\n>> Bipolar referencing signal \n');
                fprintf(1,'[');
                
                signal_ = signal_for_bip(obj.stitch_index(k):stop,:);

                % TODO - take the unneccesary bipolar stuff out of the loop

                % exclude the channels that were noisy or not SEEG
                eligible_channels_seeg = intersect(seeg_chans,obj.elec_ch_clean);
                
                if ops.doCAR & ~isempty(ecog_chans)
                    assert(isempty(intersect(eligible_channels_ecog,eligible_channels_seeg)),"Some channels are labelled as both ECoG and SEEG.")
                end

                % split labels into shanks 
                [shank_locs,~,shank_nums] = obj.extract_shanks();

                % find bipolar pairs of clean, SEEG electrodes
                % (i.e., electrodes on the same SEEG shank that are one apart and clean)
                chan_idx_for_bip = cellfun(@(x) intersect(x,eligible_channels_seeg,'stable'),shank_locs,'uni',false);
                chan_num_for_bip = cellfun(@(x) shank_nums(x),chan_idx_for_bip,'uni',false);
                % as cells 
                bipolar_diffs_idx_grp = cellfun(@(x,y) [y(find(diff(x)==1 & diff(y)==1))+1,y(find(diff(x)==1 & diff(y)==1))],chan_num_for_bip,chan_idx_for_bip,'uni',false)';
                bipolar_diffs_name_grp = cellfun(@(x) obj.elec_ch_label(x),bipolar_diffs_idx_grp,'uni',false);
                % as mats
                bipolar_diffs_idx = cell2mat(bipolar_diffs_idx_grp');
                bipolar_diffs_name = obj.elec_ch_label(bipolar_diffs_idx);
                
                bipolar_idxs = 1:size(bipolar_diffs_idx,1);
                bipolar_valid = ones(size(bipolar_diffs_idx,1),1);

                % create a biopolar version of the signal
                signal_bipolar_= double([]); 
                for bipolar_id = 1:size(bipolar_diffs_idx,1)
                    bipol_ch_1 = bipolar_diffs_idx(bipolar_id,1);
                    bipol_ch_2 = bipolar_diffs_idx(bipolar_id,2);
                    bip_ch_name1 = bipolar_diffs_name{bipolar_id,1};
                    bip_ch_name2 = bipolar_diffs_name{bipolar_id,2};
              
                    % basic assertions
                    try
                        B = cell2mat(extract(bip_ch_name2,lettersPattern));
                        A =cell2mat(extract(bip_ch_name1,lettersPattern));
                    catch
                        B = regexprep(bip_ch_name2, '\d+(?:_(?=\d))?', '');
                        A = regexprep(bip_ch_name1, '\d+(?:_(?=\d))?', '');
                    end
                    assert(all(A==B),"Some channels are not on the same shank");
                    B =str2num(cell2mat(extract(bip_ch_name2,digitsPattern)));
                    A =str2num(cell2mat(extract(bip_ch_name1,digitsPattern)));
                    assert(A-B==1, "Some channels are more than one apart");
                    
                    % subtract channel 2 from channel 1 in pair
                    signal_bipolar_(:,bipolar_id) = signal_(:,bipol_ch_1)-signal_(:,bipol_ch_2);
                    fprintf(1,'.');
                end

                bipolar_diffs_name = arrayfun(@(x) {[bipolar_diffs_name{x,1} '-' bipolar_diffs_name{x,2}]},[1:size(bipolar_diffs_name,1)])';

                signal_bipolar(obj.stitch_index(k):stop,1:size(bipolar_diffs_idx,1)) = signal_bipolar_;

                fprintf(1,'] done\n');

            elseif ops.doBipolarReferencing
                error('No SEEG channels to perform bipolar referencing on')
            end
        end

        
        % unipolar 
        if ops.doGlobalMeanRemoval || ops.doCAR || ops.doShankCSR || ops.doLaplacianReferencing
            obj.elec_data = signal';
        end

        % bipolar
        if ops.doBipolarReferencing 

            % remove unused rows in bipolar signal
            signal_bipolar = signal_bipolar(:,1:size(bipolar_diffs_idx,1));

            obj.bip_elec_data    = signal_bipolar';
            obj.bip_ch           = bipolar_idxs';
            obj.bip_ch_label     = bipolar_diffs_name; 
            obj.bip_ch_valid     = bipolar_valid;
            obj.bip_ch_grp       = bipolar_diffs_idx_grp';
            obj.bip_ch_label_grp = bipolar_diffs_name_grp';
              
            if obj.for_preproc.isPlotVisible
                % plot to make sure signal looks right
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