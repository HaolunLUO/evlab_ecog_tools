function compute_hg_power_diff_s_vs_n(obj, varargin)
    % Compute average difference in HG power between S and N conditions
    % Returns one value per electrode representing the mean HG power difference (S - N)
    % Data is already in HG band, no filtering required
    
    p = inputParser();
    addParameter(p, 'words', 1:12);
    addParameter(p, 'S_condition_flag', 'S');
    addParameter(p, 'N_condition_flag', 'N');
    addParameter(p, 'sessions', []);
    addParameter(p,'use_odd_for_inference',true); % NEW
    parse(p, varargin{:});
    ops = p.Results;
    
    % Get averaged data for both conditions
    [~, S_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.S_condition_flag, 'keep_trials', ops.sessions);
    [~, N_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.N_condition_flag, 'keep_trials', ops.sessions);
    
    for elec_type = 1:2   % 1=unipolar, 2=bipolar
        if elec_type == 1
            S_data = S_ave_cond_table.elec_data{1};
            N_data = N_ave_cond_table.elec_data{1};
        else
            S_data = S_ave_cond_table.bip_elec_data{1};
            N_data = N_ave_cond_table.bip_elec_data{1};
        end

        % --- ODD/EVEN SPLIT (dim 2 = trials) ---
        if ops.use_odd_for_inference
            S_data = S_data(:, 2:2:size(S_data,2), :);
            N_data = N_data(:, 2:2:size(N_data,2), :);
        end
        
        nElecs = size(S_data, 1);
        
        if(elec_type==1)
            strval = 'Unipolar';
        else
            strval = 'Bipolar';
        end
        fprintf('\nComputing HG power difference (S - N) for %s electrodes\n', strval);
        
        % Average across words and trials for each electrode
        % S_data and N_data are: elec x trial x words
        S_mean = mean(S_data, 3);      % Average across words: elec x trial
        S_mean = mean(S_mean, 2);      % Average across trials: elec x 1
        
        N_mean = mean(N_data, 3);      % Average across words: elec x trial
        N_mean = mean(N_mean, 2);      % Average across trials: elec x 1
        
        % Compute difference: S - N
        hg_power_diff = S_mean - N_mean;
        
        % Store results
        if elec_type == 1
            obj.hg_power_diff.results.unipolar.power_diff = hg_power_diff;
            obj.hg_power_diff.results.unipolar.power_S = S_mean;
            obj.hg_power_diff.results.unipolar.power_N = N_mean;
            obj.hg_power_diff.results.unipolar.power_S_word = squeeze(mean(S_data, 2));
            obj.hg_power_diff.results.unipolar.power_N_word = squeeze(mean(N_data, 2));
        else
            obj.hg_power_diff.results.bipolar.power_diff = hg_power_diff;
            obj.hg_power_diff.results.bipolar.power_S = S_mean;
            obj.hg_power_diff.results.bipolar.power_N = N_mean;
            obj.hg_power_diff.results.bipolar.power_S_word = squeeze(mean(S_data, 2));
            obj.hg_power_diff.results.bipolar.power_N_word = squeeze(mean(N_data, 2));
        end
        
        fprintf('Completed HG power analysis for %s (n=%d electrodes)\n', strval, nElecs);
    end
    
    obj.hg_power_diff.ops = ops;
    fprintf('\nFinished HG power difference computation.\n');
end