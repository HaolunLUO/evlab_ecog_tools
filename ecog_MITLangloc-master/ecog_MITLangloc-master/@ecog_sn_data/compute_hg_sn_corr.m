function compute_hg_sn_corr(obj, varargin)
    % Compute S-vs-N Spearman correlation per electrode using EVEN trials only
    % (Odd trials are reserved for inference in test_s_vs_n)
    % Returns one correlation value per electrode stored in obj.hg_sn_corr

    p = inputParser();
    addParameter(p, 'words',             1:12);
    addParameter(p, 'S_condition_flag',  'S');
    addParameter(p, 'N_condition_flag',  'N');
    addParameter(p, 'sessions',          []);
    addParameter(p, 'corr_type',         'Spearman');
    parse(p, varargin{:});
    ops = p.Results;

    % --- Session / trial filtering ---
    if ~isempty(ops.sessions)
        keep_trials = any(obj.session == ops.sessions, 2);
    else
        keep_trials = [];
    end

    % --- Fetch averaged condition tables ---
    [~, S_ave_cond_table] = obj.get_ave_cond_trial( ...
        'words',      ops.words, ...
        'condition',  ops.S_condition_flag, ...
        'keep_trials', keep_trials);
    [~, N_ave_cond_table] = obj.get_ave_cond_trial( ...
        'words',      ops.words, ...
        'condition',  ops.N_condition_flag, ...
        'keep_trials', keep_trials);

    for elec_type = 1:2   % 1 = unipolar, 2 = bipolar

        if elec_type == 1
            S_data = S_ave_cond_table.elec_data{1};   % [nElec x nTrials x nWords]
            N_data = N_ave_cond_table.elec_data{1};
            strval = 'Unipolar';
        else
            if isempty(S_ave_cond_table.bip_elec_data{1})
                fprintf('\nNo bipolar data found, skipping bipolar.\n');
                continue;
            end
            S_data = S_ave_cond_table.bip_elec_data{1};
            N_data = N_ave_cond_table.bip_elec_data{1};
            strval = 'Bipolar';
        end

        % --- Keep EVEN trials only (dim 2 = trials) ---
        % Odd trials are reserved for test_s_vs_n inference
        S_data = S_data(:, 2:2:size(S_data, 2), :);
        N_data = N_data(:, 2:2:size(N_data, 2), :);

        nElecs    = size(S_data, 1);
        nS_trials = size(S_data, 2);
        nN_trials = size(N_data, 2);

        fprintf('\nComputing S-vs-N %s correlation for %s electrodes (n=%d) ...\n', ...
            ops.corr_type, strval, nElecs);

        % --- Average across words (dim 3): [nElec x nTrials] ---
        S_mat = mean(S_data, 3);
        N_mat = mean(N_data, 3);

        % --- Concatenate S and N trials, build condition flag ---
        combined = [S_mat, N_mat];                              % [nElec x (nS+nN)]
        flag     = [ones(1, nS_trials), -ones(1, nN_trials)];  % +1=S, -1=N

        % --- Spearman correlation: one value per electrode ---
        obs_corr = zeros(nElecs, 1);
        for e = 1:nElecs
            obs_corr(e) = corr(combined(e,:)', flag', 'Type', ops.corr_type);
        end

        % --- Store results ---
        res.corr = obs_corr;   % [nElec x 1]

        if elec_type == 1
            obj.hg_sn_corr.results.unipolar = res;
        else
            obj.hg_sn_corr.results.bipolar  = res;
        end

        fprintf('Completed S-vs-N correlation for %s (n=%d electrodes)\n', strval, nElecs);
    end

    obj.hg_sn_corr.ops = ops;
    fprintf('\nFinished compute_hg_sn_corr.\n');
end
