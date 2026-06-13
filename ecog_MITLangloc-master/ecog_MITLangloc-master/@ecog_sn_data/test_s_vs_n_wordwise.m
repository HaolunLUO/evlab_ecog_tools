function test_s_vs_n_wordwise(obj, varargin)
        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'S_condition_flag','S');
        addParameter(p,'N_condition_flag','N');
        addParameter(p,'n_rep',1000);
        addParameter(p,'corr_type','Spearman');
        addParameter(p,'threshold',0.01);
        addParameter(p,'side','both');
        addParameter(p,'sessions',[]);
        addParameter(p,'min_consecutive',3);
        addParameter(p,'consecutiveness',true);
        addParameter(p,'use_odd_for_inference',true); % NEW
        parse(p, varargin{:});
        ops = p.Results;
        
        [~, S_ave_cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition', ...
            ops.S_condition_flag,'keep_trials',ops.sessions);
        [~, N_ave_cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition', ...
            ops.N_condition_flag,'keep_trials',ops.sessions);
        
        for elec_type = 1:2
            if elec_type == 1
                S_data = S_ave_cond_table.elec_data{1};
                N_data = N_ave_cond_table.elec_data{1};
            else
                S_data = S_ave_cond_table.bip_elec_data{1};
                N_data = N_ave_cond_table.bip_elec_data{1};
            end
        
            % --- ODD/EVEN SPLIT (dim 2 = trials) ---
            if ops.use_odd_for_inference
                S_data = S_data(:, 1:2:size(S_data,2), :);
                N_data = N_data(:, 1:2:size(N_data,2), :);
            end
        

        nElecs = size(S_data, 1);
        nWords = size(S_data, 3);
        sn_corr = zeros(nElecs, nWords);
        p_value = zeros(nElecs, nWords);

        if(elec_type==1); strval = 'Unipolar'; else strval = 'Bipolar'; end
        fprintf('\nPerforming word-level langloc analysis for %s\n', strval);
        fprintf('[');
        for w = 1:nWords
            S_word = S_data(:,:,w); % elec x S_trials
            N_word = N_data(:,:,w); % elec x N_trials
            combined_word = [S_word, N_word];
            labels = [ones(1, size(S_word,2)), -ones(1, size(N_word,2))]';
            parfor elec = 1:nElecs
                data = combined_word(elec,:)';
                sn_corr(elec,w) = corr(data, labels, 'Type', ops.corr_type);
                % Permutation test
                corr_rand = zeros(1, ops.n_rep);
                for r = 1:ops.n_rep
                    perm_labels = labels(randperm(length(labels)));
                    corr_rand(r) = corr(data, perm_labels, 'Type', ops.corr_type);
                end
                switch ops.side
                    case 'left'
                        p_value(elec,w) = mean(corr_rand < sn_corr(elec,w));
                    case 'right'
                        p_value(elec,w) = mean(corr_rand > sn_corr(elec,w));
                    case 'both'
                        p_value(elec,w) = mean(abs(corr_rand) > abs(sn_corr(elec,w)));
                end
            end
            fprintf('.');
        end
        fprintf('] done\n');

        langloc = false(nElecs, 1);
        is_sig = p_value < ops.threshold;
        for elec = 1:nElecs
            sig_vec = is_sig(elec,:);
            if ops.consecutiveness
                % require min_consecutive significant RUN
                if any(conv(double(sig_vec), ones(1, ops.min_consecutive), 'valid') == ops.min_consecutive)
                    langloc(elec) = true;
                end
            else
                % require min_consecutive anywhere in the list
                if sum(sig_vec) >= ops.min_consecutive
                    langloc(elec) = true;
                end
            end
        end

        if elec_type == 1
            obj.langloc_wordwise.results.unipolar.sn_corr = sn_corr;
            obj.langloc_wordwise.results.unipolar.p_value = p_value;
            obj.langloc_wordwise.results.unipolar.is_sig = langloc;
            obj.langloc_wordwise.results.unipolar.is_sig_vec = is_sig;
        else
            obj.langloc_wordwise.results.bipolar.sn_corr = sn_corr;
            obj.langloc_wordwise.results.bipolar.p_value = p_value;
            obj.langloc_wordwise.results.bipolar.is_sig = langloc;
            obj.langloc_wordwise.results.bipolar.is_sig_vec = is_sig;
        end
    end
    obj.langloc_wordwise.ops = ops;
    fprintf('\nFinished wordwise analysis (consecutiveness=%d).\n', ops.consecutiveness);
end
