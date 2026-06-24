classdef ecog_sn_data_seeg < ecog_data_seeg
% ECOG_SN_DATA_SEEG  Analysis class for sentence-vs-nonword (S vs N) localization.
%
% Extends ecog_data_seeg (which itself extends the advanced ieeg_pipeline engine
% @ecog_data_ieeg) with:
%   - test_s_vs_n()          configurable S/N permutation test
%   - lang_resp_plots()      timecourse + barplot pipeline
%   - get_timecourses()      extract trial-averaged timecourses per condition
%   - get_word_averages()    extract word-position averages per condition
%   - get_summary_statistics() tabular summary of localization results
%   - Plotting methods       plot_timecourse, plot_barplot, plot_s_vs_n
%
% Unlike the older MGH_utils/ecog_sn_data.m, condition flags (S_condition,
% N_condition, etc.) are configurable rather than hardcoded to 'S'/'N'.
%
% Constructor loads from a crunched .mat file written by
% brainstorm_to_mit_crunched_new. The variable in that file must be named
% 'obj' and must be an instance of ecog_data_seeg.

properties
    %% ---- LANG ELECS ----
    s_vs_n_sig              
    s_vs_n_p_ratio
    s_vs_n_corr
    s_vs_n_ops

    %% ---- EFFECT SIZE (held-out / even trials) ----
    hg_power_diff
    hg_sn_corr

    %% ---- LANG ELECS WORD LEVEL ANALYSIS ----
    langloc_wordwise

    %% ---- LANG ELECS WORD BOUNDARIES ANALYSIS (time-series) ----
    s_vs_n_wordboundaries_ops
    s_vs_n_wordboundaries_pSig
    s_vs_n_wordboundaries_pSig_bip
    s_vs_n_wordboundaries_wordBoundaries
    s_vs_n_wordboundaries_timePointsPerWord
    s_vs_n_wordboundaries_totalTimePoints
    s_vs_n_wordboundaries_sigUnipolarChannels
    s_vs_n_wordboundaries_sigBipolarChannels
    s_vs_n_wordboundaries_concatenatedEpochsSentence
    s_vs_n_wordboundaries_concatenatedEpochsNonword
    s_vs_n_wordboundaries_concatenatedEpochsSentence_bip
    s_vs_n_wordboundaries_concatenatedEpochsNonword_bip

    %% ---- PATHS ----
    langloc_save_path
    langloc_crunched_file_name
    langloc_crunched_file_path
    preproc_class_file_name
    preproc_class_file_path
end
    
methods
    %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function sn_obj = ecog_sn_data_seeg(...
            langloc_save_path,...
            langloc_crunched_file_name,...
            langloc_crunched_file_path,...
            preproc_class_file_name,...
            preproc_class_file_path)
    
        load(langloc_crunched_file_name) % loads variable 'obj'

        sn_obj@ecog_data_seeg(obj.for_preproc,...
                         obj.subject,...
                         obj.experiment,...
                         obj.crunched_file_name,...
                         obj.crunched_file_path,...
                         obj.raw_file_name,...
                         obj.raw_file_path,...
                         obj.elec_ch_label,...
                         obj.elec_ch,...
                         obj.elec_ch_prelim_deselect,...
                         obj.elec_ch_type...
        );

        sn_obj.langloc_save_path          = langloc_save_path;
        sn_obj.langloc_crunched_file_name = langloc_crunched_file_name;
        sn_obj.langloc_crunched_file_path = langloc_crunched_file_path;
        sn_obj.preproc_class_file_name    = preproc_class_file_name;
        sn_obj.preproc_class_file_path    = preproc_class_file_path;

        sn_obj.elec_data     = obj.elec_data;
        sn_obj.bip_elec_data = obj.bip_elec_data;
        sn_obj.stitch_index  = obj.stitch_index;
        sn_obj.sample_freq   = obj.sample_freq;
        sn_obj.trial_data    = obj.trial_data;

        sn_obj.trial_timing = obj.trial_timing;
        sn_obj.events_table = obj.events_table;
        sn_obj.condition    = obj.condition;
        sn_obj.session      = obj.session;
        
        sn_obj.elec_ch_with_IED      = obj.elec_ch_with_IED;
        sn_obj.elec_ch_with_noise    = obj.elec_ch_with_noise;
        sn_obj.elec_ch_user_deselect = obj.elec_ch_user_deselect;
        sn_obj.elec_ch_clean         = obj.elec_ch_clean;
        sn_obj.elec_ch_valid         = obj.elec_ch_valid;
        sn_obj.bip_ch                = obj.bip_ch;
        sn_obj.bip_ch_label          = obj.bip_ch_label;
        sn_obj.bip_ch_valid          = obj.bip_ch_valid;
        sn_obj.bip_ch_grp            = obj.bip_ch_grp;
        sn_obj.bip_ch_label_grp      = obj.bip_ch_label_grp;
    end 


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEST S VS N  (permutation correlation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function test_s_vs_n(obj,varargin)
        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'S_condition_flag','S')
        addParameter(p,'N_condition_flag','N')
        addParameter(p,'n_rep',1000);
        addParameter(p,'corr_type','Spearman');
        addParameter(p,'threshold',0.001);
        addParameter(p,'side','both');
        addParameter(p,'sessions',[]);
        addParameter(p,'do_plot',false);
        addParameter(p,'use_odd_for_inference',true); % odd trials for inference, even reserved for plotting/effect size
        parse(p, varargin{:});
        ops = p.Results;
        
        obj.s_vs_n_ops = ops;
            
        fprintf(1,'\n> Finding electrodes that respond significantly more to %s than %s ...\n',ops.S_condition_flag,ops.N_condition_flag);

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session==ops.sessions,2);
        else
            keep_trials = [];
        end

        [S_ave_tbl,S_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.S_condition_flag,'keep_trials',keep_trials);
        B = S_ave_tbl;
        S_ave = table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));

        [N_ave_tbl,N_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.N_condition_flag,'keep_trials',keep_trials);
        B = N_ave_tbl;
        N_ave = table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));

        % --- ODD/EVEN TRIAL SPLIT ---
        % Use odd trials for inference; even trials are reserved for plotting
        % (lang_resp_plots / get_word_averages) and held-out effect sizes
        % (compute_hg_power_diff_s_vs_n / compute_hg_sn_corr).
        if ops.use_odd_for_inference
            fprintf(1,'> Using ODD trials for inference (even trials reserved for plotting/effect size) ...\n');
            S_ave = cellfun(@(x) x(:,1:2:size(x,2)),S_ave,'uni',false);
            N_ave = cellfun(@(x) x(:,1:2:size(x,2)),N_ave,'uni',false);
        end

        S_N_comb = cellfun(@(x,y) horzcat(x,y),S_ave,N_ave,'uni',false);
        S_N_flag = cellfun(@(x,y) horzcat(x*0+1,y*0-1),S_ave,N_ave,'uni',false);
        S_N_corr = cellfun(@(x,y) diag(corr(x',y','Type',ops.corr_type)),S_N_comb,S_N_flag,'uni',false);
        names_no_string = ~ismember(B.Properties.VariableNames,{'string'});
        S_N_corr_table = cell2table(horzcat('s_vs_n_corr',S_N_corr),'VariableNames',B.Properties.VariableNames(names_no_string));

        fprintf(1,'\n> Running permutation analysis ...\n');
        fprintf(1,'[');

        n_rep = ops.n_rep;
        S_N_corr_rnd_all = cell(size(S_N_corr));

        for k=1:n_rep
            if ~mod(k,100) 
                fprintf(1,'.');
            end
            S_N_flag_rnd_idx = cellfun(@(x) randperm(size(x,2)),S_N_flag,'uni',false);
            S_N_flag_rnd = cellfun(@(x,y) x(:,y),S_N_flag,S_N_flag_rnd_idx,'uni',false);
            S_N_corr_rnd = cellfun(@(x,y) diag(corr(x',y','Type','Spearman')),S_N_comb,S_N_flag_rnd,'uni',false);
            S_N_corr_rnd_all = cellfun(@(x,y) horzcat(x,y),S_N_corr_rnd_all,S_N_corr_rnd,'UniformOutput',false);
        end

        fprintf(1,'] done\n');
            
        S_N_corr_rnd_table = cell2table(horzcat('s_vs_n_corr_rnd',S_N_corr_rnd_all),'VariableNames',B.Properties.VariableNames(names_no_string));
            
        switch ops.side
            case 'left'
                S_N_p_ratio = cellfun(@(x,y) sum(x<y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x) (x<(ops.threshold)),S_N_p_ratio,'uni',false);
            case 'right'
                S_N_p_ratio = cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x) (x>(1-ops.threshold)),S_N_p_ratio,'uni',false);
            case 'both'
                S_N_p_ratio = cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                S_N_p_is_sig = cellfun(@(x)  ((x<(ops.threshold)) | (x>(1-ops.threshold) )),S_N_p_ratio,'uni',false);
        end

        S_N_p_ratio_tbl = cell2table(horzcat('s_vs_n_p_ratio',S_N_p_ratio),'VariableNames',B.Properties.VariableNames(names_no_string));
            
        S_N_p_is_sig{1} = obj.elec_ch_valid & S_N_p_is_sig{1};
        s_vs_n_sig = cell2table(horzcat('s_vs_n_sig',S_N_p_is_sig),'VariableNames',B.Properties.VariableNames(names_no_string));
        obj.s_vs_n_sig = s_vs_n_sig;
        obj.s_vs_n_p_ratio = S_N_p_ratio_tbl;
        obj.s_vs_n_corr = S_N_corr;

        if ops.do_plot 
            obj.plot_s_vs_n(ops,...
                            'elec_data',...
                            obj.elec_ch_label,...
                            'S_table',S_table,...
                            'N_table',N_table,...
                            'S_N_corr_table',S_N_corr_table,...
                            'S_N_corr_rnd_table',S_N_corr_rnd_table,...
                            's_vs_n_sig',s_vs_n_sig,...
                            'S_N_p_ratio_tbl',S_N_p_ratio_tbl...
            );

            if ~isempty(obj.bip_elec_data)
                obj.plot_s_vs_n(ops,...
                                'bip_elec_data',...
                                obj.bip_ch_label,...
                                'S_table',S_table,...
                                'N_table',N_table,...
                                'S_N_corr_table',S_N_corr_table,...
                                'S_N_corr_rnd_table',S_N_corr_rnd_table,...
                                's_vs_n_sig',s_vs_n_sig,...
                                'S_N_p_ratio_tbl',S_N_p_ratio_tbl...
                );
            end
        end
    end

        
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LANG RESP PLOTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function lang_resp_plots(obj,varargin)
        p = inputParser();
        addParameter(p,'words',1:12);
        addParameter(p,'S_condition_flag','S');
        addParameter(p,'N_condition_flag','N');
        addParameter(p,'W_condition_flag',[]);
        addParameter(p,'J_condition_flag',[]);
        addParameter(p,'sessions',[]);
        addParameter(p,'subAverage',false);
        addParameter(p,'bipolarByShank',false);
        addParameter(p,'use_odd_for_inference',true); % plot even trials when odd were used for inference
        parse(p, varargin{:});
        ops = p.Results;

        if ops.subAverage
            obj.zscore_signal();
            obj.make_trials();
        end

        % Plot the held-out (even) trials when odd trials were used for inference
        if ops.use_odd_for_inference
            tc_split = 'even';
        else
            tc_split = [];
        end

        [S_values_uni,S_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                              'condition',ops.S_condition_flag,...
                                                              'sessions',ops.sessions,...
                                                              'allElecs',ops.bipolarByShank,'split',tc_split);
        [N_values_uni,N_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                              'condition',ops.N_condition_flag,...
                                                              'sessions',ops.sessions,'split',tc_split);
        if ops.W_condition_flag
            [W_values_uni,W_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.W_condition_flag,...
                                                                  'sessions',ops.sessions,'split',tc_split);
            [J_values_uni,J_values_uni_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.J_condition_flag,...
                                                                  'sessions',ops.sessions,'split',tc_split);
        end

        if ~isempty(obj.bip_elec_data)
            [S_values_bip,S_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.S_condition_flag,...
                                                                  'sessions',ops.sessions,...
                                                                  'signalType','bipolar',...
                                                                  'allElecs',ops.bipolarByShank,'split',tc_split);
            [N_values_bip,N_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                  'condition',ops.N_condition_flag,...
                                                                  'sessions',ops.sessions,...
                                                                  'signalType','bipolar',...
                                                                  'allElecs',ops.bipolarByShank,'split',tc_split);
            if ops.W_condition_flag
                [W_values_bip,W_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                      'condition',ops.W_condition_flag,...
                                                                      'sessions',ops.sessions,...
                                                                      'signalType','bipolar','split',tc_split);
                [J_values_bip,J_values_bip_sem] = obj.get_timecourses('words',ops.words,...
                                                                      'condition',ops.J_condition_flag,...
                                                                      'sessions',ops.sessions,...
                                                                      'signalType','bipolar','split',tc_split);
            end
        end

        if ~isempty(obj.bip_elec_data)
            S_values = [{S_values_uni},{S_values_bip}];
            S_values_sem = [{S_values_uni_sem},{S_values_bip_sem}];
            N_values = [{N_values_uni},{N_values_bip}];
            N_values_sem = [{N_values_uni_sem},{N_values_bip_sem}];
            W_values = []; W_values_sem = [];
            J_values = []; J_values_sem = [];
            if ops.W_condition_flag
                W_values = [{W_values_uni},{W_values_bip}];
                W_values_sem = [{W_values_uni_sem},{W_values_bip_sem}];
                J_values = [{J_values_uni},{J_values_bip}];
                J_values_sem = [{J_values_uni_sem},{J_values_bip_sem}];
            end
        else
            S_values = {S_values_uni};
            S_values_sem = {S_values_uni_sem};
            N_values = {N_values_uni};
            N_values_sem = {N_values_uni_sem};
            W_values = []; W_values_sem = [];
            J_values = []; J_values_sem = [];
            if ops.W_condition_flag
                W_values = {W_values_uni};
                W_values_sem = {W_values_uni_sem};
                J_values = {J_values_uni};
                J_values_sem = {J_values_uni_sem};
            end
        end

        clearvars S_values_uni S_values_bip S_values_uni_sem S_values_bip_sem
        clearvars N_values_uni N_values_bip N_values_uni_sem N_values_bip_sem
        clearvars W_values_uni W_values_bip W_values_uni_sem W_values_bip_sem
        clearvars J_values_uni J_values_bip J_values_uni_sem J_values_bip_sem
        
        [S_values_ave,S_values_ave_sem,nTrialsS] = obj.get_word_averages('words',ops.words,...
                                                            'condition',ops.S_condition_flag,...
                                                            'sessions',ops.sessions,...
                                                            'allElecs',ops.bipolarByShank,...
                                                            'use_odd_for_inference',ops.use_odd_for_inference);
        [N_values_ave,N_values_ave_sem,nTrialsN] = obj.get_word_averages('words',ops.words,...
                                                            'condition',ops.N_condition_flag,...
                                                            'sessions',ops.sessions,...
                                                            'allElecs',ops.bipolarByShank,...
                                                            'use_odd_for_inference',ops.use_odd_for_inference);
        W_values_ave = []; W_values_ave_sem = []; nTrialsW = [];
        J_values_ave = []; J_values_ave_sem = []; nTrialsJ = [];
        if ops.W_condition_flag
            [W_values_ave,W_values_ave_sem,nTrialsW] = obj.get_word_averages('words',ops.words,...
                                                                'condition',ops.W_condition_flag,...
                                                                'sessions',ops.sessions,...
                                                                'use_odd_for_inference',ops.use_odd_for_inference);
            [J_values_ave,J_values_ave_sem,nTrialsJ] = obj.get_word_averages('words',ops.words,...
                                                                'condition',ops.J_condition_flag,...
                                                                'sessions',ops.sessions,...
                                                                'use_odd_for_inference',ops.use_odd_for_inference);
        end

        data = struct;
        data.S_values         = S_values;
        data.W_values         = W_values;
        data.J_values         = J_values;
        data.N_values         = N_values;
        data.S_values_sem     = S_values_sem;
        data.W_values_sem     = W_values_sem;
        data.J_values_sem     = J_values_sem;
        data.N_values_sem     = N_values_sem;
        data.S_values_ave     = S_values_ave;
        data.W_values_ave     = W_values_ave;
        data.J_values_ave     = J_values_ave;
        data.N_values_ave     = N_values_ave;
        data.S_values_ave_sem = S_values_ave_sem;
        data.W_values_ave_sem = W_values_ave_sem;
        data.J_values_ave_sem = J_values_ave_sem;
        data.N_values_ave_sem = N_values_ave_sem;
        data.nTrialsS         = nTrialsS;
        data.nTrialsW         = nTrialsW;
        data.nTrialsJ         = nTrialsJ;
        data.nTrialsN         = nTrialsN;

        if ops.bipolarByShank
            obj.plot_timecourse_bip_by_shank(data,'words',ops.words);
        else
            obj.plot_timecourse(data,'subAverage',ops.subAverage,'words',ops.words);
            if ops.subAverage
                obj.plot_barplot(data,'subAverage',true,'words',ops.words);
                obj.plot_barplot(data,'subAverage',false,'words',ops.words);
            end
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET TIMECOURSES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values,values_sem]=get_timecourses(obj,varargin)
        p = inputParser();
        addParameter(p,'words',[]);
        addParameter(p,'condition',[]);
        addParameter(p,'signalType','unipolar');
        addParameter(p,'sessions',[]);
        addParameter(p,'useLangElecs',true);
        addParameter(p,'allElecs',false);
        addParameter(p,'split',[]);
        addParameter(p,'average',true);
        parse(p, varargin{:});
        ops = p.Results;
    
        if strcmp(ops.signalType,'unipolar')
            data = obj.elec_data;
            sig_elecs = logical(obj.s_vs_n_sig.elec_data{1,1});
            valid = obj.elec_ch_valid;
        elseif strcmp(ops.signalType,'bipolar')
            data = obj.bip_elec_data;
            sig_elecs = logical(obj.s_vs_n_sig.bip_elec_data{1,1});
            valid = obj.bip_ch_valid;
        end
    
        starts_n_stops = cellfun(@(x) trial_timing_word_sample_span(x, ops.words), ...
            obj.trial_timing, 'UniformOutput', false);
    
        if ops.useLangElecs && ~ops.allElecs
            values = cellfun(@(x) data(sig_elecs,x(1):x(2)),starts_n_stops,'UniformOutput',false); 
        elseif ops.allElecs
            values = cellfun(@(x) data(:,x(1):x(2)),starts_n_stops,'UniformOutput',false);
        else
            values = cellfun(@(x) data((~sig_elecs & valid),x(1):x(2)),starts_n_stops,'UniformOutput',false);
        end
    
        if ~isempty(ops.sessions)
            values = values(strcmp(obj.condition,ops.condition) & any(obj.session==ops.sessions,2));
        else
            values = values(strcmp(obj.condition,ops.condition));
        end
    
        if strcmp(ops.split,'odd')
            values = values(1:2:length(values)); 
        elseif strcmp(ops.split,'even')
            values = values(2:2:length(values));
        end
    
        if ~isempty(values)
            minLength = min(cellfun(@(x) size(x,2), values));
            values = cellfun(@(x) x(:, 1:minLength), values, 'UniformOutput', false);
        end
    
        if ops.average
            if ~isempty(values)
                values_done = zeros(size(values{1,1},1),minLength);
                values_sem = zeros(size(values{1,1},1),minLength);
                for i=1:size(values{1,1},1)
                    tmp = zeros(length(values),minLength);
                    for j=1:length(values)
                        tmp(j,:) = values{j,1}(i,1:minLength);
                    end
                    values_done(i,:) = mean(tmp,1);
                    values_sem(i,:) = std(tmp,1)/sqrt(size(tmp,1));
                end
            else
                values_done = [];
                values_sem = [];
            end
        else
            if ~isempty(values)
                values_done = zeros(length(values)*size(values{1,1},1),minLength);
                values_sem = [];
                idx_curr = 1;
                for i=1:size(values{1,1},1)
                    tmp = zeros(length(values),minLength);
                    for j=1:length(values)
                        tmp(j,:) = values{j,1}(i,1:minLength);
                    end
                    values_done(idx_curr:idx_curr+length(values)-1,:) = tmp(:,:);
                    idx_curr = idx_curr + length(values);
                end
            else
                values_done = [];
                values_sem = [];
            end
        end
        values = values_done;
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GET WORD AVERAGES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [values_ave,values_ave_sem,nTrials]=get_word_averages(obj,varargin)
        p = inputParser();
        addParameter(p,'words',[]);
        addParameter(p,'condition',[]);
        addParameter(p,'sessions',[]);
        addParameter(p,'useLangElecs',true);
        addParameter(p,'allElecs',false);
        addParameter(p,'use_odd_for_inference',true); % use even trials here if odd were used for inference
        parse(p, varargin{:});
        ops = p.Results;
    
        if ~isempty(ops.sessions)
            keep_trials = any(obj.session==ops.sessions,2);
        else
            keep_trials = [];
        end
        
        [~,cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition',ops.condition,'keep_trials',keep_trials);
        values_ave = table2cell(cond_table(:,~ismember(cond_table.Properties.VariableNames,{'key','string'})));
            
        if ~isempty(values_ave) && ~isempty(values_ave{1})
            minWordLength = min(cellfun(@(x) size(x,3), values_ave));
            values_ave = cellfun(@(x) x(:,:,1:minWordLength), values_ave, 'UniformOutput', false);
        end

        % EVEN trials for plotting (odd were reserved for inference in test_s_vs_n)
        if ops.use_odd_for_inference
            values_ave = cellfun(@(x) x(:,2:2:size(x,2),:), values_ave, 'UniformOutput', false);
        end
    
        nTrials = size(values_ave{1},2);
        values_ave_sem = cellfun(@(x) squeeze(std(x,[],2))/sqrt(nTrials),values_ave,'UniformOutput',false);
        values_ave = cellfun(@(x) squeeze(mean(x,2)),values_ave,'UniformOutput',false);
    
        sig_elecs = table2cell(obj.s_vs_n_sig(:,~ismember(obj.s_vs_n_sig.Properties.VariableNames,{'key'})));
        if ops.useLangElecs && ~ops.allElecs
            elecs_use = sig_elecs;
        elseif ops.allElecs  
            elecs_use = cellfun(@(x) logical(ones(size(x))),sig_elecs,'UniformOutput',false);
        else
            elecs_use = cellfun(@(x) ~x,sig_elecs,'UniformOutput',false);
        end
        values_ave = cellfun(@(x,y) x(y,:),values_ave,elecs_use,'UniformOutput',false);
        values_ave_sem = cellfun(@(x,y) x(y,:),values_ave_sem,elecs_use,'UniformOutput',false);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUMMARY STATISTICS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function summary=get_summary_statistics(obj,varargin)
        p = inputParser();
        addParameter(p,'sessions',[]);
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = find(any(obj.session==ops.sessions,2));
        else
            keep_trials = 1:length(obj.session);
        end

        columns = {'subject','num_sig','num_sig_left','num_sig_right',...
                   'num_clean','num_clean_left','num_clean_right',...
                   'num_prelim_deselect','num_w_IED','num_w_noise','num_user_deselect',...
                   'num_total_ecog_seeg','num_total_ecog_seeg_left','num_total_ecog_seeg_right',...
                   'num_total','percent_sig','chans_sig','names_sig','p_values_sig',...
                   'has_bipolar','num_sig_bipolar','num_total_bipolar','percent_sig_bipolar',...
                   'chans_sig_bipolar','names_sig_bipolar','p_values_sig_bipolar',...
                   'native_sample_freq','elecs_per_amp','runs_analyzed','num_runs_total',...
                   'num_words','presentation_rate','num_trials_S','num_trials_N',...
                   'num_trials_W','num_trials_J'};

        hasAnatomy = isfield(obj, 'anatomy') && isstruct(obj.anatomy) && ...
                     isfield(obj.anatomy, 'hemisphere') && isfield(obj.anatomy, 'mapping');

        if hasAnatomy
            mapped_hemi = cellfun(@(x) obj.anatomy.hemisphere(x),...
                                   obj.anatomy.mapping,'UniformOutput',false);
        else
            mapped_hemi = repmat({{'unknown'}}, length(obj.elec_ch_label), 1);
        end
        
        num_sig = {sum(obj.s_vs_n_sig.elec_data{1})};
        
        if hasAnatomy
            num_sig_left = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'left'),...
                            find(obj.s_vs_n_sig.elec_data{1}),'UniformOutput',false)))};
            num_sig_right = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'right'),...
                             find(obj.s_vs_n_sig.elec_data{1}),'UniformOutput',false)))};
        else
            num_sig_left = {[]};
            num_sig_right = {[]};
        end
        
        num_clean = {length(obj.elec_ch_clean)};
        
        if hasAnatomy
            num_clean_left = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'left'),...
                              obj.elec_ch_clean,'UniformOutput',false)))};
            num_clean_right = {sum(cell2mat(arrayfun(@(x) strcmp(mapped_hemi{x}{1},'right'),...
                               obj.elec_ch_clean,'UniformOutput',false)))};
        else
            num_clean_left = {[]};
            num_clean_right = {[]};
        end
        
        num_prelim_deselect = {length(obj.elec_ch_prelim_deselect)};
        num_w_IED = {length(obj.elec_ch_with_IED)};
        num_w_noise = {length(obj.elec_ch_with_noise)};
        num_user_deselect = {length(obj.elec_ch_user_deselect)};
        
        num_total_ecog_seeg = {sum(cell2mat(cellfun(@(x) (contains(x,'ecog') | contains(x,'seeg')),...
                                                  obj.elec_ch_type,'UniformOutput',false)))};

        if hasAnatomy
            num_total_ecog_seeg_left = {sum(cell2mat(cellfun(@(x,y) ((contains(x,'ecog') | contains(x,'seeg')) & ...
                                       strcmp(y,'left')),obj.elec_ch_type,mapped_hemi,'UniformOutput',false)))};
            num_total_ecog_seeg_right = {sum(cell2mat(cellfun(@(x,y) ((contains(x,'ecog') | contains(x,'seeg')) & ...
                                        strcmp(y,'right')),obj.elec_ch_type,mapped_hemi,'UniformOutput',false)))};
        else
            num_total_ecog_seeg_left = {[]};
            num_total_ecog_seeg_right = {[]};
        end
        
        num_total = {size(obj.elec_data,1)};
        percent_sig = {round((num_sig{1}/num_clean{1})*100,2)};
        chans_sig = {find(obj.s_vs_n_sig.elec_data{1})};
        names_sig = {obj.elec_ch_label(obj.s_vs_n_sig.elec_data{1})};
        ratio_tmp = obj.s_vs_n_p_ratio.elec_data{1};
        p_values_sig = {ratio_tmp(obj.s_vs_n_sig.elec_data{1})};

        has_bipolar = ~isempty(obj.bip_elec_data) && ...
                      ismember('bip_elec_data', obj.s_vs_n_sig.Properties.VariableNames) && ...
                      sum(obj.s_vs_n_sig.bip_elec_data{1})>0;
        if has_bipolar
            num_sig_bipolar = {sum(obj.s_vs_n_sig.bip_elec_data{1})};
            num_total_bipolar = {length(obj.bip_ch_label)};
            percent_sig_bipolar = {round((num_sig_bipolar{1}/num_total_bipolar{1})*100,2)};
            chans_sig_bipolar = {find(obj.s_vs_n_sig.bip_elec_data{1})};
            names_sig_bipolar = {obj.bip_ch_label(obj.s_vs_n_sig.bip_elec_data{1})};
            ratio_tmp = obj.s_vs_n_p_ratio.bip_elec_data{1};
            p_values_sig_bipolar = {ratio_tmp(obj.s_vs_n_sig.bip_elec_data{1})};
        else
            num_sig_bipolar = {};
            num_total_bipolar = {};
            percent_sig_bipolar = {};
            chans_sig_bipolar = {};
            names_sig_bipolar = {};
            p_values_sig_bipolar = {};
        end

        native_sample_freq = {obj.for_preproc.sample_freq_raw};
        elecs_per_amp = {obj.for_preproc.elecs_per_amp};
        if ~isempty(ops.sessions)
            runs_analyzed = {ops.sessions};
        else
            runs_analyzed = {length(obj.stitch_index)};
        end
        num_runs_total = {length(obj.stitch_index)};
        num_words = sum(cell2mat(cellfun(@(x) contains(x,'word'),obj.trial_timing{keep_trials(1),1}.key,'UniformOutput',false)));
        tt0 = obj.trial_timing{keep_trials(1),1};
        w1Row = find(strcmp(tt0.key, 'word_1'), 1);
        if isempty(w1Row)
            w1Row = 1;
        end
        presentation_rate = {(tt0.end(w1Row) - tt0.start(w1Row) + 1) / obj.sample_freq};

        % Flexible condition counting (works with any condition name mapping)
        S_flag = obj.s_vs_n_ops.S_condition_flag;
        N_flag = obj.s_vs_n_ops.N_condition_flag;
        num_trials_S = {sum(strcmp(obj.condition(keep_trials), S_flag))};
        num_trials_N = {sum(strcmp(obj.condition(keep_trials), N_flag))};

        % W and J may not exist for all tasks
        if isfield(obj.s_vs_n_ops, 'W_condition_flag') && ~isempty(obj.s_vs_n_ops.W_condition_flag)
            num_trials_W = {sum(strcmp(obj.condition(keep_trials), obj.s_vs_n_ops.W_condition_flag))};
        else
            num_trials_W = {0};
        end
        if isfield(obj.s_vs_n_ops, 'J_condition_flag') && ~isempty(obj.s_vs_n_ops.J_condition_flag)
            num_trials_J = {sum(strcmp(obj.condition(keep_trials), obj.s_vs_n_ops.J_condition_flag))};
        else
            num_trials_J = {0};
        end

        new_row = {obj.subject,...
                   num_sig, num_sig_left, num_sig_right,...
                   num_clean, num_clean_left, num_clean_right,...
                   num_prelim_deselect, num_w_IED, num_w_noise, num_user_deselect,...
                   num_total_ecog_seeg, num_total_ecog_seeg_left, num_total_ecog_seeg_right,...
                   num_total, percent_sig, chans_sig, names_sig, p_values_sig,...
                   has_bipolar,...
                   num_sig_bipolar, num_total_bipolar, percent_sig_bipolar,...
                   chans_sig_bipolar, names_sig_bipolar, p_values_sig_bipolar,...
                   native_sample_freq, elecs_per_amp, runs_analyzed, num_runs_total,...
                   num_words, presentation_rate,...
                   num_trials_S, num_trials_N, num_trials_W, num_trials_J};
        summary = cell2table(new_row,'VariableNames',columns);
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT TIMECOURSE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_timecourse(obj,varargin)
        p = inputParser();
        addRequired(p,'data');
        addParameter(p,'subAverage',false);
        addParameter(p,'elec_flag','unipolar');
        addParameter(p,'doBipolar',true);
        addParameter(p,'words',1:12);
        addParameter(p,'noDisplay',true);
        parse(p, varargin{:});
        ops = p.Results;
    
        if ~isempty(ops.data.W_values)
            conditions = {'S','W','J','N'};
            labels = {'Sentences','Words','Jabberwocky','Nonwords'};
            colors = {'r',[0.4660 0.6740 0.1880],[1 0.725 0],'b'};
        else
            conditions = {'S','N'};
            labels = {'Sentences','Nonwords'};
            colors = {'r','b'};
        end
    
        if strcmp(ops.elec_flag,'unipolar')
            chan_names = obj.elec_ch_label(logical(obj.s_vs_n_sig.elec_data{1,1}));
            chan_nums = obj.elec_ch(logical(obj.s_vs_n_sig.elec_data{1,1}));
            idx = 1;
        elseif strcmp(ops.elec_flag,'bipolar')
            chan_names = obj.bip_ch_label(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            chan_nums = obj.bip_ch(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            idx = 2;
        else
            error('Signal flag not recognized');
        end
        chan_names = cellfun(@(x) strrep(x,'_',''),chan_names,'UniformOutput',false);
        
        if ops.subAverage
            nPlots = 1;
        else
            nPlots = size(ops.data.S_values{idx},1);
        end
    
        for i=1:nPlots
            close all
            if ops.noDisplay
                set(0, 'DefaultFigureVisible', 'off')
            end
            f = figure; 
            hold on;
            set(f,'position',[1123 29 3000 1275])
    
            if ops.subAverage
                timecourse_data = ops.data.S_values{idx};
                nSamples = size(timecourse_data, 2);
                ave_data = ops.data.S_values_ave{idx};
                nWords = size(ave_data, 2);
            else
                timecourse_data = ops.data.S_values{idx}(i,:);
                nSamples = length(timecourse_data);
                ave_data = ops.data.S_values_ave{idx}(i,:);
                nWords = length(ave_data);
            end
    
            word_boundaries = (1:nSamples/nWords:nSamples)/obj.sample_freq;
            for j=1:length(word_boundaries)
                xline(word_boundaries(j),'--');
            end
            xlim([0 nSamples/obj.sample_freq]);
    
            x_ave = ((1:nWords) - 0.5) * (nSamples/nWords) / obj.sample_freq;
            xticks(x_ave);
            xticklabels(ops.words(1:nWords));
            
            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('ops.data.',conditions{j},'_values_ave{idx};'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1));
                    y = mean(y,1);
                else
                    y = eval(strcat('ops.data.',conditions{j},'_values_ave{idx}(i,:);'));
                    y_sem = eval(strcat('ops.data.',conditions{j},'_values_ave_sem{idx}(i,:);'));
                end
                x_plot = x_ave(1:length(y));
                b1(j) = obj.plot_condition(x_plot,y,y_sem,colors{j},labels{j},'isAveraged',true);
            end
    
            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('ops.data.',conditions{j},'_values{idx};'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1));
                    y = mean(y,1);
                else
                    y = eval(strcat('ops.data.',conditions{j},'_values{idx}(i,:);'));
                    y_sem = eval(strcat('ops.data.',conditions{j},'_values_sem{idx}(i,:);'));
                end
                x_timecourse = (1:length(y)) / obj.sample_freq;
                b2 = obj.plot_condition(x_timecourse,y,y_sem,colors{j},labels{j},'isAveraged',false);
            end
    
            set(gca,'FontSize',18,'LineWidth',1.5,'Box','off');
            legend(b1,'Location','northwest','NumColumns',2,'FontSize',24)
            xlabel('Word Position');
            ylabel('High Gamma Envelope (a.u.)');
            
            if ops.subAverage
                subtxt = sprintf('Average of language-responsive channels (n = %d, %s)',size(ops.data.S_values{idx},1),ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep);
                fname = sprintf('%s_average_s_v_n_%s_timecourse.png',obj.subject,ops.elec_flag);
            else
                subtxt = sprintf('%s (chan%d, nTrials = %d, %s)',chan_names{i},chan_nums(i),ops.data.nTrialsS,ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep,obj.subject,filesep,'timecourse');
                fname = sprintf('%s_%s_%s.png',obj.subject,ops.elec_flag,chan_names{i});
            end
            title(obj.subject);
            subtitle(subtxt)
    
            if ~exist(analysis_path, 'dir')
                mkdir(analysis_path);
            end   
            saveas(gcf,strcat(analysis_path,'/',fname));
            set(0, 'CurrentFigure', f);
            close(f);
        end
    
        if ~isempty(obj.bip_elec_data) && ops.doBipolar
            obj.plot_timecourse(ops.data,'subAverage',ops.subAverage,'words',ops.words,'elec_flag','bipolar','doBipolar',false);
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT BARPLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_barplot(obj,varargin)
        p = inputParser();
        addRequired(p,'data');
        addParameter(p,'subAverage',false);
        addParameter(p,'elec_flag','unipolar');
        addParameter(p,'doBipolar',true);
        addParameter(p,'words',1:12);
        addParameter(p,'noDisplay',true);
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.data.W_values)
            conditions = {'S','W','J','N'};
            labels = {'Sentences','Words','Jabberwocky','Nonwords'};
            colors = {'r',[0.4660 0.6740 0.1880],[1 0.725 0],'b'};
        else
            conditions = {'S','N'};
            labels = {'Sentences','Nonwords'};
            colors = {'r','b'};
        end

        if strcmp(ops.elec_flag,'unipolar')
            chan_names = obj.elec_ch_label(logical(obj.s_vs_n_sig.elec_data{1,1}));
            chan_nums = obj.elec_ch(logical(obj.s_vs_n_sig.elec_data{1,1}));
            idx = 1;
        elseif strcmp(ops.elec_flag,'bipolar')
            chan_names = obj.bip_ch_label(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            chan_nums = obj.bip_ch(logical(obj.s_vs_n_sig.bip_elec_data{1,1}));
            idx = 2;
        else
            error('Signal flag not recognized');
        end
        chan_names = cellfun(@(x) strrep(x,'_',''),chan_names,'UniformOutput',false);
        
        if ops.subAverage
            nPlots = 1;
        else
            nPlots = size(ops.data.S_values{idx},1);
        end

        for i=1:nPlots
            close all
            if ops.noDisplay
                set(0, 'DefaultFigureVisible', 'off')
            end
            f = figure; 
            hold on;
            set(f,'position',[1123 29 500 600])

            for j=1:length(conditions)
                if ops.subAverage
                    y = eval(strcat('mean(ops.data.',conditions{j},'_values_ave{idx},2);'));
                    y_sem = std(y,[],1) ./ sqrt(size(y,1));
                    y = mean(y,1);
                else
                    y = eval(strcat('mean(ops.data.',conditions{j},'_values_ave{idx}(i,:),2);'));
                    y_sem = eval(strcat('mean(ops.data.',conditions{j},'_values_ave_sem{idx}(i,:),2);'));
                end
                h(j) = bar(j,y,'displayname',labels{j});
                set(h(j), 'FaceColor', colors{j})
                errorbar(j,y,y_sem,'k','LineWidth',1.5,'Capsize',0);
            end
        
            ax1 = gca;                
            ax1.XAxis.Visible = 'off'; 
            ylim([-0.5,1])    
            set(gca,'FontSize',10,'LineWidth',1.5,'Box','off');
            legend(h,'FontSize',10,'Location','southwest');
            legend boxoff
            ylabel('Z-Scored High Gamma Envelope (a.u.)','FontSize',10)
            if ops.subAverage
                subtxt = sprintf('Average of language-responsive channels (n = %d, %s)',size(ops.data.S_values{idx},1),ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep);
                fname = sprintf('%s_average_s_v_n_%s_barplot.png',obj.subject,ops.elec_flag);
            else
                subtxt = sprintf('%s (chan%d, nTrials = %d, %s)',chan_names{i},chan_nums(i),ops.data.nTrialsS,ops.elec_flag);
                analysis_path = strcat(obj.langloc_save_path,'lang_electrodes',filesep,obj.subject,filesep,'barplot');
                fname = sprintf('%s_%s_%s.png',obj.subject,ops.elec_flag,chan_names{i});
            end
            title(obj.subject);
            subtitle(subtxt,'FontSize',8)
        
            if ~exist(analysis_path, 'dir')
                mkdir(analysis_path);
            end   
            saveas(gcf,strcat(analysis_path,'/',fname));
            set(0, 'CurrentFigure', f);
            close(f);
        end

        if ~isempty(obj.bip_elec_data) && ops.doBipolar
            obj.plot_barplot(ops.data,'subAverage',ops.subAverage,'words',ops.words,'elec_flag','bipolar','doBipolar',false);
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT CONDITION (shared helper)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function b=plot_condition(obj,varargin)
        p = inputParser();
        addRequired(p,'x');
        addRequired(p,'y');
        addRequired(p,'y_sem');
        addRequired(p,'c');
        addRequired(p,'label');
        addParameter(p,'isAveraged',false); 
        addParameter(p,'linewidthAveraged',9);
        addParameter(p,'markersizeAveraged',20);
        addParameter(p,'linewidthErrorAveraged',3);
        parse(p, varargin{:});
        ops = p.Results;

        if ops.isAveraged
            color           = ops.c;
            linewidth       = ops.linewidthAveraged;
            marker          = 's';
            markersize      = ops.markersizeAveraged;
            label           = ops.label;
            linestyle_err   = 'none';
            linewidth_err   = ops.linewidthErrorAveraged;
            capsize_err     = 0;
        else
            color           = ops.c;
            linewidth       = 2;
            marker          = 'default';
            markersize      = 'default';
            label           = '';
            linestyle_err   = 'none';
            facealpha_err   = 0.1;
        end

        b = plot(ops.x,ops.y,...
                        'Color',color,...
                        'LineWidth',linewidth,...
                        'Marker',marker,...
                        'MarkerSize',markersize,...
                        'MarkerFaceColor',color,...
                        'MarkerEdgeColor',color,...
                        'DisplayName',label...
        );

        if ops.isAveraged
            errorbar(ops.x,ops.y,ops.y_sem,...
                     'Color',color,...
                     'LineStyle',linestyle_err,...
                     'LineWidth',linewidth_err,...
                     'CapSize',capsize_err...
            );
        else
            patch([ops.x fliplr(ops.x)], [ops.y+ops.y_sem fliplr(ops.y-ops.y_sem)],...
                  color,...
                  'LineStyle',linestyle_err,...
                  'FaceAlpha',facealpha_err...
            );
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT S VS N (all electrode overview)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_s_vs_n(obj,prev_ops,varargin)
        p = inputParser();
        addRequired(p,'elec_flag');
        addRequired(p,'channel_labels')
        addParameter(p,'S_table',[])
        addParameter(p,'N_table',[]);
        addParameter(p,'S_N_corr_table',[]);
        addParameter(p,'S_N_corr_rnd_table',[]);
        addParameter(p,'s_vs_n_sig',[]);
        addParameter(p,'S_N_p_ratio_tbl',[]);
        addParameter(p,'noDisplay',true)
        parse(p, varargin{:});
        ops = p.Results;
    
        analysis_path = strcat(obj.langloc_save_path,'all_electrodes/');
    
        elec_flag = ops.elec_flag;
        S_dat = ops.S_table.(elec_flag){1};
        N_dat = ops.N_table.(elec_flag){1};
        SN_corr_dat = ops.S_N_corr_table.(elec_flag){1};
        SN_cor_rnd_dat = ops.S_N_corr_rnd_table.(elec_flag){1};
        SN_sig_dat = ops.s_vs_n_sig.(elec_flag){1};
        SN_p_ratio_dat = ops.S_N_p_ratio_tbl.(elec_flag){1};
    
        nElecs = size(S_dat, 1);
        nTrials_S = size(S_dat, 2);
        nTrials_N = size(N_dat, 2);
        nWords = size(S_dat, 3);
    
        num_rows = 3;
        num_columns = 3;
        nbins = 50;
        pp = 0;
    
        fprintf(1,'\n> Plotting channels from %s ...\n',elec_flag);
        fprintf(1,'[');
    
        if strcmp(ops.elec_flag,'elec_data')
            elec_flag = 'unipolar';
        elseif strcmp(ops.elec_flag,'bip_elec_data')
            elec_flag = 'bipolar';
        else
            error('Signal flag not recognized');
        end
    
        close all
        if ops.noDisplay
            set(0, 'DefaultFigureVisible', 'off')
        end
        f = figure;
        set(f,'position',[1123 29 1266 1275])
            
        for i=1:nElecs
            if nWords == 1
                s_raw = S_dat(i, :);
                n_raw = N_dat(i, :);
                s_electrode_resp = reshape(s_raw, 1, []);
                n_electrode_resp = reshape(n_raw, 1, []);
            else
                s_raw = squeeze(S_dat(i, :, :));
                n_raw = squeeze(N_dat(i, :, :));
                if size(s_raw, 1) == nTrials_S && size(s_raw, 2) == nWords
                    s_electrode_resp = s_raw';
                else
                    s_electrode_resp = s_raw;
                end
                if size(n_raw, 1) == nTrials_N && size(n_raw, 2) == nWords
                    n_electrode_resp = n_raw';
                else
                    n_electrode_resp = n_raw;
                end
            end
            
            s_n_rho = SN_corr_dat(i);
            s_n_rho_rnd = SN_cor_rnd_dat(i,:);
            is_sig = SN_sig_dat(i);
            p_ratio = SN_p_ratio_dat(i);
            
            nWords_actual = size(s_electrode_resp, 1);
            nTrials_S_actual = size(s_electrode_resp, 2);
            nTrials_N_actual = size(n_electrode_resp, 2);
            word_pos_s = repmat((1:nWords_actual)', 1, nTrials_S_actual);
            word_pos_n = repmat((1:nWords_actual)', 1, nTrials_N_actual);
              
            if isnan(s_n_rho)
                continue
            end
    
            sup_title = (strcat(ops.channel_labels{i,1}));
    
            ax = subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+1);
            h1 = histogram(s_n_rho_rnd,nbins);
            hold on
            xline(s_n_rho,'linewidth',3);
            ax.Box = 'off';
            h1.EdgeColor = 'w';
            if ~ops.noDisplay; shg; end
            ax.XAxis.LineWidth = 2;
            ax.YAxis.LineWidth = 2;
            ax.Title.String=sprintf('corr=%f,\n p_{ratio}=%0.4f sig=%d',s_n_rho,p_ratio,is_sig);
    
            ax = subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+2);
            
            x_pos_s = mean(word_pos_s, 2);
            x_pos_n = mean(word_pos_n, 2);
            y_mean_s = mean(s_electrode_resp, 2);
            y_mean_n = mean(n_electrode_resp, 2);
            x_pos_s = x_pos_s(:); x_pos_n = x_pos_n(:);
            y_mean_s = y_mean_s(:); y_mean_n = y_mean_n(:);
            
            b1 = plot(x_pos_s + 0.1, y_mean_s, 'color', [1,.5,.5], 'linewidth', 2, 'marker', 'o', ...
                      'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0], 'displayname', 'S');
            hold on
            b2 = plot(x_pos_n - 0.1, y_mean_n, 'color', [.5,.5,1], 'linewidth', 2, 'marker', 'o', ...
                      'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', [0,0,1], 'displayname', 'N');
            
            y_sem_s = std(s_electrode_resp, [], 2) ./ sqrt(nTrials_S_actual);
            y_sem_s = y_sem_s(:);
            bl = errorbar(x_pos_s + 0.1, y_mean_s, y_sem_s);
            bl.LineStyle = 'none'; bl.Color = [1,.5,.5]; bl.LineWidth = 2; bl.CapSize = 2;
            hAnnotation = get(bl,'Annotation'); hLegendEntry = get(hAnnotation,'LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off');
            
            y_sem_n = std(n_electrode_resp, [], 2) ./ sqrt(nTrials_N_actual);
            y_sem_n = y_sem_n(:);
            bl = errorbar(x_pos_n - 0.1, y_mean_n, y_sem_n);
            bl.LineStyle = 'none'; bl.Color = [.5,.5,1]; bl.LineWidth = 2; bl.CapSize = 2;
            hAnnotation = get(bl,'Annotation'); hLegendEntry = get(hAnnotation,'LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off');
                    
            ax.XAxis.Visible = 'on';
            ax.XTick = 1:nWords_actual;
            ax.XLim = [0, nWords_actual + 1];
            all_points = [s_electrode_resp(:); n_electrode_resp(:)];
            y_quantile = quantile(all_points, 10);
            ax.FontSize = 12;
            set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
            ah = get(ax,'children');
            arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
            ah(3).DisplayName = 'N'; ah(4).DisplayName = 'S';
            set(ax,'children',ah);
            legend(ah(3:4),'Location','northwest','NumColumns',2)
            xlabel('word position');
            ax.YLabel.String='High Gamma (a.u.)';
            ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
            ax.Title.String = erase(sup_title,'_');
    
            if ~mod(i,num_rows) || i==nElecs
                pp = pp+1;
                if ~exist(strcat(analysis_path,obj.subject), 'dir')
                    mkdir(strcat(analysis_path,obj.subject));
                end
                set(gcf,'PaperPosition',[.25 .25 8 6])
                set(gcf,'PaperOrientation','landscape');
                fname = sprintf('%s_s_v_n_words-%d-%d_p_%0.2f_%s_%s_%s.pdf',...
                                obj.subject,...
                                min(prev_ops.words),...
                                max(prev_ops.words),...
                                prev_ops.threshold,...
                                prev_ops.side,...
                                elec_flag,...
                                num2str(pp)...
                );
                print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.subject,'/',fname));
                set(0, 'CurrentFigure', f);
                clf reset;
            end    
            fprintf(1,'.');
        end
                
        close(f);
        fprintf(1,'] done\n')
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WORD-WISE S VS N  (per-word permutation + consecutiveness criterion)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function test_s_vs_n_wordwise(obj, varargin)
        % Per-word-position permutation test. A channel is called
        % language-responsive if it shows significance in a run of
        % min_consecutive words (consecutiveness=true) or in at least
        % min_consecutive words anywhere (consecutiveness=false).
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
        addParameter(p,'use_odd_for_inference',true);
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session==ops.sessions,2);
        else
            keep_trials = [];
        end

        [~, S_ave_cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition', ...
            ops.S_condition_flag,'keep_trials',keep_trials);
        [~, N_ave_cond_table] = obj.get_ave_cond_trial('words',ops.words,'condition', ...
            ops.N_condition_flag,'keep_trials',keep_trials);

        for elec_type = 1:2
            if elec_type == 1
                S_data = S_ave_cond_table.elec_data{1};
                N_data = N_ave_cond_table.elec_data{1};
            else
                if isempty(S_ave_cond_table.bip_elec_data{1})
                    fprintf('\nNo bipolar data found, skipping bipolar.\n');
                    continue;
                end
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

            if(elec_type==1); strval = 'Unipolar'; else; strval = 'Bipolar'; end
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
                    if any(conv(double(sig_vec), ones(1, ops.min_consecutive), 'valid') == ops.min_consecutive)
                        langloc(elec) = true;
                    end
                else
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


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WORD-BOUNDARIES S VS N  (cluster-based time-series permutation test)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function test_s_vs_n_wordboundaries(obj, varargin)
        % Cluster-based time-series permutation test (timePermCluster) run on
        % an epoch around each word position. Requires timePermCluster.m on
        % the path (provided in kumar_ieeg_utils/).
        %
        % NOTE: trial_timing row 1 is fixation (key = fix); word position N
        % maps to key word_N (MGH / ieeg_pipeline convention).
        p = inputParser();
        addParameter(p, 'S_condition_flag', 'sentence');
        addParameter(p, 'N_condition_flag', 'nonword');
        addParameter(p, 'n_rep', 1000);
        addParameter(p, 'threshold', 0.05);
        addParameter(p, 'sessions', []);
        addParameter(p, 'do_plot', false);
        addParameter(p, 'epoch_range', [-0.5, 0.5]);
        addParameter(p, 'num_words', 12);
        addParameter(p, 'use_odd_for_inference', true);
        parse(p, varargin{:});
        ops = p.Results;

        obj.s_vs_n_wordboundaries_ops = ops;

        fprintf(1, '\n> Running word boundaries time-series permutation test for %s vs %s ...\n', ...
                ops.S_condition_flag, ops.N_condition_flag);

        availableConditions = unique(obj.condition);
        if ~ismember(ops.S_condition_flag, availableConditions) || ~ismember(ops.N_condition_flag, availableConditions)
            error('Required conditions "%s" and "%s" not found. Available conditions: %s', ...
                  ops.S_condition_flag, ops.N_condition_flag, strjoin(availableConditions, ', '));
        end

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session == ops.sessions, 2);
        else
            keep_trials = true(size(obj.condition));
        end

        pSig = cell(size(obj.elec_ch_label, 1), ops.num_words);
        pSig_bip = cell(size(obj.bip_ch_label, 1), ops.num_words);

        concatenatedEpochsSentence = [];
        concatenatedEpochsNonword = [];
        concatenatedEpochsSentence_bip = [];
        concatenatedEpochsNonword_bip = [];

        fprintf(1, '> Processing word boundaries for %s experiment...\n', obj.experiment);

        for wordPos = 1:ops.num_words
            fprintf(1, '  Processing word position %d/%d...\n', wordPos, ops.num_words);

            [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', ops.epoch_range, ...
                                                                  'key', sprintf('word_%d', wordPos));

            sentence_trials_idx = strcmp(obj.condition, ops.S_condition_flag) & keep_trials;
            nonword_trials_idx = strcmp(obj.condition, ops.N_condition_flag) & keep_trials;

            % --- ODD/EVEN SPLIT ---
            if ops.use_odd_for_inference
                s_idx = find(sentence_trials_idx);
                n_idx = find(nonword_trials_idx);
                sentence_trials_idx = false(size(sentence_trials_idx));
                nonword_trials_idx  = false(size(nonword_trials_idx));
                sentence_trials_idx(s_idx(1:2:end)) = true;
                nonword_trials_idx(n_idx(1:2:end))  = true;
            end

            sentenceTrials = epochData(:, sentence_trials_idx, :);
            nonwordTrials = epochData(:, nonword_trials_idx, :);

            if isempty(concatenatedEpochsSentence)
                concatenatedEpochsSentence = sentenceTrials;
                concatenatedEpochsNonword = nonwordTrials;
            else
                concatenatedEpochsSentence = cat(3, concatenatedEpochsSentence, sentenceTrials);
                concatenatedEpochsNonword = cat(3, concatenatedEpochsNonword, nonwordTrials);
            end

            fprintf(1, '    Running permutation tests for unipolar channels...\n');
            parfor iChan = 1:size(epochData, 1)
                aTrialData = squeeze(epochData(iChan, sentence_trials_idx, :));
                bTrialData = squeeze(epochData(iChan, nonword_trials_idx, :));

                if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                    try
                        pSig{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'nPerm', ops.n_rep, 'statstype','corr','pThresh',ops.threshold);
                    catch ME
                        warning('TimePermCluster failed for unipolar channel %d, word %d: %s', ...
                               iChan, wordPos, ME.message);
                        pSig{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData, 3)), ...
                                                     'p_val', ones(1, size(epochData, 3)));
                    end
                else
                    pSig{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData, 3)), ...
                                                 'p_val', ones(1, size(epochData, 3)));
                end
            end

            if ~isempty(epochData_bip)
                fprintf(1, '    Running permutation tests for bipolar channels...\n');

                sentenceTrials_bip = epochData_bip(:, sentence_trials_idx, :);
                nonwordTrials_bip = epochData_bip(:, nonword_trials_idx, :);

                if isempty(concatenatedEpochsSentence_bip)
                    concatenatedEpochsSentence_bip = sentenceTrials_bip;
                    concatenatedEpochsNonword_bip = nonwordTrials_bip;
                else
                    concatenatedEpochsSentence_bip = cat(3, concatenatedEpochsSentence_bip, sentenceTrials_bip);
                    concatenatedEpochsNonword_bip = cat(3, concatenatedEpochsNonword_bip, nonwordTrials_bip);
                end

                parfor iChan = 1:size(epochData_bip, 1)
                    aTrialData = squeeze(epochData_bip(iChan, sentence_trials_idx, :));
                    bTrialData = squeeze(epochData_bip(iChan, nonword_trials_idx, :));

                    if size(aTrialData, 1) > 1 && size(bTrialData, 1) > 1
                        try
                            pSig_bip{iChan, wordPos} = timePermCluster(aTrialData, bTrialData, 'nPerm', ops.n_rep, 'statstype','corr','pThresh',ops.threshold);
                        catch ME
                            warning('TimePermCluster failed for bipolar channel %d, word %d: %s', ...
                                   iChan, wordPos, ME.message);
                            pSig_bip{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData_bip, 3)), ...
                                                             'p_val', ones(1, size(epochData_bip, 3)));
                        end
                    else
                        pSig_bip{iChan, wordPos} = struct('h_sig_05', zeros(1, size(epochData_bip, 3)), ...
                                                         'p_val', ones(1, size(epochData_bip, 3)));
                    end
                end
            end
        end

        timePointsPerWord = size(concatenatedEpochsSentence, 3) / ops.num_words;
        totalTimePoints = timePointsPerWord * ops.num_words;
        wordBoundaries = 0:timePointsPerWord:totalTimePoints;

        obj.s_vs_n_wordboundaries_pSig = pSig;
        obj.s_vs_n_wordboundaries_wordBoundaries = wordBoundaries;
        obj.s_vs_n_wordboundaries_timePointsPerWord = timePointsPerWord;
        obj.s_vs_n_wordboundaries_totalTimePoints = totalTimePoints;

        if ~isempty(concatenatedEpochsSentence_bip)
            obj.s_vs_n_wordboundaries_pSig_bip = pSig_bip;
        end

        sigChannels = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), pSig'));
        sigUnipolarChannels = find(sigChannels > 0);

        fprintf(1, '> Found %d/%d unipolar channels with significant word boundary clusters\n', ...
                length(sigUnipolarChannels), size(concatenatedEpochsSentence, 1));

        if ~isempty(concatenatedEpochsSentence_bip)
            sigChannels_bip = sum(cellfun(@(x) any(cell2mat(arrayfun(@(y) any(y.h_sig_05), x, 'UniformOutput', false))), pSig_bip'));
            sigBipolarChannels = find(sigChannels_bip > 0);

            fprintf(1, '> Found %d/%d bipolar channels with significant word boundary clusters\n', ...
                    length(sigBipolarChannels), size(concatenatedEpochsSentence_bip, 1));

            obj.s_vs_n_wordboundaries_sigBipolarChannels = sigBipolarChannels;
        end

        obj.s_vs_n_wordboundaries_sigUnipolarChannels = sigUnipolarChannels;

        obj.s_vs_n_wordboundaries_concatenatedEpochsSentence = concatenatedEpochsSentence;
        obj.s_vs_n_wordboundaries_concatenatedEpochsNonword = concatenatedEpochsNonword;

        if ~isempty(concatenatedEpochsSentence_bip)
            obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip = concatenatedEpochsSentence_bip;
            obj.s_vs_n_wordboundaries_concatenatedEpochsNonword_bip = concatenatedEpochsNonword_bip;
        end

        if ops.do_plot
            obj.plot_s_vs_n_wordboundaries(ops);
        end

        fprintf(1, '> Word boundaries time-series permutation test completed\n');
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT WORD-BOUNDARIES RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_s_vs_n_wordboundaries(obj, ops)
        concatenatedEpochsSentence = obj.s_vs_n_wordboundaries_concatenatedEpochsSentence;
        concatenatedEpochsNonword = obj.s_vs_n_wordboundaries_concatenatedEpochsNonword;
        pSig = obj.s_vs_n_wordboundaries_pSig;
        wordBoundaries = obj.s_vs_n_wordboundaries_wordBoundaries;
        timePointsPerWord = obj.s_vs_n_wordboundaries_timePointsPerWord;
        totalTimePoints = obj.s_vs_n_wordboundaries_totalTimePoints;

        numChanBlock = 5;
        totChanBlock = ceil(size(concatenatedEpochsSentence, 1) / numChanBlock);

        colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];

        x = 1:totalTimePoints;
        numWords = size(pSig, 2);

        for iChanBlock = 0:(totChanBlock-1)
            f = figure('Visible', 'off', 'Position', [100 100 1200 1400], 'Renderer', 'painters');

            for iChan = 1:min(numChanBlock, size(concatenatedEpochsSentence, 1) - iChanBlock*numChanBlock)
                iChan2 = iChanBlock*numChanBlock + iChan;
                subplot(numChanBlock, 1, iChan);
                hold on;

                title(sprintf('%s', obj.elec_ch_label{iChan2}), 'Interpreter', 'none', ...
                      'FontSize', 10, 'FontWeight', 'bold');

                trialData_sentence = squeeze(concatenatedEpochsSentence(iChan2, :, :));
                trialMean_sentence = nanmean(trialData_sentence, 1);
                trialSEM_sentence = nanstd(trialData_sentence, 0, 1) / sqrt(size(trialData_sentence, 1));
                maxval_sentence = max(trialMean_sentence + trialSEM_sentence);

                trialData_nonword = squeeze(concatenatedEpochsNonword(iChan2, :, :));
                trialMean_nonword = nanmean(trialData_nonword, 1);
                trialSEM_nonword = nanstd(trialData_nonword, 0, 1) / sqrt(size(trialData_nonword, 1));
                maxval_nonword = max(trialMean_nonword + trialSEM_nonword);

                plot(x, trialMean_sentence, 'Color', colors(1,:), 'LineWidth', 2);
                patch([x, fliplr(x)], [trialMean_sentence+trialSEM_sentence, ...
                      fliplr(trialMean_sentence-trialSEM_sentence)], ...
                      colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                plot(x, trialMean_nonword, 'Color', colors(2,:), 'LineWidth', 2);
                patch([x, fliplr(x)], [trialMean_nonword+trialSEM_nonword, ...
                      fliplr(trialMean_nonword-trialSEM_nonword)], ...
                      colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                for boundary = wordBoundaries
                    if boundary > 0 && boundary <= totalTimePoints
                        xline(boundary, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
                    end
                end

                if numWords <= 12
                    midPoints = (wordBoundaries(1:end-1) + wordBoundaries(2:end)) / 2;
                    wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:numWords, 'UniformOutput', false);
                    set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                else
                    midPoints = (wordBoundaries(1:2:end-1) + wordBoundaries(2:2:end)) / 2;
                    wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:2:numWords, 'UniformOutput', false);
                    set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                end

                sigTimePoints = [];
                for wordPos = 1:numWords
                    if wordPos <= size(pSig, 2) && ~isempty(pSig{iChan2, wordPos})
                        wordSigPoints = find(pSig{iChan2, wordPos}.h_sig_05) + (wordPos-1)*timePointsPerWord;
                        sigTimePoints = [sigTimePoints, wordSigPoints];
                    end
                end

                maxval = max(maxval_sentence, maxval_nonword);
                ylim([-1.5, maxval + 0.5]);

                if ~isempty(sigTimePoints)
                    scatter(sigTimePoints, -1.25*ones(size(sigTimePoints)), 10, 'r', 'filled', 'LineWidth', 0.5);
                end

                yline(0, 'k-.', 'LineWidth', 1, 'Alpha', 0.5);

                xlabel('Word Position', 'FontSize', 8);
                ylabel('High Gamma Power', 'FontSize', 8);
                set(gca, 'FontSize', 10);
                grid on; grid minor;

                text(0.02, 0.98, sprintf('Sentence trials: %d\nNonword trials: %d', ...
                    size(trialData_sentence, 1), size(trialData_nonword, 1)), ...
                    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
                    'BackgroundColor', 'white', 'EdgeColor', 'black');

                hold off;
            end

            legendAxes = axes('Position', [0.1, 0.01, 0.8, 0.08], 'Visible', 'off');
            hold(legendAxes, 'on');
            scatter(legendAxes, NaN, NaN, 10, 'r', 'filled');
            plot(legendAxes, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2);
            plot(legendAxes, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2);
            plot(legendAxes, NaN, NaN, 'k--', 'LineWidth', 1.5);
            plot(legendAxes, NaN, NaN, 'k-.', 'LineWidth', 1);
            legendLabels = {'Significant Sentence vs Nonword', ops.S_condition_flag, ops.N_condition_flag, ...
                           'Word Boundaries', 'Baseline'};
            legend(legendAxes, legendLabels, 'Orientation', 'horizontal', ...
                   'Location', 'southoutside', 'FontSize', 8);
            hold(legendAxes, 'off');

            sgtitle(sprintf('%s - Word Boundaries Analysis (UNIPOLAR) - Block %d/%d', ...
                   obj.experiment, iChanBlock+1, totChanBlock), ...
                   'FontSize', 8, 'FontWeight', 'bold');

            if ops.do_plot
                set(f, 'Visible', 'on');
            end
        end

        if isprop(obj, 's_vs_n_wordboundaries_concatenatedEpochsSentence_bip') && ...
           ~isempty(obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip)

            concatenatedEpochsSentence_bip = obj.s_vs_n_wordboundaries_concatenatedEpochsSentence_bip;
            concatenatedEpochsNonword_bip = obj.s_vs_n_wordboundaries_concatenatedEpochsNonword_bip;
            pSig_bip = obj.s_vs_n_wordboundaries_pSig_bip;

            totChanBlock_bip = ceil(size(concatenatedEpochsSentence_bip, 1) / numChanBlock);

            for iChanBlock = 0:(totChanBlock_bip-1)
                f = figure('Visible', 'off', 'Position', [100 100 1200 1400], 'Renderer', 'painters');

                for iChan = 1:min(numChanBlock, size(concatenatedEpochsSentence_bip, 1) - iChanBlock*numChanBlock)
                    iChan2 = iChanBlock*numChanBlock + iChan;
                    subplot(numChanBlock, 1, iChan);
                    hold on;

                    title(sprintf('%s', obj.bip_ch_label{iChan2}), 'Interpreter', 'none', ...
                          'FontSize', 10, 'FontWeight', 'bold');

                    trialData_sentence = squeeze(concatenatedEpochsSentence_bip(iChan2, :, :));
                    trialMean_sentence = nanmean(trialData_sentence, 1);
                    trialSEM_sentence = nanstd(trialData_sentence, 0, 1) / sqrt(size(trialData_sentence, 1));
                    maxval_sentence = max(trialMean_sentence + trialSEM_sentence);

                    trialData_nonword = squeeze(concatenatedEpochsNonword_bip(iChan2, :, :));
                    trialMean_nonword = nanmean(trialData_nonword, 1);
                    trialSEM_nonword = nanstd(trialData_nonword, 0, 1) / sqrt(size(trialData_nonword, 1));
                    maxval_nonword = max(trialMean_nonword + trialSEM_nonword);

                    plot(x, trialMean_sentence, 'Color', colors(1,:), 'LineWidth', 2);
                    patch([x, fliplr(x)], [trialMean_sentence+trialSEM_sentence, ...
                          fliplr(trialMean_sentence-trialSEM_sentence)], ...
                          colors(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                    plot(x, trialMean_nonword, 'Color', colors(2,:), 'LineWidth', 2);
                    patch([x, fliplr(x)], [trialMean_nonword+trialSEM_nonword, ...
                          fliplr(trialMean_nonword-trialSEM_nonword)], ...
                          colors(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

                    for boundary = wordBoundaries
                        if boundary > 0 && boundary <= totalTimePoints
                            xline(boundary, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
                        end
                    end

                    if numWords <= 12
                        midPoints = (wordBoundaries(1:end-1) + wordBoundaries(2:end)) / 2;
                        wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:numWords, 'UniformOutput', false);
                        set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                    else
                        midPoints = (wordBoundaries(1:2:end-1) + wordBoundaries(2:2:end)) / 2;
                        wordLabels = arrayfun(@(x) sprintf('Wd%d', x), 1:2:numWords, 'UniformOutput', false);
                        set(gca, 'XTick', midPoints, 'XTickLabel', wordLabels);
                    end

                    sigTimePoints = [];
                    for wordPos = 1:numWords
                        if wordPos <= size(pSig_bip, 2) && ~isempty(pSig_bip{iChan2, wordPos})
                            wordSigPoints = find(pSig_bip{iChan2, wordPos}.h_sig_05) + (wordPos-1)*timePointsPerWord;
                            sigTimePoints = [sigTimePoints, wordSigPoints];
                        end
                    end

                    maxval = max(maxval_sentence, maxval_nonword);
                    ylim([-1.5, maxval + 0.5]);

                    if ~isempty(sigTimePoints)
                        scatter(sigTimePoints, -1.25*ones(size(sigTimePoints)), 10, 'r', 'filled', 'LineWidth', 0.5);
                    end

                    yline(0, 'k-.', 'LineWidth', 1, 'Alpha', 0.5);

                    xlabel('Word Position', 'FontSize', 8);
                    ylabel('High Gamma Power', 'FontSize', 8);
                    set(gca, 'FontSize', 9);
                    grid on; grid minor;

                    text(0.02, 0.98, sprintf('Sentence trials: %d\nNonword trials: %d', ...
                        size(trialData_sentence, 1), size(trialData_nonword, 1)), ...
                        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
                        'BackgroundColor', 'white', 'EdgeColor', 'black');

                    hold off;
                end

                legendAxes = axes('Position', [0.1, 0.01, 0.8, 0.08], 'Visible', 'off');
                hold(legendAxes, 'on');
                scatter(legendAxes, NaN, NaN, 30, 'r', 'filled');
                plot(legendAxes, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2);
                plot(legendAxes, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2);
                plot(legendAxes, NaN, NaN, 'k--', 'LineWidth', 1.5);
                plot(legendAxes, NaN, NaN, 'k-.', 'LineWidth', 1);
                legendLabels = {'Significant Sentence vs Nonword', ops.S_condition_flag, ops.N_condition_flag, ...
                               'Word Boundaries', 'Baseline'};
                legend(legendAxes, legendLabels, 'Orientation', 'horizontal', ...
                       'Location', 'southoutside', 'FontSize', 20);
                hold(legendAxes, 'off');

                sgtitle(sprintf('%s - Word Boundaries Analysis (BIPOLAR) - Block %d/%d', ...
                       obj.experiment, iChanBlock+1, totChanBlock_bip), ...
                       'FontSize', 20, 'FontWeight', 'bold');

                if ops.do_plot
                    set(f, 'Visible', 'on');
                end
            end
        end

        fprintf(1, '> Unipolar analysis completed: %d channels, %d significant channels\n', ...
               size(concatenatedEpochsSentence, 1), length(obj.s_vs_n_wordboundaries_sigUnipolarChannels));

        if ~isempty(obj.s_vs_n_wordboundaries_sigBipolarChannels)
            fprintf(1, '> Bipolar analysis completed: %d significant channels\n', ...
                   length(obj.s_vs_n_wordboundaries_sigBipolarChannels));
        end
    end


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HG POWER DIFFERENCE (S - N) effect size on held-out (even) trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function compute_hg_power_diff_s_vs_n(obj, varargin)
        p = inputParser();
        addParameter(p, 'words', 1:12);
        addParameter(p, 'S_condition_flag', 'S');
        addParameter(p, 'N_condition_flag', 'N');
        addParameter(p, 'sessions', []);
        addParameter(p, 'use_odd_for_inference', true);
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session == ops.sessions, 2);
        else
            keep_trials = [];
        end

        [~, S_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.S_condition_flag, 'keep_trials', keep_trials);
        [~, N_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.N_condition_flag, 'keep_trials', keep_trials);

        for elec_type = 1:2
            if elec_type == 1
                S_data = S_ave_cond_table.elec_data{1};
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

            % --- EVEN trials only (dim 2 = trials); odd reserved for inference ---
            if ops.use_odd_for_inference
                S_data = S_data(:, 2:2:size(S_data,2), :);
                N_data = N_data(:, 2:2:size(N_data,2), :);
            end

            nElecs = size(S_data, 1);
            fprintf('\nComputing HG power difference (S - N) for %s electrodes\n', strval);

            S_mean = mean(mean(S_data, 3), 2);   % elec x 1
            N_mean = mean(mean(N_data, 3), 2);   % elec x 1
            hg_power_diff = S_mean - N_mean;

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


    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HG S-vs-N CORRELATION effect size on held-out (even) trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function compute_hg_sn_corr(obj, varargin)
        p = inputParser();
        addParameter(p, 'words',             1:12);
        addParameter(p, 'S_condition_flag',  'S');
        addParameter(p, 'N_condition_flag',  'N');
        addParameter(p, 'sessions',          []);
        addParameter(p, 'corr_type',         'Spearman');
        parse(p, varargin{:});
        ops = p.Results;

        if ~isempty(ops.sessions)
            keep_trials = any(obj.session == ops.sessions, 2);
        else
            keep_trials = [];
        end

        [~, S_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.S_condition_flag, 'keep_trials', keep_trials);
        [~, N_ave_cond_table] = obj.get_ave_cond_trial('words', ops.words, 'condition', ops.N_condition_flag, 'keep_trials', keep_trials);

        for elec_type = 1:2
            if elec_type == 1
                S_data = S_ave_cond_table.elec_data{1};
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

            % --- EVEN trials only; odd reserved for test_s_vs_n inference ---
            S_data = S_data(:, 2:2:size(S_data, 2), :);
            N_data = N_data(:, 2:2:size(N_data, 2), :);

            nElecs    = size(S_data, 1);
            nS_trials = size(S_data, 2);
            nN_trials = size(N_data, 2);

            fprintf('\nComputing S-vs-N %s correlation for %s electrodes (n=%d) ...\n', ...
                ops.corr_type, strval, nElecs);

            S_mat = mean(S_data, 3);
            N_mat = mean(N_data, 3);

            combined = [S_mat, N_mat];
            flag     = [ones(1, nS_trials), -ones(1, nN_trials)];

            obs_corr = zeros(nElecs, 1);
            for e = 1:nElecs
                obs_corr(e) = corr(combined(e,:)', flag', 'Type', ops.corr_type);
            end

            res = struct();
            res.corr = obs_corr;

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


end

end
