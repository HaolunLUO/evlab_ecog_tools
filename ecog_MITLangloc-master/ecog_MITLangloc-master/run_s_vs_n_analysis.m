
% Script to load crunched HG_zscore data, run S vs N analysis, and save updated ecog_sn object
function run_s_vs_n_analysis()
    % Define base directory path
    base_dir = '/Volumes/disk/nese/LangLoc/data/derivatives/visual/'; % MODIFY THIS PATH
    save_dir = '/Volumes/disk/nese/LangLoc/data/crunched/visual/';
    % List all subject folders in base directory
    dirs = dir(fullfile(base_dir, 'sub-EM*'));
    subject_ids = {dirs([dirs.isdir]).name};

    % Loop over each subject
    for i = 1:length(subject_ids)
        subject_id = subject_ids{i};
        fprintf('/n Processing subject: %s', subject_id);
        try
            subject_dir = fullfile(base_dir, subject_id);
            crunched_dir = fullfile(subject_dir, 'preproc', 'crunched');

            % Find all crunched HG ZScore files
            files = dir(fullfile(crunched_dir, sprintf('%s*crunched_HG_ZScore.mat', subject_id)));
            if isempty(files)
                fprintf('/n No crunched_HG_ZScore files found for %s, skipping...', subject_id);
                continue;
            end

            % Process each file
            for f = 1:length(files)
                crunched_file = fullfile(crunched_dir, files(f).name);

                % Create ecog_sn_data object
                langloc_save_path = fullfile(subject_dir, 'analysis');
                langloc_crunched_file_name = crunched_file;
                langloc_crunched_file_path = crunched_dir;
                preproc_class_file_name = '/Users/dsuseendar/repos/ieeg_pipeline/@ecog_data/ecog_data.m';
                preproc_class_file_path = '/Users/dsuseendar/repos/ieeg_pipeline/@ecog_data';

                fprintf('/n  Loading file: %s', files(f).name);
                sn_obj = ecog_sn_data(langloc_save_path, ...
                                      langloc_crunched_file_name, ...
                                      langloc_crunched_file_path, ...
                                      preproc_class_file_name, ...
                                      preproc_class_file_path);

                % cut into trials if it hasn't already been done
                if isempty(sn_obj.trial_data)
                    sn_obj.make_trials();
                end

                % Run S vs N analysis
                sn_obj.test_s_vs_n('word', 1:12, ...
                                  'S_condition_flag', 'sentence', ...
                                  'N_condition_flag', 'nonword', ...
                                  'n_rep', 1000, ...
                                  'corr_type', 'Spearman', ...
                                  'threshold', 0.05, ...
                                  'side', 'right', ...
                                  'sessions', [], ...
                                  'do_plot', false);

                sn_obj.compute_hg_power_diff_s_vs_n('words', 1:12, ...
                        'S_condition_flag', 'sentence', ...
                        'N_condition_flag', 'nonword', ...
                        'sessions', []);

                sn_obj.compute_hg_sn_corr('words', 1:12, ...
                        'S_condition_flag', 'sentence', ...
                        'N_condition_flag', 'nonword', ...
                        'sessions', []);

                % Run S vs N test
                sn_obj.test_s_vs_n_wordwise('word', 1:12, ...
                    'S_condition_flag', 'sentence', ...
                    'N_condition_flag', 'nonword', ...
                    'n_rep', 1000, ...
                    'corr_type', 'Spearman', ...
                    'threshold', 0.05, ...
                    'side', 'right', ...
                    'sessions', []);
                % Run word boundaries analysis
                % sn_obj.test_s_vs_n_wordboundaries('S_condition_flag', 'sentence', ...
                %                   'N_condition_flag', 'nonword', ...
                %                   'n_rep', 1000, ...
                %                   'epoch_range', [-0.25 0.25], ...
                %                   'num_words', 12, ...
                %                   'do_plot', false);

                % Prepare output directory
                new_crunched_dir = save_dir;
                if ~exist(new_crunched_dir, 'dir')
                    mkdir(new_crunched_dir);
                end

                % Save result
                [~, name, ~] = fileparts(files(f).name);
                output_filename = fullfile(new_crunched_dir, sprintf('%s_ecog_sn_analyzed.mat', name));
                save(output_filename, 'sn_obj', '-v7.3');
                fprintf('  Saved analyzed object: %s', output_filename);
            end

            fprintf('Successfully processed %s', subject_id);
        catch ME
            fprintf('Error processing %s: %s', subject_id, ME.message);
            continue;
        end
    end
    fprintf('S vs N analysis completed.');
end
