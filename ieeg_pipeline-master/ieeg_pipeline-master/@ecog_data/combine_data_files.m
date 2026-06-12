function combine_data_files(obj)
    % Combines trial info from multiple data files
    %
    % This should be the last step in preprocessing after trial info is 
    % added if there is more than one data file.

    if isa(obj.trial_timing{1},'table') % should be a cell if multiple files
        fprintf(1,'No need to combine data files!\n')
        return
    end

    % output
    trial_timing_raw = [];
    trial_timing_dec = [];
    condition        = [];
    session          = [];
    events_table     = [];

    fprintf(1, '\n> Combining trial info from data files \n');
    fprintf(1,'[');

    for i=1:length(obj.stitch_index) % number of separate data files with signal
        
        samples_to_add_raw = obj.for_preproc.stitch_index_raw(i)-1;
        samples_to_add_dec = obj.for_preproc.stitch_index_dec(i)-1;
        trial_timing_raw_  = obj.for_preproc.trial_timing_raw{i}; % all trial timing tables for one file
        trial_timing_dec_  = obj.for_preproc.trial_timing_dec{i}; % all trial timing tables for one file

        % go through all trials in file
        for j=1:size(trial_timing_raw_,1)
            trial_timing_raw__ = trial_timing_raw_{j}; % one trial timing table
            trial_timing_dec__ = trial_timing_dec_{j}; % one trial timing table

            trial_timing_raw__(:,'start') = table(trial_timing_raw__.start + samples_to_add_raw); 
            trial_timing_raw__(:,'end')   = table(trial_timing_raw__.end + samples_to_add_raw); 
            trial_timing_dec__(:,'start') = table(trial_timing_dec__.start + samples_to_add_dec); 
            trial_timing_dec__(:,'end')   = table(trial_timing_dec__.end + samples_to_add_dec); 

            trial_timing_raw_{j} = trial_timing_raw__;
            trial_timing_dec_{j} = trial_timing_dec__;

            fprintf(1,'.');

        end

        trial_timing_raw = [trial_timing_raw; trial_timing_raw_];
        trial_timing_dec = [trial_timing_dec; trial_timing_dec_];

        % no need to update timing information
        condition    = [condition; obj.condition{i}];
        session      = [session; obj.session{i}];
        events_table = [events_table; obj.events_table{i}]; % based on within trial timing

    end

    fprintf(1,'] done\n');
        
    % add both versions of trial timing (raw & ds) 
    obj.for_preproc.trial_timing_raw = trial_timing_raw;
    obj.for_preproc.trial_timing_dec = trial_timing_dec;

    % update obj.trial_timing (trial timing for current version of signal)
    if obj.sample_freq == obj.for_preproc.decimation_freq % downsampling in preproc
        obj.trial_timing = trial_timing_dec;
    else % no downsampling in preproc
        obj.trial_timing = trial_timing_raw;
    end

    obj.condition    = condition;
    obj.session      = session;
    obj.events_table = events_table;

end 
