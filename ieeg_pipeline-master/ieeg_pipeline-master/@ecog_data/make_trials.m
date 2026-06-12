function make_trials(obj)
    % TODO - description

    if ~isempty(obj.bip_elec_data)
        trial_keys = {'key','string','elec_data','bip_elec_data'};
    else % no bipolar data
        trial_keys = {'key','string','elec_data'};
    end

    fprintf(1, '\n> Cutting signal into trial data ... \n');
    fprintf(1,'[');

    % go through each trial timing table
    for k = 1:size(obj.trial_timing,1)
                trial_time_tbl = obj.trial_timing{k};

        trial_elec_data = arrayfun(@(x) obj.elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
    
        if ~isempty(obj.bip_elec_data)
            trial_bip_elec_data = arrayfun(@(x) obj.bip_elec_data(:,trial_time_tbl(x,:).start:trial_time_tbl(x,:).end),1:size(trial_time_tbl,1),'uni',false)';
            obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,trial_bip_elec_data,'VariableNames',trial_keys);
        else % no bipolar data 
            obj.trial_data{k,1} = table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data,'VariableNames',trial_keys);
        end

        fprintf(1,'.');

    end
end