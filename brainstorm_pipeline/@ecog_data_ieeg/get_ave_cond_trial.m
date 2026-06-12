function [output_tbl,cond_table]=get_ave_cond_trial(obj,ops)
    % Averages signal for given condition across all words.    
    arguments
        obj ecog_data_ieeg
        ops.words double = 1:12;
        ops.condition = [];
        ops.keep_trials = [];
    end
    % p = inputParser();
    % addParameter(p,'words',1:12);
    % addParameter(p,'condition',[]);
    % addParameter(p,'keep_trials',[])
    % parse(p, varargin{:});
    % ops = p.Results;
        
    func = @(x) cell2mat(permute(x,[3,2,1])); % format : electrode*trial_id*state/word
        
    % get trials condition
    if ops.condition
        condition_flag = ops.condition;
        cond_data = obj.get_cond_resp(condition_flag,'keep_trials',ops.keep_trials);
    else 
        if isempty(obj.trial_data)
            obj.make_trials();
        end
        condition_flag = 'all';
        cond_data = obj.trial_data;
    end
    
    cond_data_ave = obj.get_average(cond_data);
    word_data = obj.get_value(cond_data_ave,'key','word','type','contain');
        
    B = obj.combine_trial_cond(word_data);
    [keys,strings,values] = obj.get_columns(B);
    values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
    cond_table = cell2table(horzcat(condition_flag,{strings.string},values_comb),'VariableNames',B.Properties.VariableNames);

    % select how many words are selected : default 1:12 and average over words 
    % create the comparision between the two condition
    func_1 = @(x) x{1}(:,:,ops.words); % format : electrode*trial_id
    func_2 = @(x) nanmean(x,3); % format : electrode*trial_id
    
    B = cond_table;
    [keys,strings,values] = obj.get_columns(B);
    condition_ave = cellfun(@(X) func_2(func_1(values.(X))),values.Properties.VariableNames,'uni',false);
    
    output_tbl = cell2table(horzcat(condition_flag,strings.string,condition_ave),'VariableNames',B.Properties.VariableNames);
    
end 