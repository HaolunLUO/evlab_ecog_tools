function output_d=combine_trial_cond(obj,input_d)
    % Combines separate trial tables (e.g., for a condition) into one 
    % table.
    
    all_keys = cellfun(@(x) x.key,input_d,'uni',false);

    % make sure the keys are the same for all trials 
    [X,Y] = ndgrid(1:numel(all_keys));
    Z = tril(true(numel(all_keys)),-1);
    assert(all(arrayfun(@(x,y) isequal(all_keys{x},all_keys{y}),X(Z),Y(Z))));
        
    % function to combine trial data tables
    func = @(x) cell2mat(reshape(vertcat(x),1,[]));

    output_d = table();

    % go through keys in trial data table
    for k=1:numel(all_keys{1})
        
        each_key = all_keys{1}{k};
        temp = (cellfun(@(x) x(ismember(x.key,each_key),:),input_d,'uni',false));
        cond_tbl = vertcat(temp{:});
        [~,strings,values] = obj.get_columns(cond_tbl);
        values_comb = cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
        temp_table = cell2table(horzcat(each_key,{strings.string},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
        
        output_d = [output_d; temp_table];

    end 
        
end