function [cond_id]=get_cond_id(obj,condition,ops)
    % Extracts trial data of given condition. 
    arguments
        obj ecog_data_ieeg
        condition = [];
        ops.keep_trials = []; 
    end
    % p = inputParser();
    % addParameter(p,'keep_trials',[]);
    % parse(p, varargin{:});
    % ops = p.Results;

    
    if ~isempty(ops.keep_trials)
        assert(length(obj.condition)==length(ops.keep_trials),'Logical array of trials to keep has incorrect dimensions');
        cond_id = (cell2mat(arrayfun(@(x) (strcmp(obj.condition{x},condition) && ops.keep_trials(x)),1:length(obj.condition),'UniformOutput',false)));
    else
        cond_id = (cell2mat(arrayfun(@(x) strcmp(obj.condition{x},condition),1:length(obj.condition),'UniformOutput',false)));
    end

    
    
end