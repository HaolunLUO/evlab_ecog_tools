function [cond_data, cond_id]=get_cond_resp(obj,condition,ops)
    % Extracts trial data of given condition. 
    arguments
        obj ecog_data
        condition char
        ops.keep_trials = []; 
    end
    % p = inputParser();
    % addParameter(p,'keep_trials',[]);
    % parse(p, varargin{:});
    % ops = p.Results;

    if isempty(obj.trial_data)
        obj.make_trials();
    end
    
    cond_id = obj.get_cond_id(condition,'keep_trials',ops.keep_trials);

    cond_data = obj.trial_data(cond_id);
    
end