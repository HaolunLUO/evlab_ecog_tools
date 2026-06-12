 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function output_d=get_value(obj,input_d,ops)
    % Extracts events with label ops.key (e.g., 'word') from trial
    % data. Returns cell array of trial timing tables with only 
    % events of interest.
    %
    % NOTE - this does not operate directly on obj.trial_data.
    arguments
        obj ecog_data
        input_d 
        ops.key = 'word';
        ops.type = 'match';
    end
    % p =inputParser();
    % addParameter(p, 'key', 'word');
    % addParameter(p, 'type', 'match'); % match or contain
    % parse(p, varargin{:});
    % ops = p.Results;

    if strcmp(ops.type,'match')
        func = @(x,y) ismember(x,y);
    else
        func = @(x,y) contains(x,y);
    end

    output_d = input_d;

    % go through trials of input
    for k = 1:size(input_d,1)
        B = input_d{k};
        output_d{k} = B(func(B.key,ops.key),:);
    end

end