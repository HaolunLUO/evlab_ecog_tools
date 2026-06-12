%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AVERAGE TRIAL DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function input_d_ave=get_average(obj,input_d)
        % Averages trial data signal. Returns cell array of trial timing
        % tables with the average signal value for each event (e.g., word)
        % and signal type (e.g., unipolar, bipolar)
        %
        % NOTE - it does not operate on obj.trial_data.
        input_d_ave = input_d;

        % go through trials of input
        for k=1:size(input_d,1)
            B = input_d{k};

            [keys,strings,values] = obj.get_columns(B);
            
            values_ave = varfun( @(x) cellfun(@(y) nanmean(y,2), x,'uni',false), values,'OutputFormat','table');
            values_ave.Properties.VariableNames = values.Properties.VariableNames;
            input_d_ave{k} = [keys,strings,values_ave];

        end

    end