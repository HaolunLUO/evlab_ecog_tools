function [keys,strings,values]=get_columns(obj,B)
    % Pulls out columns of provided trial data table.

    keys = B(:,ismember(B.Properties.VariableNames,'key'));
    strings = B(:,ismember(B.Properties.VariableNames,'string'));
    values = B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
    
end
