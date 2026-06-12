function [locs,ch_tags,ch_num]=extract_shanks(obj)
        
    ch_num = cellfun(@str2num,extract(obj.elec_ch_label,digitsPattern)); % nums
    try
        ch_tags = extract(obj.elec_ch_label,lettersPattern); % letters from label
    catch
        ch_tags = cellfun(@(x) regexprep(x, '\d+(?:_(?=\d))?', ''),...
            obj.elec_ch_label, 'UniformOutput', false);
    end
       
    % ch_tags = cellfun(@(x) [x '_'],ch_tags,'uni',false); % append underscore 
    [tags,~,~] = unique(ch_tags,'stable');
    locs = cellfun(@(x) find(ismember(ch_tags,x))',tags,'uni',false);
end