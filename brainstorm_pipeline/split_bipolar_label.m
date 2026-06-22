function parts = split_bipolar_label(bipLabel)
% SPLIT_BIPOLAR_LABEL  Split a bipolar name into two contact labels.
%
%   'M5-M4' and 'LPC_2-LPC_1' -> two contact strings (split on hyphen only).

parts = {};
if isempty(bipLabel) || ~ischar(bipLabel) && ~isstring(bipLabel)
    return;
end

bipLabel = char(strtrim(bipLabel));
tok = regexp(bipLabel, '^(.+)-(.+)$', 'tokens', 'once');
if isempty(tok)
    return;
end

parts = {strtrim(tok{1}), strtrim(tok{2})};

end
