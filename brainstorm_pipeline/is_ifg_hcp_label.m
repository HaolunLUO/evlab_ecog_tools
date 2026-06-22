function tf = is_ifg_hcp_label(label)
% IS_IFG_HCP_LABEL  True when an HCP-MMP / anatomy string is in IFG.
%
%   Matches HCP-MMP areas 44, 45, 47*, IFS*, FOP*, and common IFG synonyms.
%   Label format is case-insensitive; hemisphere prefixes are ignored.

if isempty(label) || (isstring(label) && strlength(label) == 0)
    tf = false;
    return;
end

norm = lower(regexprep(char(string(label)), '[^a-z0-9]', ''));
if isempty(norm) || strcmp(norm, 'na')
    tf = false;
    return;
end

if contains(norm, 'ifg') || contains(norm, 'inferiorfrontal') ...
        || contains(norm, 'parsopercularis') || contains(norm, 'parstriangularis') ...
        || contains(norm, 'parsorbitalis') || contains(norm, 'broca')
    tf = true;
    return;
end

ifgCodes = {'44', '45', '47l', '47m', '47s', '47r', 'p47r', 'a47r', ...
    'ifsa', 'ifsp', 'fop1', 'fop2', 'fop3', 'fop4', 'fop5'};
for i = 1:numel(ifgCodes)
    code = ifgCodes{i};
    if strcmp(norm, code) || endsWith(norm, code)
        tf = true;
        return;
    end
end

tf = false;
end
