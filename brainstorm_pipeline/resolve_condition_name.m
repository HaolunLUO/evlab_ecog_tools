function name = resolve_condition_name(avail, candidates, label)
% RESOLVE_CONDITION_NAME  Match a condition flag to sn_obj.condition names.
%
%   name = resolve_condition_name(avail, {'Sentences','SENTENCES'}, 'S')
%
%   Matching is case-insensitive; hyphens, spaces, and underscores are ignored
%   when comparing (e.g. 'Nonword-lists' matches 'NONWORDS').

if nargin < 3 || isempty(label)
    label = 'condition';
end

avail = cellstr(avail);
if isempty(avail)
    error('No conditions available in object for %s.', label);
end

candidates = cellstr(candidates);
for c = 1:numel(candidates)
    hit = find(strcmpi(strtrim(avail), strtrim(candidates{c})), 1, 'first');
    if ~isempty(hit)
        name = avail{hit};
        return;
    end
end

norm = @(s) regexprep(lower(strtrim(string(s))), '[^a-z0-9]', '');
for c = 1:numel(candidates)
    target = norm(candidates{c});
    if strlength(target) == 0
        continue;
    end
    for i = 1:numel(avail)
        if norm(avail{i}) == target
            name = avail{i};
            return;
        end
    end
end

error('Could not resolve %s condition. Tried: %s. Available: %s', ...
    label, strjoin(candidates, ', '), strjoin(avail, ', '));

end
