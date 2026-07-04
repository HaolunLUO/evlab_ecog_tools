function [obj, taskConfig] = map_and_resolve_task_conditions(obj, taskConfig)
% MAP_AND_RESOLVE_TASK_CONDITIONS  Apply conditionMap and resolve S/N flags.
%
%   Handles Brainstorm event labels such as 'sentence_tp1' by stripping an
%   optional '_tpN' suffix before lookup in taskConfig.conditionMap.

fprintf('Original conditions: %s\n', strjoin(unique(obj.condition), ', '));

if isfield(taskConfig, 'conditionMap') && ~isempty(taskConfig.conditionMap)
    condMap = taskConfig.conditionMap;
    mapKeys = condMap.keys;
    for i = 1:numel(obj.condition)
        key = lookup_condition_map_key(obj.condition{i}, mapKeys);
        if ~isempty(key)
            obj.condition{i} = condMap(key);
        end
    end
end

avail = unique(obj.condition, 'stable');
fprintf('Mapped conditions:   %s\n', strjoin(avail, ', '));

sCandidates = build_condition_candidates(taskConfig, 'S_condition');
nCandidates = build_condition_candidates(taskConfig, 'N_condition');

taskConfig.S_condition = resolve_condition_name(avail, sCandidates, 'S');
taskConfig.N_condition = resolve_condition_name(avail, nCandidates, 'N');

fprintf('Resolved contrast: %s vs %s\n', taskConfig.S_condition, taskConfig.N_condition);

end


function key = lookup_condition_map_key(raw, mapKeys)
key = '';
if any(strcmp(raw, mapKeys))
    key = raw;
    return;
end

base = regexprep(raw, '_tp\d+$', '', 'ignorecase');
if any(strcmp(base, mapKeys))
    key = base;
    return;
end

hit = find(strcmpi(raw, mapKeys), 1);
if ~isempty(hit)
    key = mapKeys{hit};
    return;
end

hit = find(strcmpi(base, mapKeys), 1);
if ~isempty(hit)
    key = mapKeys{hit};
end
end


function candidates = build_condition_candidates(taskConfig, fieldName)
candidates = {};
if isfield(taskConfig, fieldName) && ~isempty(taskConfig.(fieldName))
    candidates{end + 1} = taskConfig.(fieldName); %#ok<AGROW>
end

if ~isfield(taskConfig, 'conditionMap') || isempty(taskConfig.conditionMap)
    return;
end

condMap = taskConfig.conditionMap;
mapKeys = condMap.keys;
target = '';
if isfield(taskConfig, fieldName)
    target = taskConfig.(fieldName);
end

for i = 1:numel(mapKeys)
    mapped = condMap(mapKeys{i});
    if ~isempty(target) && strcmpi(mapped, target)
        candidates{end + 1} = mapKeys{i}; %#ok<AGROW>
        candidates{end + 1} = mapped; %#ok<AGROW>
    end
end

candidates = unique(candidates, 'stable');
end
