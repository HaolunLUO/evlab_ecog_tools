function [obj, taskConfig, report] = normalize_mitlangloc_events(obj, taskConfig)
% NORMALIZE_MITLANGLOC_EVENTS  Align conditions/events for MITLangloc analysis.
%
%   [obj, taskConfig, report] = normalize_mitlangloc_events(obj, taskConfig)
%
%   - Fills obj.condition from events_table when missing
%   - Applies taskConfig.conditionMap aliases
%   - Resolves S/N condition names against what is actually in the file
%     (MGH/MITLangloc often uses Sentences vs Jabberwocky, not Nonword-lists)

report = struct();
report.changedCondition = false;
report.changedS = false;
report.changedN = false;

if isempty(obj.condition) && isprop(obj, 'events_table') && istable(obj.events_table) ...
        && ismember('condition', obj.events_table.Properties.VariableNames)
    c = obj.events_table.condition;
    if iscell(c) && ~isempty(c) && iscell(c{1})
        obj.condition = cellfun(@(x) strtrim(char(x{1})), c, 'UniformOutput', false);
    else
        obj.condition = cellstr(c);
    end
    report.changedCondition = true;
    fprintf('Populated obj.condition from events_table (%d trials).\n', numel(obj.condition));
end

if isempty(obj.condition)
    error(['obj.condition is empty and could not be inferred from events_table. ' ...
        'Cannot run S-vs-N localization.']);
end

if isfield(taskConfig, 'conditionMap') && ~isempty(taskConfig.conditionMap)
    condMap = taskConfig.conditionMap;
    mapKeys = condMap.keys;
    for i = 1:numel(obj.condition)
        raw = obj.condition{i};
        if isKey(condMap, raw)
            obj.condition{i} = condMap(raw);
            continue;
        end
        hit = find(strcmpi(raw, mapKeys), 1);
        if ~isempty(hit)
            obj.condition{i} = condMap(mapKeys{hit});
        end
    end
end

avail = unique(obj.condition, 'stable');
report.availableConditions = avail;
fprintf('Trial conditions in file: %s\n', strjoin(avail, ', '));

if ~isempty(obj.trial_timing)
    tt = obj.trial_timing{1};
    if istable(tt) && ismember('key', tt.Properties.VariableNames)
        report.sampleTimingKeys = tt.key;
        nWordKeys = sum(contains(tt.key, 'word', 'IgnoreCase', true));
        fprintf('Trial 1 timing: %d rows, %d with ''word'' in key.\n', height(tt), nWordKeys);
        if nWordKeys == 0
            warning(['No trial_timing keys contain ''word''. get_ave_cond_trial expects ' ...
                'keys like word_1..word_12.']);
        end
    end
end

[taskConfig, report] = resolve_langloc_contrast(taskConfig, avail, obj.condition, report);

end


function [taskConfig, report] = resolve_langloc_contrast(taskConfig, avail, allConditions, report)

sAliases = {'Sentences', 'sentence', 'SENTENCES', 'S'};
nAliasesJabber = {'Jabberwocky', 'jabberwocky', 'J'};
nAliasesLists  = {'Nonword-lists', 'non-word lists', 'nonword', 'NONWORDS', 'nonword lists'};

sResolved = resolve_condition_name(taskConfig.S_condition, avail, sAliases);
nResolved = resolve_condition_name(taskConfig.N_condition, avail, [nAliasesLists, nAliasesJabber]);

if isempty(sResolved)
    error('S condition ''%s'' not found. Available: %s', ...
        taskConfig.S_condition, strjoin(avail, ', '));
end

if isempty(nResolved)
    nResolved = resolve_condition_name('Jabberwocky', avail, nAliasesJabber);
    if isempty(nResolved)
        nResolved = resolve_condition_name('Nonword-lists', avail, nAliasesLists);
    end
    if isempty(nResolved)
        error(['N condition ''%s'' not found. Available: %s\n' ...
            'For MGH/MITLangloc crunched files, set N_condition = ''Jabberwocky''.'], ...
            taskConfig.N_condition, strjoin(avail, ', '));
    end
    fprintf(['NOTE: N condition ''%s'' not in file; using ''%s'' instead ' ...
        '(standard MGH MITLangloc contrast).\n'], taskConfig.N_condition, nResolved);
    report.changedN = true;
end

if ~strcmp(taskConfig.S_condition, sResolved)
    fprintf('NOTE: S condition resolved to ''%s''.\n', sResolved);
    report.changedS = true;
end
if ~strcmp(taskConfig.N_condition, nResolved) && ~report.changedN
    fprintf('NOTE: N condition resolved to ''%s''.\n', nResolved);
    report.changedN = true;
end

taskConfig.S_condition = sResolved;
taskConfig.N_condition = nResolved;
report.S_condition = sResolved;
report.N_condition = nResolved;

fprintf('S vs N contrast: %s (%d trials) vs %s (%d trials)\n', ...
    sResolved, sum(strcmp(allConditions, sResolved)), ...
    nResolved, sum(strcmp(allConditions, nResolved)));

% Drop optional W/J conditions when this file only has S+N (typical MGH MITLangloc).
optionalFields = {'W_condition', 'J_condition'};
for i = 1:numel(optionalFields)
    fn = optionalFields{i};
    if isfield(taskConfig, fn) && ~isempty(taskConfig.(fn)) ...
            && ~any(strcmp(avail, taskConfig.(fn)))
        fprintf('NOTE: %s ''%s'' not in file; plots will use S vs N only.\n', ...
            fn, taskConfig.(fn));
        taskConfig.(fn) = '';
    end
end

end


function name = resolve_condition_name(requested, avail, aliases)
name = '';
if any(strcmp(avail, requested))
    name = requested;
    return;
end
for i = 1:numel(aliases)
    hit = find(strcmpi(avail, aliases{i}), 1);
    if ~isempty(hit)
        name = avail{hit};
        return;
    end
end
end
