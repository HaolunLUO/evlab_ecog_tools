function n = trial_timing_max_word_index(tbl)
% TRIAL_TIMING_MAX_WORD_INDEX  Largest word index from keys like word_1..word_N.
n = 0;
if ~istable(tbl) || ~ismember('key', tbl.Properties.VariableNames)
    return;
end
for i = 1:height(tbl)
    tok = regexp(strtrim(tbl.key{i}), '^word_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        n = max(n, str2double(tok{1}));
    end
end
end
