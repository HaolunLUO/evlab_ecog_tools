function tbl = prepend_fixation_row(tbl, fixationDurationSec, sampleRate)
% PREPEND_FIXATION_ROW  Insert MGH-style fix row before word_1.
if trial_timing_has_fixation(tbl)
    return;
end
w1Row = find(strcmp(tbl.key, 'word_1'), 1);
if isempty(w1Row)
    error('prepend_fixation_row: word_1 not found in trial_timing keys.');
end
word1Start = tbl.start(w1Row);
fixEnd = word1Start - 1;
fixStart = word1Start - round(fixationDurationSec * sampleRate);
fixStart = max(1, fixStart);
if fixEnd < fixStart
    fixStart = max(1, fixEnd);
end
fixRow = table({'fix'}, {''}, fixStart, fixEnd, ...
    'VariableNames', tbl.Properties.VariableNames);
tbl = [fixRow; tbl];
end
