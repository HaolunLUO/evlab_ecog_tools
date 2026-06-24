function span = trial_timing_word_sample_span(tbl, words)
% TRIAL_TIMING_WORD_SAMPLE_SPAN  [start end] sample indices spanning word keys.
if nargin < 2 || isempty(words)
    words = 1:trial_timing_max_word_index(tbl);
end
wStart = NaN;
wEnd = NaN;
for w = words(:)'
    row = find(strcmp(tbl.key, sprintf('word_%d', w)), 1);
    if isempty(row)
        continue;
    end
    if isnan(wStart)
        wStart = tbl.start(row);
    end
    wEnd = tbl.end(row);
end
if isnan(wStart) || isnan(wEnd)
    error('trial_timing_word_sample_span: no word rows found for indices %s', mat2str(words));
end
span = [wStart, wEnd];
end
