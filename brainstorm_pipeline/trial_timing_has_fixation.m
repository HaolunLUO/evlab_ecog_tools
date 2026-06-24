function tf = trial_timing_has_fixation(tbl)
% TRIAL_TIMING_HAS_FIXATION  True when row 1 is the MGH-style fixation epoch.
if ~istable(tbl) || ~ismember('key', tbl.Properties.VariableNames) || isempty(tbl)
    tf = false;
    return;
end
tf = strcmpi(strtrim(tbl.key{1}), 'fix');
end
