function tbl = load_subject_anatomy_tsv(subject, anatomyDir)
% LOAD_SUBJECT_ANATOMY_TSV  Load per-contact anatomy from <subject>.tsv.
%
%   Expected columns include Channel, tissues_cat12, tissues_cat12_prob.
%   File path: <anatomyDir>/<subject>.tsv

if nargin < 2 || isempty(anatomyDir)
    error('anatomyDir is required.');
end
if isempty(subject)
    error('subject is required.');
end

subject = char(string(subject));
tsvPath = fullfile(anatomyDir, [subject '.tsv']);
if ~isfile(tsvPath)
    error('Anatomy TSV not found: %s', tsvPath);
end

opts = detectImportOptions(tsvPath, 'FileType', 'text', 'Delimiter', '\t');
opts = setvaropts(opts, 'Channel', 'Type', 'string');
tbl = readtable(tsvPath, opts);

if ~ismember('Channel', tbl.Properties.VariableNames)
    error('Column ''Channel'' missing in %s', tsvPath);
end
if ~ismember('tissues_cat12', tbl.Properties.VariableNames) ...
        || ~ismember('tissues_cat12_prob', tbl.Properties.VariableNames)
    error(['Columns ''tissues_cat12'' and ''tissues_cat12_prob'' required in %s'], tsvPath);
end

tbl.channel_norm = normalize_anatomy_label(tbl.Channel);
end


function out = normalize_anatomy_label(labels)
out = upper(string(labels));
out = regexprep(out, '[^A-Z0-9]', '');
end
