function [mask, probs] = get_grey_matter_channel_mask(sn_obj, signalType, minProb, varargin)
% GET_GREY_MATTER_CHANNEL_MASK  Grey-matter probability channel filter.
%
%   [mask, probs] = get_grey_matter_channel_mask(sn_obj, signalType, minProb, ...
%       'subject', subjectId, 'anatomyDir', anatomyDir)
%
%   Reads contact-level CAT12 tissue probabilities from
%   <anatomyDir>/<subject>.tsv (columns tissues_cat12, tissues_cat12_prob).
%   For bipolar channels, mask is true when at least one contact has
%   grey-matter probability > minProb.
%
%   Grey-matter probability:
%     tissues_cat12 == Gray  -> tissues_cat12_prob
%     tissues_cat12 == White -> 1 - tissues_cat12_prob
%     other tissues          -> 0

p = inputParser();
addParameter(p, 'subject', '');
addParameter(p, 'anatomyDir', '');
parse(p, varargin{:});
subject = char(string(p.Results.subject));
anatomyDir = p.Results.anatomyDir;

if nargin < 2 || isempty(signalType)
    signalType = 'bipolar';
end
if nargin < 3 || isempty(minProb)
    minProb = 0.5;
end
if isempty(subject) || isempty(anatomyDir)
    error('get_grey_matter_channel_mask:MissingInput', ...
        'Both ''subject'' and ''anatomyDir'' are required.');
end

anatomyTbl = load_subject_anatomy_tsv(subject, anatomyDir);

labels = channel_labels(sn_obj, signalType);
nChan = numel(labels);
mask = false(nChan, 1);
probs = nan(nChan, 1);

for i = 1:nChan
    if strcmpi(signalType, 'bipolar')
        [mask(i), probs(i)] = bipolar_grey_matter_ok(sn_obj, i, minProb, anatomyTbl);
    else
        [mask(i), probs(i)] = unipolar_grey_matter_ok(sn_obj, i, minProb, anatomyTbl);
    end
end
end


function labels = channel_labels(sn_obj, signalType)
if strcmpi(signalType, 'bipolar')
    labels = sn_obj.bip_ch_label(:);
else
    labels = sn_obj.elec_ch_label(:);
end
end


function [ok, prob] = unipolar_grey_matter_ok(sn_obj, uniIdx, minProb, anatomyTbl)
ok = false;
prob = NaN;

elecLabels = get_object_cellstr(sn_obj, 'elec_ch_label');
if uniIdx > numel(elecLabels)
    return;
end

prob = tsv_grey_matter_prob(anatomyTbl, elecLabels{uniIdx});
ok = isfinite(prob) && prob > minProb;
end


function [ok, prob] = bipolar_grey_matter_ok(sn_obj, bipIdx, minProb, anatomyTbl)
ok = false;
prob = NaN;

bipLabels = get_object_cellstr(sn_obj, 'bip_ch_label');
pairLabels = get_bipolar_pair_labels(sn_obj, bipIdx, bipLabels);
if isempty(pairLabels)
    return;
end

contactProbs = nan(numel(pairLabels), 1);
for i = 1:numel(pairLabels)
    contactProbs(i) = tsv_grey_matter_prob(anatomyTbl, pairLabels{i});
end
contactProbs = contactProbs(isfinite(contactProbs));
if isempty(contactProbs)
    return;
end

prob = max(contactProbs);
ok = any(contactProbs > minProb);
end


function prob = tsv_grey_matter_prob(anatomyTbl, contactLabel)
prob = NaN;

hit = find_anatomy_row(anatomyTbl, contactLabel);
if isempty(hit)
    return;
end

tissue = anatomyTbl.tissues_cat12(hit);
tissueProb = anatomyTbl.tissues_cat12_prob(hit);
prob = grey_matter_prob_from_tissue(tissue, tissueProb);
end


function prob = grey_matter_prob_from_tissue(tissue, tissueProb)
prob = NaN;
p = parse_percent_str(tissueProb);
if isnan(p)
    return;
end

tissueNorm = lower(strtrim(char(string(tissue))));
if strcmp(tissueNorm, 'gray') || strcmp(tissueNorm, 'grey')
    prob = p;
elseif strcmp(tissueNorm, 'white')
    prob = 1 - p;
else
    prob = 0;
end
end


function p = parse_percent_str(value)
p = NaN;
if isempty(value)
    return;
end
if ismissing(value)
    return;
end

s = strtrim(char(string(value)));
if isempty(s) || strcmpi(s, 'N/A') || strcmpi(s, 'NA')
    return;
end

s = strrep(s, '%', '');
p = str2double(s);
if ~isfinite(p)
    return;
end
if p > 1
    p = p / 100;
end
end


function hit = find_anatomy_row(anatomyTbl, label)
target = normalize_anatomy_label(label);
hit = find(anatomyTbl.channel_norm == target, 1, 'first');
end


function out = normalize_anatomy_label(labels)
out = upper(string(labels));
out = regexprep(out, '[^A-Z0-9]', '');
end


function pairLabels = get_bipolar_pair_labels(sn_obj, bipIdx, labels)
pairLabels = {};
if ~isempty(labels) && bipIdx <= numel(labels)
    parts = split_bipolar_label(labels{bipIdx});
    if numel(parts) >= 2
        pairLabels = parts(1:2);
        return;
    end
end
if isprop(sn_obj, 'bip_ch_label_valid') && ~isempty(sn_obj.bip_ch_label_valid)
    bvl = sn_obj.bip_ch_label_valid;
    if size(bvl, 2) >= 2 && bipIdx <= size(bvl, 1)
        pairLabels = {bvl{bipIdx, 1}, bvl{bipIdx, 2}};
    end
end
end


function labels = get_object_cellstr(obj, name)
labels = {};
if isprop(obj, name)
    val = obj.(name);
    if iscell(val) && ~isempty(val)
        labels = val;
    end
end
end
