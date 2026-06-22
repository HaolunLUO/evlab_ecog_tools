function [sigMask, info] = get_lang_sig_mask(sn_obj, signalType, sigSource)
% GET_LANG_SIG_MASK  Language-responsive channel mask for comparison scripts.
%
%   Bipolar (default for cross-cohort compare):
%     sig  -> sn_obj.s_vs_n_sig.bip_elec_data{1,1}
%     labels -> sn_obj.bip_ch_label  (same index order)
%
%   sigSource:
%     's_vs_n'   - trial-averaged Sentences vs N (test_s_vs_n), default
%     'wordwise' - per-word permutation + >=3 consecutive significant words

if nargin < 2 || isempty(signalType)
    signalType = 'bipolar';
end
if nargin < 3 || isempty(sigSource)
    sigSource = 's_vs_n';
end

info = struct('source', sigSource, 'signalType', signalType, ...
    'description', '');

switch lower(sigSource)
    case 'wordwise'
        sigMask = wordwise_sig_mask(sn_obj, signalType);
        info.description = 'wordwise (>=3 consecutive significant words)';
    case 's_vs_n'
        sigMask = svN_sig_mask(sn_obj, signalType);
        info.description = 's_vs_n_sig.bip_elec_data{1,1} / elec_data{1,1}';
    otherwise
        error('Unknown sigSource: %s (use ''s_vs_n'' or ''wordwise'')', sigSource);
end

end


function sigMask = svN_sig_mask(sn_obj, signalType)
if ~istable(sn_obj.s_vs_n_sig)
    error('sn_obj.s_vs_n_sig missing; run test_s_vs_n first.');
end

if strcmpi(signalType, 'bipolar')
    if ~ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
        error('No bipolar s_vs_n_sig.bip_elec_data in object.');
    end
    sigMask = logical(sn_obj.s_vs_n_sig.bip_elec_data{1, 1});
    validate_bipolar_label_alignment(sn_obj, sigMask);
else
    if ~ismember('elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
        error('No unipolar s_vs_n_sig.elec_data in object.');
    end
    sigMask = logical(sn_obj.s_vs_n_sig.elec_data{1, 1});
end
sigMask = sigMask(:);
end


function validate_bipolar_label_alignment(sn_obj, sigMask)
if ~isprop(sn_obj, 'bip_ch_label') || isempty(sn_obj.bip_ch_label)
    error('bip_ch_label required; must align with s_vs_n_sig.bip_elec_data{1,1}.');
end
if numel(sn_obj.bip_ch_label) ~= numel(sigMask)
    error(['bip_ch_label (%d) and s_vs_n_sig.bip_elec_data{1,1} (%d) ' ...
        'length mismatch for %s.'], numel(sn_obj.bip_ch_label), numel(sigMask), ...
        get_subject_id(sn_obj));
end
end


function id = get_subject_id(sn_obj)
id = 'unknown';
if isprop(sn_obj, 'subject') && ~isempty(sn_obj.subject)
    id = sn_obj.subject;
end
end


function sigMask = wordwise_sig_mask(sn_obj, signalType)
if ~isprop(sn_obj, 'langloc_wordwise') || isempty(sn_obj.langloc_wordwise) ...
        || ~isfield(sn_obj.langloc_wordwise, 'results')
    error('langloc_wordwise results missing; run test_s_vs_n_wordwise first.');
end

res = sn_obj.langloc_wordwise.results;
if strcmpi(signalType, 'bipolar')
    if ~isfield(res, 'bipolar') || ~isfield(res.bipolar, 'is_sig')
        error('No bipolar wordwise results in object.');
    end
    sigMask = logical(res.bipolar.is_sig(:));
    validate_bipolar_label_alignment(sn_obj, sigMask);
else
    if ~isfield(res, 'unipolar') || ~isfield(res.unipolar, 'is_sig')
        error('No unipolar wordwise results in object.');
    end
    sigMask = logical(res.unipolar.is_sig(:));
end
end
