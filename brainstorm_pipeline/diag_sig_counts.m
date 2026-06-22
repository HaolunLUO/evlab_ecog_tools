function diag_sig_counts()
%DIAG_SIG_COUNTS  Compare significance definitions for MGH crunched files.

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

crunchedDir = 'F:\seeg\luohong\analysisEV\crunched\MITLangloc';
outputDir   = 'F:\seeg\luohong\analysisEV\output\MITLangloc';
subs = {'AMC092', 'BJH006', 'BJH007', 'BJH011'};

fprintf('\n=== MGH significance counts ===\n');
fprintf('%-8s  %-6s  %-6s  %-6s  %-6s  %-6s  %-6s\n', ...
    'Subject', 'uniSn', 'bipSn', 'bipWw', 'bipMni', 'nSigGr', 'nSigBGr');
fprintf('%s\n', repmat('-', 1, 60));

for i = 1:numel(subs)
    subj = subs{i};
    matFile = fullfile(crunchedDir, sprintf('%s_MITLangloc_crunched_HG_ZScore.mat', subj));
    if ~isfile(matFile)
        fprintf('%s: file missing\n', subj);
        continue;
    end
    sn = load(matFile, 'sn_obj');
    sn = sn.sn_obj;

    sigU = sum(logical(sn.s_vs_n_sig.elec_data{1, 1}));
    sigB = 0;
    if ismember('bip_elec_data', sn.s_vs_n_sig.Properties.VariableNames)
        sigB = sum(logical(sn.s_vs_n_sig.bip_elec_data{1, 1}));
    end

    sigW = NaN;
    if isprop(sn, 'langloc_wordwise') && ~isempty(sn.langloc_wordwise) ...
            && isfield(sn.langloc_wordwise, 'results') ...
            && isfield(sn.langloc_wordwise.results, 'bipolar')
        sigW = sum(sn.langloc_wordwise.results.bipolar.is_sig);
    end

    sigBipMni = 0;
    if sigB > 0
        idx = find(logical(sn.s_vs_n_sig.bip_elec_data{1, 1}));
        for k = 1:numel(idx)
            [~, ok] = get_lang_channel_mni_coords(sn, idx(k), 'bipolar');
            sigBipMni = sigBipMni + ok;
        end
    end

    nGr = NaN;
    nBGr = NaN;
    grFile = fullfile(outputDir, sprintf('%s_MITLangloc_groupResult.mat', subj));
    if isfile(grFile)
        G = load(grFile, 'groupResult');
        if isfield(G.groupResult, 'nSig')
            nGr = G.groupResult.nSig;
        end
        if isfield(G.groupResult, 'nSigBip')
            nBGr = G.groupResult.nSigBip;
        end
    end

    fprintf('%-8s  %6d  %6d  %6.0f  %6d  %6.0f  %6.0f\n', ...
        subj, sigU, sigB, sigW, sigBipMni, nGr, nBGr);
end

fprintf(['\nColumns:\n' ...
    '  uniSn   = test_s_vs_n unipolar (trial-averaged S vs N, p<0.05 right)\n' ...
    '  bipSn   = test_s_vs_n bipolar (same test, used by compare_lang_profiles_mni)\n' ...
    '  bipWw   = test_s_vs_n_wordwise bipolar (>=3 consecutive words sig)\n' ...
    '  bipMni  = bipSn channels with valid MNI\n' ...
    '  nSigGr  = groupResult.nSig (unipolar, saved by run script)\n' ...
    '  nSigBGr = groupResult.nSigBip (if saved)\n']);
end
