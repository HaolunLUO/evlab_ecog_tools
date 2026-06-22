function diag_compare_pipeline()
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

subject = 'Subject01';
crunchedDir = 'F:\seeg\luohong\analysisEV';
taskType = 'MITSWJNTask';
auxSearchDirs = {crunchedDir, fullfile(crunchedDir, 'crunched', 'MITLangloc'), 'F:\iEEG_evlab'};

fprintf('=== load_sn_obj_for_comparison ===\n');
try
    [sn_obj, meta] = load_sn_obj_for_comparison(subject, crunchedDir, repoRoot, ...
        taskType, auxSearchDirs, 'brainstorm');
    fprintf('loaded: %s\n', meta.file);
catch ME
    fprintf('LOAD FAILED: %s\n', ME.message);
    return;
end

sig = find(sn_obj.s_vs_n_sig.bip_elec_data{1, 1});
fprintf('sig bipolar: %d\n', numel(sig));
for k = 1:min(3, numel(sig))
    [m, ok] = get_lang_channel_mni_coords(sn_obj, sig(k), 'bipolar');
    fprintf('  %s ok=%d %s\n', sn_obj.bip_ch_label{sig(k)}, ok, mat2str(m, 3));
end

fprintf('=== build_lang_channel_profiles (no timecourses) ===\n');
args = {'cohort', 'local', 'subject', meta.subject, ...
    'signalType', 'bipolar', 'words', 1:12, ...
    'S_condition', 'Sentences', 'sigOnly', true, ...
    'use_odd_for_inference', false};
try
    P = build_lang_channel_profiles(sn_obj, args{:});
    fprintf('profiles with MNI: %d\n', numel(P));
catch ME
    fprintf('BUILD FAILED: %s\n', ME.message);
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
end
