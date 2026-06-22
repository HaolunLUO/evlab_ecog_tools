function compare_lang_profiles_mni()
% COMPARE_LANG_PROFILES_MNI  Match lang channels by MNI and compare HG profiles.
%
%   Pairs language-responsive bipolar channels (s_vs_n_sig.bip_elec_data{1,1},
%   labels in bip_ch_label) at similar MNI coordinates between
%   your Brainstorm pipeline data and other-lab MITLangloc results, then
%   compares bipolar (or unipolar) HG response profiles.
%
%   Local (Brainstorm) naming from complete_mit_pipeline_brainstorm.m:
%     <workingDir>/<Subject>_<localTaskType>_crunched.mat
%     legacy: <Subject>_MITLangloc_crunched.mat
%
%   Other-lab naming from run_mitlangloc_from_broadband_crunched.m:
%     <workingDir>/crunched/MITLangloc/<Subject>_<otherTaskType>_crunched_HG_ZScore.mat
%
%   Run:
%     compare_lang_profiles_mni

%% USER SETTINGS  (aligned with complete_mit_pipeline_brainstorm.m)
repoRoot   = fileparts(fileparts(mfilename('fullpath')));
workingDir = 'F:\seeg\luohong\analysisEV';

% File naming differs between cohorts (see complete_mit_pipeline_brainstorm.m):
localTaskType  = 'MITSWJNTask';   % Brainstorm: Subject01_MITSWJNTask_crunched.mat
otherTaskType  = 'MITLangloc';    % Other lab:  AMC092_MITLangloc_crunched_HG_ZScore.mat

% --- Local (Brainstorm pipeline) ---
localCohortName  = 'local';
localCrunchedDir = workingDir;   % same as complete_mit_pipeline taskCrunchedFile folder

% --- Other lab (broadband -> HG pipeline) ---
otherCohortName  = 'MGH';
otherCrunchedDir = fullfile(workingDir, 'crunched', 'MITLangloc');

% Leave empty {} to auto-detect all subjects with finished analysis in each folder
localSubjects = {'Subject07','Subject09','Subject06','Subject11','Subject04','Subject01'};
otherSubjects = {'AMC092','BJH007','BJH006','BJH011'};

maxMniDistMm  = 12;
signalType    = 'bipolar';       % 'bipolar' | 'unipolar'
sigSource     = 's_vs_n';        % s_vs_n_sig.bip_elec_data{1,1} + bip_ch_label
words         = 1:12;
S_condition   = 'Sentences';
N_condition   = '';              % auto: Jabberwocky / Nonword-lists

outputDir = fullfile(workingDir, 'output', 'MITLangloc', ...
    sprintf('cross_cohort_mni_%s', lower(signalType)));
doPlots   = true;
auxSearchDirs = {workingDir, otherCrunchedDir, 'F:\iEEG_evlab'};

%% SETUP
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

if isempty(localSubjects)
    localSubjects = discover_cohort_subjects(localCrunchedDir, localTaskType, 'brainstorm');
end
if isempty(otherSubjects)
    otherSubjects = discover_cohort_subjects(otherCrunchedDir, otherTaskType, 'broadband');
end

localSubjects = subjects_to_cellstr(localSubjects);
otherSubjects = subjects_to_cellstr(otherSubjects);

if isempty(localSubjects)
    error(['No local subjects found in %s\n' ...
        'Expected: <Subject>_%s_crunched.mat (or legacy <Subject>_MITLangloc_crunched.mat)\n' ...
        'with sn_obj from complete_mit_pipeline_brainstorm.'], ...
        localCrunchedDir, localTaskType);
end
if isempty(otherSubjects)
    error(['No other-lab subjects found in %s\n' ...
        'Expected: <Subject>_%s_crunched_HG_ZScore.mat (run run_mitlangloc_from_broadband_crunched).'], ...
        otherCrunchedDir, otherTaskType);
end

if ~exist(outputDir, 'dir'); mkdir(outputDir); end

fprintf('\n=== CROSS-COHORT LANG CHANNEL COMPARISON (MNI) ===\n');
fprintf('Signal type: %s\n', signalType);
fprintf('Sig source:  %s\n', sigSource);
fprintf('Local dir:   %s\n', localCrunchedDir);
fprintf('Other dir:   %s\n', otherCrunchedDir);
fprintf('Local (n=%d):  %s\n', numel(localSubjects), strjoin(localSubjects, ', '));
fprintf('Other (n=%d):  %s\n', numel(otherSubjects), strjoin(otherSubjects, ', '));
fprintf('Max MNI distance: %.1f mm\n', maxMniDistMm);

profilesLocal = load_cohort_profiles(localSubjects, localCrunchedDir, localCohortName, ...
    repoRoot, localTaskType, signalType, sigSource, words, S_condition, N_condition, ...
    auxSearchDirs, 'brainstorm');
profilesOther = load_cohort_profiles(otherSubjects, otherCrunchedDir, otherCohortName, ...
    repoRoot, otherTaskType, signalType, sigSource, words, S_condition, N_condition, ...
    auxSearchDirs, 'broadband');

fprintf('Local lang channels w/ MNI:  %d\n', numel(profilesLocal));
fprintf('Other lang channels w/ MNI: %d\n', numel(profilesOther));

if isempty(profilesLocal) || isempty(profilesOther)
    diagnose_missing_mni(localSubjects, localCrunchedDir, localTaskType, ...
        otherSubjects, otherCrunchedDir, otherTaskType, repoRoot, signalType, auxSearchDirs);
end

pairs = match_lang_channels_by_mni(profilesLocal, profilesOther, maxMniDistMm);
fprintf('Matched pairs (<= %.1f mm): %d\n', maxMniDistMm, numel(pairs));

if isempty(pairs)
    warning('No MNI-matched pairs found. Try increasing maxMniDistMm or check bipolar MNI coords.');
    return;
end

pairTable = pairs_to_table(pairs);
writetable(pairTable, fullfile(outputDir, 'mni_matched_pairs.csv'));
save(fullfile(outputDir, 'cross_cohort_comparison.mat'), ...
    'pairs', 'profilesLocal', 'profilesOther', 'pairTable', ...
    'localSubjects', 'otherSubjects', 'signalType', ...
    'localTaskType', 'otherTaskType', 'maxMniDistMm', '-v7.3');

disp(pairTable(:, {'local_subject', 'local_label', 'other_subject', 'other_label', 'distMm'}));

if doPlots
    plot_matched_profiles(pairs, outputDir, words, signalType);
end

fprintf('\nSaved results to: %s\n', outputDir);

end


function out = subjects_to_cellstr(subjects)
% SUBJECTS_TO_CELLSTR  Normalize subject IDs for strjoin and {s} indexing.
if ischar(subjects)
    out = {subjects};
elseif isstring(subjects)
    out = cellstr(subjects(:))';
elseif iscell(subjects)
    out = cellfun(@char, subjects, 'UniformOutput', false);
    out = out(:)';
else
    error('compare_lang_profiles_mni:InvalidSubjects', ...
        'Subject list must be a char vector, string array, or cell array.');
end
end


function profiles = load_cohort_profiles(subjects, crunchedDir, cohortName, repoRoot, ...
    taskType, signalType, sigSource, words, S_condition, N_condition, auxSearchDirs, cohortKind)

profiles = struct([]);
for s = 1:numel(subjects)
    try
        [sn_obj, meta] = load_sn_obj_for_comparison(subjects{s}, crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, cohortKind);
        sigMask = get_lang_sig_mask(sn_obj, signalType, sigSource);
        nSig = sum(sigMask);
        args = {'cohort', cohortName, 'subject', meta.subject, ...
            'signalType', signalType, 'sigSource', sigSource, ...
            'words', words, 'S_condition', S_condition, 'sigOnly', true, ...
            'use_odd_for_inference', false};
        if ~isempty(N_condition)
            args = [args, {'N_condition', N_condition}]; %#ok<AGROW>
        end
        P = build_lang_channel_profiles(sn_obj, args{:});
        fprintf('  %s (%s): %d sig %s (%s), %d w/ MNI\n  %s\n', ...
            meta.subject, cohortName, nSig, signalType, sigSource, numel(P), meta.file);
        profiles = append_profiles(profiles, P);
    catch ME
        warning('compare_lang_profiles_mni:SkipSubject', ...
            'Skipping %s: %s', subjects{s}, ME.message);
    end
end

end


function profiles = append_profiles(profiles, P)
if isempty(P)
    return;
end
if isempty(profiles)
    profiles = P(:);
else
    profiles = [profiles(:); P(:)]; %#ok<AGROW>
end
end


function subjects = discover_cohort_subjects(crunchedDir, taskType, cohortKind)
% DISCOVER_COHORT_SUBJECTS  Auto-find subjects with sn_obj analysis files.

if ~isfolder(crunchedDir)
    subjects = {};
    return;
end

switch lower(cohortKind)
    case 'brainstorm'
        patterns = {
            fullfile(crunchedDir, sprintf('*%s_crunched.mat', taskType))
            fullfile(crunchedDir, '*_MITLangloc_crunched.mat')
        };
    case 'broadband'
        patterns = {
            fullfile(crunchedDir, sprintf('*%s_crunched_HG_ZScore.mat', taskType))
            fullfile(crunchedDir, sprintf('*%s_crunched_HG.mat', taskType))
            fullfile(crunchedDir, '*_MITLangloc_crunched_HG_ZScore.mat')
            fullfile(crunchedDir, '*_MITLangloc_crunched_HG.mat')
        };
    otherwise
        error('Unknown cohortKind: %s', cohortKind);
end

d = [];
for p = 1:numel(patterns)
    d = [d; dir(patterns{p})]; %#ok<AGROW>
end

if isempty(d)
    subjects = {};
    return;
end

[~, ia] = unique({d.name}, 'stable');
d = d(ia);

subjects = {};
for i = 1:numel(d)
    if ~has_sn_obj(fullfile(d(i).folder, d(i).name))
        continue;
    end
    subj = subject_from_filename(d(i).name, taskType, cohortKind);
    if ~isempty(subj)
        subjects{end+1} = subj; %#ok<AGROW>
    end
end
subjects = unique(subjects, 'stable');
end


function ok = has_sn_obj(matFile)
ok = false;
try
    vars = whos('-file', matFile);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function subj = subject_from_filename(fname, taskType, cohortKind)
subj = '';
switch lower(cohortKind)
    case 'brainstorm'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched\.mat$'], 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
            return;
        end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
        end
    case 'broadband'
        tok = regexp(fname, ['^(.+)_' regexptranslate('escape', taskType) '_crunched_HG(_ZScore)?\.mat$'], 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
            return;
        end
        tok = regexp(fname, '^(.+)_MITLangloc_crunched_HG(_ZScore)?\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            subj = tok{1};
        end
end
end


function T = pairs_to_table(pairs)
n = numel(pairs);
localSubj = strings(n, 1);
localLab  = strings(n, 1);
localMni  = zeros(n, 3);
otherSubj = strings(n, 1);
otherLab  = strings(n, 1);
otherMni  = zeros(n, 3);
distMm    = zeros(n, 1);
localP    = zeros(n, 1);
otherP    = zeros(n, 1);

for i = 1:n
    a = pairs(i).profileA;
    b = pairs(i).profileB;
    localSubj(i)  = string(a.subject);
    localLab(i)   = string(a.label);
    localMni(i,:) = a.mni;
    otherSubj(i)  = string(b.subject);
    otherLab(i)   = string(b.label);
    otherMni(i,:) = b.mni;
    distMm(i)     = pairs(i).distMm;
    localP(i)     = a.p_ratio;
    otherP(i)     = b.p_ratio;
end

T = table(localSubj, localLab, localMni(:,1), localMni(:,2), localMni(:,3), ...
    otherSubj, otherLab, otherMni(:,1), otherMni(:,2), otherMni(:,3), ...
    distMm, localP, otherP, ...
    'VariableNames', {'local_subject', 'local_label', ...
    'local_mni_x', 'local_mni_y', 'local_mni_z', ...
    'other_subject', 'other_label', ...
    'other_mni_x', 'other_mni_y', 'other_mni_z', ...
    'distMm', 'local_p_ratio', 'other_p_ratio'});
end


function plot_matched_profiles(pairs, outputDir, words, signalType)
nPlot = min(numel(pairs), 24);
nCol = 3;
nRow = ceil(nPlot / nCol);

fig = figure('Color', 'w', 'Position', [50 50 1400 900], 'Visible', 'off');
tiled = tiledlayout(fig, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:nPlot
    ax = nexttile(tiled);
    a = pairs(i).profileA;
    b = pairs(i).profileB;

    if ~isempty(a.word_ave_S) && ~isempty(b.word_ave_S)
        tA = profile_word_time_axis(a);
        tB = profile_word_time_axis(b);
        plot(ax, tA, a.word_ave_S, 'b', 'LineWidth', 1.2); hold(ax, 'on');
        plot(ax, tA, a.word_ave_N, 'b--', 'LineWidth', 1);
        plot(ax, tB, b.word_ave_S, 'r', 'LineWidth', 1.2);
        plot(ax, tB, b.word_ave_N, 'r--', 'LineWidth', 1);
        xlabel(ax, 'Time (s)'); ylabel(ax, 'HG (z)');
        titleStr = sprintf('%s %s \\leftrightarrow %s %s (%.1f mm)', ...
            a.subject, a.label, b.subject, b.label, pairs(i).distMm);
        if numel(tA) ~= numel(tB)
            titleStr = sprintf('%s\n(%d vs %d samples)', titleStr, numel(tA), numel(tB));
        end
        title(ax, titleStr, 'Interpreter', 'none', 'FontSize', 9);
        legend(ax, {sprintf('%s S', a.cohort), sprintf('%s N', a.cohort), ...
            sprintf('%s S', b.cohort), sprintf('%s N', b.cohort)}, ...
            'Location', 'best', 'FontSize', 7);
        grid(ax, 'on');
    else
        text(ax, 0.5, 0.5, 'No word-average data', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
    end
end

saveas(fig, fullfile(outputDir, sprintf('matched_pair_word_averages_%s.png', signalType)));
close(fig);

fig2 = figure('Color', 'w', 'Position', [50 50 1400 900], 'Visible', 'off');
tiled2 = tiledlayout(fig2, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nPlot
    ax = nexttile(tiled2);
    a = pairs(i).profileA;
    b = pairs(i).profileB;
    if ~isempty(a.word_ave_S) && ~isempty(b.word_ave_S)
        tA = profile_word_time_axis(a);
        tB = profile_word_time_axis(b);
        diffA = profile_s_minus_n(a);
        diffB = profile_s_minus_n(b);
        plot(ax, tA, diffA, 'b', 'LineWidth', 1.3); hold(ax, 'on');
        plot(ax, tB, diffB, 'r', 'LineWidth', 1.3);
        yline(ax, 0, 'k:');
        xlabel(ax, 'Time (s)'); ylabel(ax, 'S - N (z)');
        title(ax, sprintf('%s vs %s (%.1f mm)', a.label, b.label, pairs(i).distMm), ...
            'Interpreter', 'none', 'FontSize', 9);
        legend(ax, {a.cohort, b.cohort}, 'Location', 'best');
        grid(ax, 'on');
    end
end
saveas(fig2, fullfile(outputDir, sprintf('matched_pair_S_minus_N_%s.png', signalType)));
close(fig2);

fprintf('Saved profile comparison figures (%s).\n', signalType);
end


function t = profile_word_time_axis(profile)
nSamp = size(profile.word_ave_S, 2);
fs = profile.sample_freq;
if isempty(fs) || fs <= 0
    fs = 1;
end
t = (0:nSamp - 1) / fs;
end


function diff = profile_s_minus_n(profile)
n = min(numel(profile.word_ave_S), numel(profile.word_ave_N));
diff = profile.word_ave_S(1:n) - profile.word_ave_N(1:n);
end


function diagnose_missing_mni(localSubjects, localDir, localTaskType, ...
    otherSubjects, otherDir, otherTaskType, repoRoot, signalType, auxSearchDirs)

fprintf('\n--- MNI diagnostics (why profiles may be empty) ---\n');

if ~isempty(localSubjects)
    subj = localSubjects{1};
    try
        [sn, meta] = load_sn_obj_for_comparison(subj, localDir, repoRoot, ...
            localTaskType, auxSearchDirs, 'brainstorm');
        report_sn_mni_status(sn, signalType, meta.file);
    catch ME
        fprintf('Local sample %s: %s\n', subj, ME.message);
    end
end

if ~isempty(otherSubjects)
    subj = otherSubjects{1};
    try
        [sn, meta] = load_sn_obj_for_comparison(subj, otherDir, repoRoot, ...
            otherTaskType, auxSearchDirs, 'broadband');
        report_sn_mni_status(sn, signalType, meta.file);
    catch ME
        fprintf('Other sample %s: %s\n', subj, ME.message);
    end
end

fprintf(['Tips:\n' ...
    '  Local:  Subject_%s_crunched.mat (Brainstorm pipeline).\n' ...
    '  Other:  Subject_%s_crunched_HG_ZScore.mat (broadband pipeline).\n'], ...
    localTaskType, otherTaskType);
end


function report_sn_mni_status(sn, signalType, filePath)
fprintf('\nFile: %s\n', filePath);

if isprop(sn, 'condition') && ~isempty(sn.condition)
    fprintf('  conditions: %s\n', strjoin(unique(sn.condition, 'stable'), ', '));
    try
        sName = resolve_condition_name(unique(sn.condition, 'stable'), ...
            {'Sentences', 'SENTENCES', 'sentence'}, 'S');
        nName = resolve_condition_name(unique(sn.condition, 'stable'), ...
            {'Jabberwocky', 'JABBERWOCKY', 'Nonword-lists', 'NONWORDS', 'Nonwords'}, 'N');
        fprintf('  resolved S/N: %s vs %s\n', sName, nName);
    catch ME
        fprintf('  condition resolve failed: %s\n', ME.message);
    end
end

if strcmpi(signalType, 'bipolar')
    sig = find(get_lang_sig_mask(sn, signalType, 's_vs_n'));
    nPos = 0;
    if isprop(sn, 'bip_ch_pos_mni') && ~isempty(sn.bip_ch_pos_mni)
        nPos = size(sn.bip_ch_pos_mni, 1);
    end
    fprintf('  sig bipolar: %d | bip_ch_pos_mni rows: %d | bip labels: %d\n', ...
        numel(sig), nPos, numel(sn.bip_ch_label));
    if numel(sig) > 0
        [m, ok] = get_lang_channel_mni_coords(sn, sig(1), 'bipolar');
        fprintf('  first sig bipolar MNI ok=%d (%s): %s\n', ok, sn.bip_ch_label{sig(1)}, mat2str(m,3));
        parts = split_bipolar_label(sn.bip_ch_label{sig(1)});
        if numel(parts) >= 2
            fprintf('  parsed contacts: %s | %s\n', parts{1}, parts{2});
        end
    end
else
    sig = find(get_lang_sig_mask(sn, 'unipolar', 's_vs_n'));
    nPos = 0;
    if isprop(sn, 'elec_ch_pos_mni') && ~isempty(sn.elec_ch_pos_mni)
        nPos = size(sn.elec_ch_pos_mni, 1);
    end
    hasAnat = isprop(sn, 'anatomy') && ~isempty(sn.anatomy);
    fprintf('  sig unipolar: %d | elec_ch_pos_mni rows: %d | anatomy: %d\n', ...
        numel(sig), nPos, hasAnat);
end
end
