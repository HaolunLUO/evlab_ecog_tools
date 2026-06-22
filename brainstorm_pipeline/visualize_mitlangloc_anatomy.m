function visualize_mitlangloc_anatomy()
% VISUALIZE_MITLANGLOC_ANATOMY  MNI anatomy plots for other-lab MITLangloc.
%
%   Visualizes language-responsive electrode locations for the MGH / other-lab
%   MITLangloc cohort (broadband -> HG crunched files), aligned with settings
%   in compare_lang_profiles_mni.m (otherTaskType = 'MITLangloc').
%
%   Other-lab files:
%     <workingDir>/crunched/MITLangloc/<Subject>_<otherTaskType>_crunched_HG_ZScore.mat
%
%   Outputs per-subject MNI scatter/GIF and a cohort-level combined figure.
%
%   Run:
%     visualize_mitlangloc_anatomy

%% USER SETTINGS  (aligned with compare_lang_profiles_mni.m)
repoRoot   = fileparts(fileparts(mfilename('fullpath')));
workingDir = 'F:\seeg\luohong\analysisEV';

otherTaskType   = 'MITLangloc';
otherCohortName = 'MGH';
otherCrunchedDir = fullfile(workingDir, 'crunched', 'MITLangloc');

% Leave {} to auto-detect all subjects with finished analysis in otherCrunchedDir
otherSubjects = {};

signalType = 'bipolar';        % 'bipolar' | 'unipolar'
sigSource  = 's_vs_n';         % 's_vs_n' | 'wordwise'
auxSearchDirs = {workingDir, otherCrunchedDir, 'F:\iEEG_evlab'};

doPerSubject       = true;
doGroupPlot        = true;
plotAllElectrodes  = true;     % grey = all channels with MNI
doGIF              = false;
angle              = 270;
isPlotVisible      = false;

outputDir = fullfile(workingDir, 'output', 'MITLangloc', ...
    sprintf('anatomy_distribution_%s', lower(signalType)));

%% SETUP
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(scriptDir));
addpath(fullfile(repoRoot, 'ieeg_pipeline-master', 'ieeg_pipeline-master'));

if isempty(otherSubjects)
    otherSubjects = discover_broadband_subjects(otherCrunchedDir, otherTaskType);
end
otherSubjects = subjects_to_cellstr(otherSubjects);

if isempty(otherSubjects)
    error(['No other-lab subjects found in %s\n' ...
        'Expected: <Subject>_%s_crunched_HG_ZScore.mat\n' ...
        'Run run_mitlangloc_from_broadband_crunched first.'], ...
        otherCrunchedDir, otherTaskType);
end

if ~exist(outputDir, 'dir'); mkdir(outputDir); end

fprintf('\n=== MITLANGLOC ANATOMY DISTRIBUTION (%s) ===\n', otherCohortName);
fprintf('Task:        %s\n', otherTaskType);
fprintf('Directory:   %s\n', otherCrunchedDir);
fprintf('Subjects:    %s\n', strjoin(otherSubjects, ', '));
fprintf('Signal:      %s (%s)\n', signalType, sigSource);
fprintf('Output:      %s\n', outputDir);

groupRows = table();

for i = 1:numel(otherSubjects)
    subj = otherSubjects{i};
    subjDir = fullfile(outputDir, subj);
    if doPerSubject && ~exist(subjDir, 'dir')
        mkdir(subjDir);
    end

    try
        [sn_obj, meta] = load_sn_obj_for_comparison(subj, otherCrunchedDir, ...
            repoRoot, otherTaskType, auxSearchDirs, 'broadband');
    catch ME
        warning('Skipping %s: %s', subj, ME.message);
        continue;
    end

    nSig = sum(get_lang_sig_mask(sn_obj, signalType, sigSource));
    fprintf('\n%s (%s): %d sig %s channels | %s\n', ...
        subj, otherCohortName, nSig, signalType, meta.file);

    plotOut = struct('nSig', 0, 'sigCoords', zeros(0, 3), 'sigLabels', {{}});

    if doPerSubject
        plotOut = plot_lang_sig_mni_distribution(sn_obj, subjDir, ...
            'subjectName', subj, 'signalType', signalType, ...
            'sigSource', sigSource, 'plotAllElectrodes', plotAllElectrodes, ...
            'doGIF', doGIF, 'angle', angle, 'isPlotVisible', isPlotVisible);
    else
        plotOut = plot_lang_sig_mni_distribution(sn_obj, outputDir, ...
            'subjectName', subj, 'signalType', signalType, ...
            'sigSource', sigSource, 'saveFigure', false, ...
            'isPlotVisible', false);
    end

    for k = 1:plotOut.nSig
        newRow = table(string(subj), string(otherCohortName), string(plotOut.sigLabels{k}), ...
            plotOut.sigCoords(k, 1), plotOut.sigCoords(k, 2), plotOut.sigCoords(k, 3), ...
            'VariableNames', {'subject', 'cohort', 'label', 'mni_x', 'mni_y', 'mni_z'});
        groupRows = [groupRows; newRow]; %#ok<AGROW>
    end
end

if isempty(groupRows) || height(groupRows) == 0
    warning('No language-responsive channels with MNI coordinates found.');
    return;
end

T = groupRows;
writetable(T, fullfile(outputDir, sprintf('%s_lang_sig_mni_coords.csv', otherCohortName)));
fprintf('\nCohort total: %d sig channels with MNI across %d subjects\n', ...
    height(T), numel(unique(T.subject)));

if doGroupPlot
    plot_cohort_mni_distribution(T, otherSubjects, otherCohortName, ...
        otherCrunchedDir, repoRoot, otherTaskType, auxSearchDirs, ...
        signalType, sigSource, outputDir, angle, isPlotVisible);
end

fprintf('\nSaved anatomy figures to: %s\n', outputDir);

end


function plot_cohort_mni_distribution(T, subjects, cohortName, crunchedDir, ...
    repoRoot, taskType, auxSearchDirs, signalType, sigSource, outputDir, angle, isPlotVisible)

templateSn = [];
for i = 1:numel(subjects)
    try
        [sn, ~] = load_sn_obj_for_comparison(subjects{i}, crunchedDir, ...
            repoRoot, taskType, auxSearchDirs, 'broadband');
        if isprop(sn, 'anatomy') && isstruct(sn.anatomy) ...
                && isfield(sn.anatomy, 'template_brain') ...
                && ~isempty(sn.anatomy.template_brain)
            templateSn = sn;
            break;
        end
    catch
    end
end

cmap = lines(numel(subjects));
subjList = cellstr(unique(T.subject, 'stable'));

fig = figure('Color', 'w', 'Position', [80 80 1300 900], ...
    'Visible', ternary(isPlotVisible, 'on', 'off'));
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

if ~isempty(templateSn)
    try_plot_cortex_from_sn(ax, templateSn, 0.3);
    axis(ax, 'off');
else
    xlabel(ax, 'MNI X (mm)'); ylabel(ax, 'MNI Y (mm)'); zlabel(ax, 'MNI Z (mm)');
    grid(ax, 'on');
end

hLeg = gobjects(numel(subjList), 1);
for s = 1:numel(subjList)
    subj = subjList{s};
    mask = strcmp(string(T.subject), string(subj));
    coords = [T.mni_x(mask), T.mni_y(mask), T.mni_z(mask)];
    if isempty(coords)
        continue;
    end
    hLeg(s) = scatter3(ax, coords(:,1), coords(:,2), coords(:,3), ...
        100, cmap(s,:), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
end

axis(ax, 'equal');
view(ax, angle, 0);
if ~isempty(templateSn)
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
end
title(ax, sprintf('%s %s - Lang Sig MNI (%s, %s, n=%d)', ...
    cohortName, taskType, signalType, sigSource, height(T)));
legend(hLeg, subjList, 'Location', 'eastoutside', 'Interpreter', 'none');
drawnow;

pngName = fullfile(outputDir, sprintf('%s_cohort_langSigMNI_%s.png', cohortName, signalType));
saveas(fig, pngName);
close(fig);
fprintf('Saved cohort figure: %s\n', pngName);
end


function tf = try_plot_cortex_from_sn(ax, sn_obj, alphaVal)
tf = false;
try
    tmpl = sn_obj.anatomy.template_brain;
    if exist('plot3DModel', 'file') == 2 && isfield(tmpl, 'cortex')
        plot3DModel(ax, tmpl.cortex, [], alphaVal);
        colormap(ax, gray);
        tf = true;
        return;
    end
    if isfield(tmpl, 'cortex')
        cortex = tmpl.cortex;
        if isfield(cortex, 'faces') && isfield(cortex, 'vertices')
            F = cortex.faces; V = cortex.vertices;
        elseif isfield(cortex, 'tri') && isfield(cortex, 'vert')
            F = cortex.tri; V = cortex.vert;
        else
            return;
        end
        patch(ax, 'Faces', F, 'Vertices', V, ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', alphaVal);
        tf = true;
    end
catch
end
end


function subjects = discover_broadband_subjects(crunchedDir, taskType)
patterns = {
    fullfile(crunchedDir, sprintf('*%s_crunched_HG_ZScore.mat', taskType))
    fullfile(crunchedDir, sprintf('*%s_crunched_HG.mat', taskType))
    fullfile(crunchedDir, '*_MITLangloc_crunched_HG_ZScore.mat')
    fullfile(crunchedDir, '*_MITLangloc_crunched_HG.mat')
};

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
    if ~has_sn_obj_in_file(fullfile(d(i).folder, d(i).name))
        continue;
    end
    subj = subject_from_broadband_name(d(i).name, taskType);
    if ~isempty(subj)
        subjects{end+1} = subj; %#ok<AGROW>
    end
end
subjects = unique(subjects, 'stable');
end


function ok = has_sn_obj_in_file(matFile)
ok = false;
try
    vars = whos('-file', matFile);
    ok = any(strcmp({vars.name}, 'sn_obj'));
catch
end
end


function subj = subject_from_broadband_name(fname, taskType)
subj = '';
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


function out = subjects_to_cellstr(subjects)
if ischar(subjects)
    out = {subjects};
elseif isstring(subjects)
    out = cellstr(subjects(:))';
elseif iscell(subjects)
    out = cellfun(@char, subjects, 'UniformOutput', false);
    out = out(:)';
else
    error('Subject list must be char, string, or cell array.');
end
end


function v = ternary(cond, a, b)
if cond
    v = a;
else
    v = b;
end
end
