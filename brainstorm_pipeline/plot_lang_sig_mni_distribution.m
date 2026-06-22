function out = plot_lang_sig_mni_distribution(sn_obj, outputDir, varargin)
% PLOT_LANG_SIG_MNI_DISTRIBUTION  Plot language-responsive channels in MNI.
%
%   Uses get_lang_channel_mni_coords (bip_ch_pos_mni, clinical info, or
%   anatomy mapping). Works on other-lab MITLangloc crunched files that lack
%   full VERA anatomy.template_brain.
%
%   plot_lang_sig_mni_distribution(sn_obj, outputDir, Name, Value, ...)
%
%   Name-Value (optional):
%     subjectName       - title / filename prefix (default 'subject')
%     signalType        - 'bipolar' (default) or 'unipolar'
%     sigSource         - 's_vs_n' (default) or 'wordwise'
%     plotAllElectrodes - show all channels with valid MNI (default false)
%     doGIF             - rotate and save GIF (default false)
%     angle             - azimuth degrees (default 270)
%     brainAlpha        - cortex transparency (default 0.35)
%     sigColor          - RGB for sig channels (default red)
%     allColor          - RGB for non-sig (default grey)
%     sigSize, allSize  - marker sizes (default 120, 25)
%     isPlotVisible     - logical (default true)
%     saveFigure        - logical (default true)
%
%   Returns struct out with fields sigCoords (Nx3), sigLabels, nSig.

p = inputParser();
addRequired(p, 'sn_obj');
addRequired(p, 'outputDir');
addParameter(p, 'subjectName', 'subject');
addParameter(p, 'signalType', 'bipolar');
addParameter(p, 'sigSource', 's_vs_n');
addParameter(p, 'plotAllElectrodes', false);
addParameter(p, 'doGIF', false);
addParameter(p, 'angle', 270);
addParameter(p, 'brainAlpha', 0.35);
addParameter(p, 'allColor', [0.75 0.75 0.75]);
addParameter(p, 'sigColor', [1 0 0]);
addParameter(p, 'allSize', 25);
addParameter(p, 'sigSize', 120);
addParameter(p, 'isPlotVisible', true);
addParameter(p, 'saveFigure', true);
parse(p, sn_obj, outputDir, varargin{:});
ops = p.Results;

if ops.saveFigure && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

[sigCoords, sigLabels, allCoords] = collect_mni_coords(sn_obj, ops.signalType, ops.sigSource);

out = struct();
out.sigCoords = sigCoords;
out.sigLabels = sigLabels;
out.nSig = size(sigCoords, 1);

if out.nSig == 0
    warning('plot_lang_sig_mni_distribution:NoSigMni', ...
        '%s: no significant %s channels with valid MNI.', ops.subjectName, ops.signalType);
    if ~ops.saveFigure
        return;
    end
end

if ops.doGIF
    angles = ops.angle:(ops.angle + 360);
else
    angles = ops.angle;
end

fig = figure('Color', 'w', 'Position', [100 100 1200 800], ...
    'Visible', ternary(ops.isPlotVisible, 'on', 'off'));
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

im = cell(numel(angles), 1);
hasCortex = try_plot_cortex(ax, sn_obj, ops.brainAlpha);

for i = 1:numel(angles)
    if i > 1
        cla(ax);
        hasCortex = try_plot_cortex(ax, sn_obj, ops.brainAlpha);
    end

    if ops.plotAllElectrodes && ~isempty(allCoords)
        scatter3(ax, allCoords(:,1), allCoords(:,2), allCoords(:,3), ...
            ops.allSize, ops.allColor, 'filled', 'MarkerEdgeColor', 'k');
    end

    if ~isempty(sigCoords)
        scatter3(ax, sigCoords(:,1), sigCoords(:,2), sigCoords(:,3), ...
            ops.sigSize, ops.sigColor, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end

    axis(ax, 'equal');
    if hasCortex
        axis(ax, 'off');
    else
        xlabel(ax, 'MNI X (mm)'); ylabel(ax, 'MNI Y (mm)'); zlabel(ax, 'MNI Z (mm)');
        grid(ax, 'on');
    end
    view(ax, angles(i), 0);
    if hasCortex
        camlight(ax, 'headlight');
        lighting(ax, 'gouraud');
    end
    title(ax, sprintf('%s - Lang Sig (%s, %s, n=%d)', ...
        ops.subjectName, ops.signalType, ops.sigSource, out.nSig));
    drawnow;

    frame = getframe(fig);
    im{i} = frame2im(frame);
end

if ops.saveFigure
    baseName = fullfile(outputDir, sprintf('%s_langSigMNI_%s', ...
        ops.subjectName, ops.signalType));
    if ops.doGIF
        gifName = [baseName '.gif'];
        for i = 1:numel(im)
            [A, map] = rgb2ind(im{i}, 256);
            if i == 1
                imwrite(A, map, gifName, 'gif', 'LoopCount', Inf, 'DelayTime', 0.05);
            else
                imwrite(A, map, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
            end
        end
        fprintf('Saved GIF: %s\n', gifName);
    else
        pngName = [baseName sprintf('_angle%d.png', ops.angle)];
        imwrite(im{1}, pngName);
        fprintf('Saved PNG: %s\n', pngName);
    end
end

close(fig);
end


function [sigCoords, sigLabels, allCoords] = collect_mni_coords(sn_obj, signalType, sigSource)
sigMask = get_lang_sig_mask(sn_obj, signalType, sigSource);
if strcmpi(signalType, 'bipolar')
    labels = sn_obj.bip_ch_label(:);
    nCh = numel(labels);
else
    labels = sn_obj.elec_ch_label(:);
    nCh = numel(labels);
end

allCoords = zeros(0, 3);
allLabels = {};
for ci = 1:nCh
    [mni, ok] = get_lang_channel_mni_coords(sn_obj, ci, signalType);
    if ok
        allCoords(end+1, :) = mni; %#ok<AGROW>
        allLabels{end+1} = labels{ci}; %#ok<AGROW>
    end
end

sigIdx = find(sigMask);
sigCoords = zeros(0, 3);
sigLabels = {};
for k = 1:numel(sigIdx)
    ci = sigIdx(k);
    [mni, ok] = get_lang_channel_mni_coords(sn_obj, ci, signalType);
    if ok
        sigCoords(end+1, :) = mni; %#ok<AGROW>
        sigLabels{end+1} = labels{ci}; %#ok<AGROW>
    end
end
end


function tf = try_plot_cortex(ax, sn_obj, alphaVal)
tf = false;
if ~isprop(sn_obj, 'anatomy') || isempty(sn_obj.anatomy) || ~isstruct(sn_obj.anatomy)
    return;
end
if ~isfield(sn_obj.anatomy, 'template_brain') || isempty(sn_obj.anatomy.template_brain)
    return;
end

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


function v = ternary(cond, a, b)
if cond
    v = a;
else
    v = b;
end
end
