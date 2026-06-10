% PLOT_ANATOMY_UTILS  Anatomy visualization helpers for the Brainstorm pipeline.
%
% Standalone helper functions called by complete_mit_pipeline_brainstorm.m.
%
%   plot_sig_electrodes_anatomy  - scatter sig electrodes on brain mesh
%   get_electrode_coords         - extract coordinate array from anatomy struct
%   plot_template_cortex         - render the template brain surface


function plot_sig_electrodes_anatomy(sn_obj, outputDir, varargin)
% PLOT_SIG_ELECTRODES_ANATOMY  Plot significant/ROI electrodes on brain.
%
%   plot_sig_electrodes_anatomy(sn_obj, outputDir, Name, Value, ...)
%
%   Name-Value Parameters (all optional):
%     subjectName       - string for figure title / filename
%     space             - 'mni' (default) or 'subject'
%     signalType        - 'unipolar' (default) or 'bipolar'
%     plotAllElectrodes - logical (default false)
%     doGIF             - logical, create animated GIF (default false)
%     angle             - starting azimuth in degrees (default 270)
%     brainAlpha        - cortex transparency (default 0.35)
%     allColor          - RGB for non-sig electrodes (default grey)
%     sigColor          - RGB for sig electrodes (default red)
%     allSize           - marker size for all (default 25)
%     sigSize           - marker size for sig (default 120)
%     isPlotVisible     - logical (default true)

p = inputParser();
addRequired(p, 'sn_obj');
addRequired(p, 'outputDir');
addParameter(p, 'subjectName', 'subject');
addParameter(p, 'space', 'mni');
addParameter(p, 'signalType', 'unipolar');
addParameter(p, 'plotAllElectrodes', false);
addParameter(p, 'doGIF', false);
addParameter(p, 'angle', 270);
addParameter(p, 'brainAlpha', 0.35);
addParameter(p, 'allColor', [0.75 0.75 0.75]);
addParameter(p, 'sigColor', [1 0 0]);
addParameter(p, 'allSize', 25);
addParameter(p, 'sigSize', 120);
addParameter(p, 'isPlotVisible', true);
parse(p, sn_obj, outputDir, varargin{:});
ops = p.Results;

if ~exist(outputDir, 'dir'); mkdir(outputDir); end

if ~istable(sn_obj.s_vs_n_sig) || ~ismember('elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
    error('sn_obj.s_vs_n_sig is missing or not in expected format.');
end

% Choose unipolar or bipolar significant channels
if strcmp(ops.signalType, 'bipolar') && ...
        ismember('bip_elec_data', sn_obj.s_vs_n_sig.Properties.VariableNames)
    sigChans = find(sn_obj.s_vs_n_sig.bip_elec_data{1});
else
    sigChans = find(sn_obj.s_vs_n_sig.elec_data{1});
end

coords = get_electrode_coords(sn_obj, ops.space);

% Map data channel indices to anatomy electrode indices
mappedCells = sn_obj.anatomy.mapping(sigChans);
ok = ~cellfun(@isempty, mappedCells);
mappedIdx = cell2mat(mappedCells(ok));
mappedIdx = mappedIdx(mappedIdx >= 1 & mappedIdx <= size(coords,1));
mappedIdx = unique(mappedIdx(:));
sigCoords = coords(mappedIdx, :);

if ops.doGIF
    angles = ops.angle:(ops.angle+360);
else
    angles = ops.angle;
end

fig = figure('Color','w', 'Position',[100 100 1200 800], 'Visible', ops.isPlotVisible);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

im = cell(numel(angles),1);

for i = 1:numel(angles)
    cla(ax);

    plot_template_cortex(ax, sn_obj, ops.brainAlpha);

    if ops.plotAllElectrodes
        scatter3(ax, coords(:,1), coords(:,2), coords(:,3), ...
            ops.allSize, ops.allColor, 'filled', 'MarkerEdgeColor','k');
    end

    if ~isempty(sigCoords)
        scatter3(ax, sigCoords(:,1), sigCoords(:,2), sigCoords(:,3), ...
            ops.sigSize, ops.sigColor, 'filled', 'MarkerEdgeColor','k', 'LineWidth', 1.5);
    end

    axis(ax, 'equal'); axis(ax, 'off');
    view(ax, angles(i), 0);
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    title(ax, sprintf('%s - ROI/Sig Electrodes (%s space, %s)', ...
          ops.subjectName, upper(ops.space), ops.signalType));
    drawnow;

    frame = getframe(fig);
    im{i} = frame2im(frame);
end

baseName = fullfile(outputDir, sprintf('%s_sigElecs_%s_%s', ops.subjectName, ops.space, ops.signalType));

if ops.doGIF
    gifName = [baseName '.gif'];
    for i = 1:numel(im)
        [A,map] = rgb2ind(im{i}, 256);
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

close(fig);
end


function coords = get_electrode_coords(sn_obj, space)
% GET_ELECTRODE_COORDS  Extract electrode XYZ from anatomy struct.
space = lower(space);
switch space
    case 'mni'
        S = sn_obj.anatomy.mni_space;
    case 'subject'
        S = sn_obj.anatomy.subject_space;
    otherwise
        error('space must be ''mni'' or ''subject''.');
end

if isfield(S, 'tala') && isfield(S.tala, 'electrodes')
    coords = S.tala.electrodes;
elseif isfield(S, 'vera_mat_minimal') && isfield(S.vera_mat_minimal, 'tala') ...
        && isfield(S.vera_mat_minimal.tala, 'electrodes')
    coords = S.vera_mat_minimal.tala.electrodes;
else
    error('Could not find electrode coordinates in sn_obj.anatomy.%s_space', space);
end
end


function plot_template_cortex(ax, sn_obj, alphaVal)
% PLOT_TEMPLATE_CORTEX  Render brain surface on axis ax.
tmpl = sn_obj.anatomy.template_brain;

% Prefer the custom plot3DModel helper if available
if exist('plot3DModel', 'file') == 2 && isfield(tmpl, 'cortex')
    plot3DModel(ax, tmpl.cortex, [], alphaVal);
    colormap(ax, gray);
    return;
end

if isfield(tmpl, 'cortex')
    cortex = tmpl.cortex;
else
    error('template_brain missing cortex field and plot3DModel not available.');
end

if isfield(cortex, 'faces') && isfield(cortex, 'vertices')
    F = cortex.faces; V = cortex.vertices;
elseif isfield(cortex, 'tri') && isfield(cortex, 'vert')
    F = cortex.tri; V = cortex.vert;
else
    error('Unknown cortex mesh format. Check your anatomy struct.');
end

patch(ax, 'Faces', F, 'Vertices', V, ...
    'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', alphaVal);
end
