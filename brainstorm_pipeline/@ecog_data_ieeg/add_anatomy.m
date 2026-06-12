function add_anatomy(obj,anatomy_path,varargin)
    % Adds anatomy files to object, including a mapping from the channels
    % in the anatomy file to the channels in the object.
    %
    % Assumes particular naming convention for files 
    % TODO: change to be more general

    p = inputParser();
    addParameter(p,'veraFolder',true); % if anatomy_path refers to directory with vera folders or directly to anatomical files
    addParameter(p,'subdirName','VERA_'); % subject name is appended to end of this prefix
    addParameter(p,'templateName','cvs_avg35_inMNI152'); % same folder and file name (with .mat appended)
    parse(p, varargin{:});
    ops = p.Results;

    % subdirectory name (if applicable)
    if ops.veraFolder
        folder = [ops.subdirName obj.subject filesep];
        template_folder = [ops.templateName filesep];
    else
        folder = ''; % path refers to folder with all anatomical files
        template_folder = '';
    end

    % subject-specific space
    filename = [anatomy_path folder obj.subject '_brain.mat'];
    subject_space = load(filename);
    obj.anatomy.subject_space = subject_space;

    % MNI space
    filename = [anatomy_path folder obj.subject '_MNI_brain.mat'];
    mni_space = load(filename);
    obj.anatomy.mni_space = mni_space;

    % template brain
    filename = [anatomy_path template_folder ops.templateName '.mat'];
    template_brain = load(filename);
    obj.anatomy.template_brain = template_brain;

    % add mapping from anatomical files to object
    [mapping,labels] = obj.channel_mapping_anatomical(subject_space);
    obj.anatomy.mapping = mapping;
    obj.anatomy.labels = labels;
    
    % add hemisphere labels
    hemisphere = cell(size(subject_space.tala.electrodes,1),1);
    right_idxs = find(subject_space.tala.electrodes(:,1) > 0);
    hemisphere(right_idxs,1) = {'right'};
    left_idxs = find(subject_space.tala.electrodes(:,1) < 0);
    hemisphere(left_idxs,1) = {'left'};
    obj.anatomy.hemisphere = hemisphere;

end