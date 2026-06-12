
function saveUpdatedObject(obj)
    % Create the crunched folder if it doesn't exist
    crunchedFolder = fullfile(obj.crunched_file_path);
    if ~exist(crunchedFolder, 'dir')
        mkdir(crunchedFolder);
    end

    % Generate the filename
    filename = fullfile(crunchedFolder, [obj.subject '_' obj.experiment '_crunched_HG_ZScore.mat']);

    % Save the object
    save(filename, 'obj', '-v7.3');
    
    fprintf('Updated object saved as: %s\n', filename);
end