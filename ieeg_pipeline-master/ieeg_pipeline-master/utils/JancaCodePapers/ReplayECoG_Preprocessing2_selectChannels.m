%%%%%% Preprocessing2 - automatic channels selection %%%%%%
% JB Eichenlaub - Final version replay ECoG 01.12.17
% Intput:
%   - rawDataIni: Preprocessed data from Preprocessing1. fieldtrip structure - 1 cell per segment (empty if no data)
%   - nameChans
% Output [MGxx_Preprocessing2.mat]:
%   - rawData_LFP_autoSelec: fieldtrip structure - 1 cell per segment (empty if no data) after automatic channel selection
%   - detectionIEDs:         output from automatic IEDs detection (but after having deleted segments info)
%   - nameChans:             name of channels
%   - infoFile:              keep track of the main info (path, ...)
% Functions:
%   - automaticSpikeDetection_UsingJancaMethod
%   - selectChannels_UsingJancaMethod

%% Define ALL Patients and paths
clearvars
% allPatientsID = {'BW26';'BW27';'MG72';'MG76';'MG79';'MG83';'MG85';'MG88';'MG91'};
allPatientsID = {'MG72'};

pathProject = 'Y:\JBE\replayECoG\';
segmentsID  = {'sleepPre'; 'rest'; 'learning'; 'sleepPost'; 'retest'};
addPaths_ReplayECoG(pathProject); % add paths to Fiedltrip, NPMK and home-made functions

for p=1:numel(allPatientsID)
    
    clearvars -except allPatientsID pathProject segmentsID p
    patientID = allPatientsID(p);
    
    %% Load rawdata from Preprocessing1 (LFP data)
    pathFile = [pathProject 'data\' patientID{1} '\analysisSession\'];
    nameFile = [patientID{1} '_Preprocessing1.mat'];
    load([pathFile nameFile], 'rawDataIni', 'nameChans');
    
    %% Automatic detection, per channel, of interictal epileptic discharges (IEDs) using Janca et al., Brain Topogr, 2015
    detectionIEDs          = []; % output from Janca et al. script - 1 cell array per segment
    detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1'; % if you change "-dec 200" here, do not forget to change in selectChannels_Using ... below
    detectionIEDs.segments = [];
    
    for seg = 1:numel(segmentsID)
        if ~isempty(rawDataIni{seg})
            d  = []; d  = transpose(rawDataIni{seg}.trial{1}); % d  ... signal (time x channel)
            fs = []; fs = rawDataIni{seg}.fsample;             % fs ... sampling frequency (Hz)
            detectionIEDs.segments{seg} = automaticSpikeDetection_UsingJancaMethod(d, fs, detectionIEDs.settings);
        else
            detectionIEDs.segments{seg} = [];
        end
    end
    
    % From this automatic assessment, selection of the final pool of channels
    detectionIEDs.tableChanSelection = [];  % info regarding chan selection
    detectionIEDs.threshold          = 6.5; % channels with IEDs higher than threshold are removed
    
    detectionIEDs.tableChanSelection = selectChannels_UsingJancaMethod(detectionIEDs.segments, detectionIEDs.threshold, ...
        nameChans.Ini, detectionIEDs.settings, patientID, segmentsID);
    
    nameChans.rejectedAutomatic = setdiff(nameChans.Ini, detectionIEDs.tableChanSelection.nameChansSelected, 'stable');
    nameChans.selectedAutomatic = detectionIEDs.tableChanSelection.nameChansSelected;
    
    %% rawdata after exclusion of automatically rejected channels
    rawData_LFP_autoSelec = [];
    
    if ~isempty(nameChans.rejectedAutomatic)    % if few channels ARE REJECTED
        
        for seg = 1:numel(segmentsID)
            if ~isempty(rawDataIni{seg})
                cfg                        = [];
                cfg.channel                = nameChans.selectedAutomatic;
                rawData_LFP_autoSelec{seg} = ft_preprocessing(cfg, rawDataIni{seg});
            else
                rawData_LFP_autoSelec{seg} = [];
            end
        end
        
    else                                        % if NO rejection
        rawData_LFP_autoSelec = rawDataIni;
    end
    
    %% Save data
    % save info file to keep track
    infoFile             = [];
    infoFile.patientID   = patientID;
    infoFile.segmentsID  = segmentsID;
    infoFile.pathProject = pathProject;
    % but delete detectionIEDs.segments (really big)
    detectionIEDs.segments = [];
    % save
    pathOutput = [pathProject 'data\' patientID{1} '\analysisSession\'];
    fileOutput = [pathOutput patientID{1}, '_Preprocessing2.mat'];
    save(fileOutput, 'rawData_LFP_autoSelec', 'detectionIEDs' ,'nameChans', 'infoFile', '-v7.3');
    
end
