function run_langloc_report(obj, taskType, taskConfig, reportOutputDir, subjectName, useV2)
% RUN_LANGLOC_REPORT  Generate the langloc PDF via generateReportLangloc*.
%
%   run_langloc_report(obj, taskType, taskConfig, reportOutputDir, subjectName, useV2)

langlocLikeTasks = {'MITSWJNTask', 'MITLangloc'};
if ~ismember(taskType, langlocLikeTasks)
    fprintf('Skipping langloc report for task type: %s\n', taskType);
    return;
end

if isempty(obj.stats) || ~isstruct(obj.stats) ...
        || ~isfield(obj.stats, 'sig_hg_channel') ...
        || ~isfield(obj.stats, 'time_series')
    error(['obj.stats is incomplete. Run Step 7.5 ' ...
        '(extract_significant_channel / extract_time_significance) first.']);
end

reportObj = prepare_obj_for_langloc_report(obj, taskType, taskConfig, reportOutputDir);
reportName = sprintf('%s_%s_langloc_report', subjectName, taskType);
pdfPath = fullfile(reportOutputDir, [reportName, '.pdf']);

if exist(pdfPath, 'file')
    delete(pdfPath);
end

try
    if useV2
        generateReportLangloc_v2(reportObj, reportName);
    else
        generateReportLangloc(reportObj, reportName);
    end
catch ME
    if exist(pdfPath, 'file')
        fInfo = dir(pdfPath);
        if fInfo.bytes < 1024
            delete(pdfPath);
            fprintf('Removed incomplete PDF (%d bytes): %s\n', fInfo.bytes, pdfPath);
        end
    end
    rethrow(ME);
end

fInfo = dir(pdfPath);
if isempty(fInfo) || fInfo.bytes < 1024
    error('Langloc report PDF is missing or empty: %s', pdfPath);
end

fprintf('Langloc report saved (%d KB): %s\n', round(fInfo.bytes / 1024), pdfPath);
end


function reportObj = prepare_obj_for_langloc_report(obj, taskType, taskConfig, reportOutputDir)

reportObj = obj;
reportObj.experiment = 'MITLangloc';
reportObj.crunched_file_path = reportOutputDir;

condReportMap = build_report_condition_map(taskType, taskConfig);
condKeys = keys(condReportMap);
for i = 1:numel(reportObj.condition)
    c = reportObj.condition{i};
    matchIdx = find(strcmp(condKeys, c), 1);
    if ~isempty(matchIdx)
        reportObj.condition{i} = condReportMap(condKeys{matchIdx});
    end
end

reportObj.events_table = ensure_events_table(reportObj);

if isfield(taskConfig, 'nWordPositions') && ~isempty(taskConfig.nWordPositions)
    reportObj.report_numWords = taskConfig.nWordPositions;
elseif isfield(taskConfig, 'testWords') && ~isempty(taskConfig.testWords)
    reportObj.report_numWords = numel(taskConfig.testWords);
end
end


function condReportMap = build_report_condition_map(taskType, taskConfig)

switch taskType
    case 'MITLangloc'
        condReportMap = containers.Map( ...
            {taskConfig.S_condition, taskConfig.N_condition}, ...
            {'Sentences', 'Jabberwocky'});
    case 'MITSWJNTask'
        condReportMap = containers.Map( ...
            {'SENTENCES', 'NONWORDS'}, ...
            {'Sentences', 'Jabberwocky'});
    otherwise
        condReportMap = containers.Map( ...
            {taskConfig.S_condition, taskConfig.N_condition}, ...
            {'Sentences', 'Jabberwocky'});
end
end


function events_table = ensure_events_table(obj)

nTrials = numel(obj.condition);
if nTrials == 0
    error('prepare_obj_for_langloc_report: obj.condition is empty.');
end

if isprop(obj, 'events_table') && istable(obj.events_table) ...
        && height(obj.events_table) == nTrials ...
        && ismember('accuracy', obj.events_table.Properties.VariableNames) ...
        && ismember('RT', obj.events_table.Properties.VariableNames)
    events_table = obj.events_table;
    return;
end

events_table = table(ones(nTrials, 1), nan(nTrials, 1), ...
    'VariableNames', {'accuracy', 'RT'});
end
