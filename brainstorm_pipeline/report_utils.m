% REPORT_UTILS  Helpers for syncing analysis results and generating langloc PDFs.
%
% Called by complete_mit_pipeline_brainstorm.m.

function obj = sync_obj_from_sn_obj(obj, sn_obj)
% SYNC_OBJ_FROM_SN_OBJ  Copy processed data and stats from sn_obj onto obj.
%
%   obj = sync_obj_from_sn_obj(obj, sn_obj)
%
%   After Step 7.5+, all preprocessing and engine stats live on sn_obj.
%   This keeps the saved `obj` variable in the crunched .mat aligned with
%   the analysis object.

obj.elec_data     = sn_obj.elec_data;
obj.bip_elec_data = sn_obj.bip_elec_data;
obj.stitch_index  = sn_obj.stitch_index;
obj.sample_freq   = sn_obj.sample_freq;
obj.trial_data    = sn_obj.trial_data;

obj.trial_timing = sn_obj.trial_timing;
obj.condition    = sn_obj.condition;
obj.session      = sn_obj.session;

if isprop(sn_obj, 'events_table')
    obj.events_table = sn_obj.events_table;
end
if isprop(sn_obj, 'for_preproc')
    obj.for_preproc = sn_obj.for_preproc;
end

obj.elec_ch_with_IED      = sn_obj.elec_ch_with_IED;
obj.elec_ch_with_noise    = sn_obj.elec_ch_with_noise;
obj.elec_ch_user_deselect = sn_obj.elec_ch_user_deselect;
obj.elec_ch_clean         = sn_obj.elec_ch_clean;
obj.elec_ch_valid         = sn_obj.elec_ch_valid;

obj.bip_ch           = sn_obj.bip_ch;
obj.bip_ch_label     = sn_obj.bip_ch_label;
obj.bip_ch_valid     = sn_obj.bip_ch_valid;
obj.bip_ch_grp       = sn_obj.bip_ch_grp;
obj.bip_ch_label_grp = sn_obj.bip_ch_label_grp;

if isprop(sn_obj, 'anatomy') && ~isempty(sn_obj.anatomy)
    obj.anatomy = sn_obj.anatomy;
end

obj.stats = sn_obj.stats;
end


function reportObj = prepare_obj_for_langloc_report(obj, taskType, taskConfig, reportOutputDir)
% PREPARE_OBJ_FOR_LANGLOC_REPORT  Adapt object fields for generateReportLangloc*.
%
%   reportObj = prepare_obj_for_langloc_report(obj, taskType, taskConfig, reportOutputDir)
%
%   Returns a copy with:
%     - experiment name recognized by generateReportLangloc
%     - S/N condition labels mapped to report conventions
%     - a minimal events_table when behavioral logs are missing

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
end


function condReportMap = build_report_condition_map(taskType, taskConfig)
% Map pipeline condition names to labels expected by generateReportLangloc.

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
% Minimal behavioral table so the report script does not error.

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


function run_langloc_report(obj, taskType, taskConfig, reportOutputDir, subjectName, useV2)
% RUN_LANGLOC_REPORT  Generate the langloc PDF (non-fatal on missing toolbox).

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

if useV2
    generateReportLangloc_v2(reportObj, reportName);
else
    generateReportLangloc(reportObj, reportName);
end

fprintf('Langloc report saved: %s\n', fullfile(reportOutputDir, [reportName, '.pdf']));
end
