function taskConfig = get_mit_task_config(taskType)
% GET_MIT_TASK_CONFIG  Task metadata shared by MIT pipeline scripts.
%
%   taskConfig = get_mit_task_config('MITSWJNTask')

if nargin < 1 || isempty(taskType)
    error('taskType is required.');
end

switch taskType
    case 'MITLangloc'
        taskConfig = struct();
        taskConfig.nWordPositions = 12;
        taskConfig.wordDuration   = 0.45;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map(...
            {'sentence', 'word lists', 'Jabberwocky', 'non-word lists'}, ...
            {'Sentences', 'Word-lists', 'Jabberwocky', 'Nonword-lists'});
        taskConfig.S_condition = 'Sentences';
        taskConfig.N_condition = 'Nonword-lists';
        taskConfig.W_condition = 'Word-lists';
        taskConfig.J_condition = 'Jabberwocky';
        taskConfig.testWords   = 1:12;
        taskConfig.subAverage  = true;

    case 'MITSWJNTask'
        taskConfig = struct();
        taskConfig.nWordPositions = 8;
        taskConfig.wordDuration   = 0.70;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map(...
            {'sentence', 'word lists', 'Jabberwocky', 'non-word lists'}, ...
            {'SENTENCES', 'WORDS', 'JABBERWOCKY', 'NONWORDS'});
        taskConfig.S_condition = 'SENTENCES';
        taskConfig.N_condition = 'NONWORDS';
        taskConfig.W_condition = 'WORDS';
        taskConfig.J_condition = 'JABBERWOCKY';
        taskConfig.testWords   = 1:8;
        taskConfig.subAverage  = false;
        taskConfig.addFixationRow = false;
        taskConfig.fixationDuration = 0.5;

    case 'Auditory'
        taskConfig = struct();
        taskConfig.nWordPositions = 1;
        taskConfig.wordDuration   = 18;
        taskConfig.eventPattern   = '%s';
        taskConfig.conditionMap = containers.Map({'Intact', 'Degraded'}, {'Intact', 'Degraded'});
        taskConfig.S_condition = 'Intact';
        taskConfig.N_condition = 'Degraded';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1;
        taskConfig.subAverage  = true;

    case {'WM', 'Math', 'vWM'}
        taskConfig = struct();
        taskConfig.nWordPositions = 5;
        taskConfig.wordDuration   = 1;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map({'hard', 'easy'}, {'Hard', 'Easy'});
        taskConfig.S_condition = 'Hard';
        taskConfig.N_condition = 'Easy';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1:5;
        taskConfig.subAverage  = false;

    case {'MSIT', 'vMSIT'}
        taskConfig = struct();
        taskConfig.nWordPositions = 1;
        taskConfig.wordDuration   = 1.5;
        taskConfig.eventPattern   = '%s_tp%d';
        taskConfig.conditionMap = containers.Map({'hard', 'easy'}, {'Hard', 'Easy'});
        taskConfig.S_condition = 'Hard';
        taskConfig.N_condition = 'Easy';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1;
        taskConfig.subAverage  = false;

    case 'Naturalistic'
        taskConfig = struct();
        taskConfig.markers = {'NA', 'NA2', 'NA3'};
        taskConfig.segmentDurations = struct( ...
            'NA', 564.9824, ...
            'NA2', 642.6646, ...
            'NA3', 643.703);
        taskConfig.preBuffer  = 0;
        taskConfig.postBuffer = 0;
        taskConfig.segmentMode = 'duration';  % 'duration' | 'nextmarker'
        taskConfig.nWordPositions = 1;
        taskConfig.wordDuration   = 1;
        taskConfig.S_condition = '';
        taskConfig.N_condition = '';
        taskConfig.W_condition = '';
        taskConfig.J_condition = '';
        taskConfig.testWords   = 1;
        taskConfig.subAverage  = false;

    otherwise
        error('Unknown task type: %s', taskType);
end
end
