%=========================================================================%
% ecog_:plus  ––  Concatenate two recordings end-to-end
%=========================================================================%
function obj3 = plus(obj1, obj2)
    % Validation
    assert(isa(obj1,'ecog_data_ieeg') && isa(obj2,'ecog_data_ieeg'), ...
        'ecog_plus:InvalidInput','Both must be ecog_data_ieeg objects');
    assert(strcmp(obj1.subject,obj2.subject), ...
        'ecog_plus:SubjectMismatch','Subject IDs differ');
    assert(abs(obj1.sample_freq-obj2.sample_freq)<1e-6, ...
        'ecog_plus:SamplingFrequencyMismatch','Fs mismatch');
    assert(size(obj1.elec_data,1)==size(obj2.elec_data,1), ...
        'ecog_plus:UnipolarChannelMismatch','Unipolar channels differ');
    if isprop(obj1,'bip_elec_data') || isprop(obj2,'bip_elec_data')
        assert(isprop(obj1,'bip_elec_data') && isprop(obj2,'bip_elec_data'), ...
            'ecog_plus:BipolarPresenceMismatch','Bipolar presence mismatch');
        assert(size(obj1.bip_elec_data,1)==size(obj2.bip_elec_data,1), ...
            'ecog_plus:BipolarChannelMismatch','Bipolar channels differ');
    end

    % Offsets
    n1      = size(obj1.elec_data,2);
    offSamp = n1;
    offSec  = offSamp / obj1.sample_freq;

    % % Copy every *public* property
    % for p = string(properties(obj1))'
    %     fn = p{1};           % char name
    %     obj3.(fn) = obj1.(fn);
    % end
    % 
    % % Clone obj1
    % ctorArgs = {obj1.for_preproc};      % REQUIRED constructor arg #1
    % % Append here any other mandatory args, e.g. ctorArgs{2}=obj1.cfg;
    % 
    % obj3 = feval(class(obj1), ctorArgs{:});   % brand-new instance
    % 
    objtemp = obj1;
    obj3 = objtemp;

    % Concatenate signals
    obj3.elec_data = [obj1.elec_data, obj2.elec_data];
    obj3.for_preproc.elec_data = [obj1.for_preproc.elec_data, obj2.for_preproc.elec_data];
    if isprop(obj1,'bip_elec_data')
        obj3.bip_elec_data = [obj1.bip_elec_data, obj2.bip_elec_data];
        obj3.for_preproc.bip_elec_data = [obj1.for_preproc.bip_elec_data, obj2.for_preproc.bip_elec_data];
    end

    % trial_timing (samples only)
    obj3.trial_timing = [ toCell(obj1.trial_timing);
                          cellfun(@(T) shiftSamples(T,offSamp), ...
                                  toCell(obj2.trial_timing), ...
                                  'UniformOutput',false) ];
    obj3.for_preproc.trial_timing = [ toCell(obj1.for_preproc.trial_timing);
                          cellfun(@(T) shiftSamples(T,offSamp), ...
                                  toCell(obj2.for_preproc.trial_timing), ...
                                  'UniformOutput',false) ];

    % events_table (mixed units with updated rules)
    obj3.events_table = mergeEvents(obj1.events_table, obj2.events_table, offSamp, offSec);
    obj3.for_preproc.events_table = mergeEvents(obj1.for_preproc.events_table, obj2.for_preproc.events_table, offSamp, offSec);
    % Merge condition & session
    if isprop(obj1,'condition')
        obj3.condition = [obj1.condition; obj2.condition];
        obj3.for_preproc.condition = [obj1.for_preproc.condition; obj2.for_preproc.condition];
    end
    if isprop(obj1,'session')
        obj3.session   = [obj1.session; obj2.session + max(obj1.session)];
        obj3.for_preproc.session   = [obj1.for_preproc.session; obj2.for_preproc.session + max(obj1.for_preproc.session)];
    end

    % Stitch index & history
    if isprop(obj1,'stitch_index')
        obj3.stitch_index = [obj1.stitch_index, offSamp];
        obj3.for_preproc.stitch_index = [obj1.for_preproc.stitch_index, offSamp];
    else
        obj3.stitch_index = offSamp;
        obj3.for_preproc.stitch_index = offSamp;
    end

    % Duration
    if isprop(obj1,'recording_duration')
        obj3.recording_duration = size(obj3.elec_data,2)/obj3.sample_freq;
        obj3.for_preproc.recording_duration = size(obj3.for_preproc.elec_data,2)/obj3.for_preproc.sample_freq;
    end

    % Log
    if ~isprop(obj3,'processing_history'), obj3.processing_history={}; end
    obj3.processing_history{end+1} = struct(...
        'operation','concatenate','timestamp',datestr(now),...
        'offset_samples',offSamp,'offset_seconds',offSec);
end

%======================== Helper Functions ===============================%

function C = toCell(x)
    if iscell(x), C = x(:); else C = {x}; end
end

function T = shiftSamples(T,off)
    if istable(T)
        numCols = varfun(@isnumeric,T,'OutputFormat','uniform');
        T{:,numCols} = T{:,numCols} + off;
    end
end

function ET = mergeEvents(E1,E2,offS,offT)
    if isempty(E1), ET = shiftEvents(E2,offS,offT); return; end
    if isempty(E2), ET = E1;               return; end
    ET = [E1; shiftEvents(E2,offS,offT)];
end

function T = shiftEvents(T,offS,offT)
    % Define column unit categories
    sampleCols = ["time_info"];              % leave 'list' unchanged
    secondCols = [ ...
        "*_natus"                             % any column ending with _natus
        "trial_onset"                         % seconds
        "audio_ended" 
        "RT"                                  % keep unchanged? user said keep RT same => we will NOT shift RT
    ];

    % Build dynamic list of nat­us columns
    allVars = string(T.Properties.VariableNames);
    natusCols = allVars(endsWith(allVars,"_natus"));
    secondCols = unique([setdiff(secondCols,"RT"); natusCols']);

    for v = allVars
        if ~isnumeric(T.(v)), continue; end
        if strcmp(v,"trial") || strcmp(v,"RT")
            % keep these unchanged
        elseif ismember(v, sampleCols)
            T.(v) = T.(v) + offS;
        elseif any(strcmp(v,secondCols))
            T.(v) = T.(v) + offT;
        end
    end
end
