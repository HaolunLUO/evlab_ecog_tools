function diag_mni()
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

files = {
    'F:\seeg\luohong\analysisEV\Subject01_MITSWJNTask_crunched.mat'
    'F:\seeg\luohong\analysisEV\crunched\MITLangloc\AMC092_MITLangloc_crunched_HG_ZScore.mat'
};

for f = 1:numel(files)
    matFile = files{f};
    fprintf('\n=== %s ===\n', matFile);
    if ~isfile(matFile)
        fprintf('MISSING\n');
        continue;
    end
    w = whos('-file', matFile);
    fprintf('vars: %s\n', strjoin({w.name}, ', '));
    S = load(matFile, 'sn_obj', 'obj');
    sn = S.sn_obj;
    fprintf('sn class: %s\n', class(sn));
    inspect_obj(sn, 'sn_obj');
    if isfield(S, 'obj') && ~isempty(S.obj)
        fprintf('obj class: %s\n', class(S.obj));
        inspect_obj(S.obj, 'obj');
    end
    if isprop(sn, 'anatomy') && ~isempty(sn.anatomy)
        if isfield(sn.anatomy, 'mapping')
            fprintf('mapping len: %d\n', numel(sn.anatomy.mapping));
        end
        if isfield(sn.anatomy, 'mni_space')
            ms = sn.anatomy.mni_space;
            if isfield(ms, 'tala')
                fprintf('mni tala electrodes: %s\n', mat2str(size(ms.tala.electrodes)));
            elseif isfield(ms, 'vera_mat_minimal')
                fprintf('mni vera electrodes: %s\n', mat2str(size(ms.vera_mat_minimal.tala.electrodes)));
            end
        end
    end
    sig = find(sn.s_vs_n_sig.bip_elec_data{1, 1});
    fprintf('sig bipolar: %d\n', numel(sig));
    if ~isempty(sig)
        [m, ok] = get_lang_channel_mni_coords(sn, sig(1), 'bipolar');
        fprintf('first sig (%s) mni ok=%d: %s\n', sn.bip_ch_label{sig(1)}, ok, mat2str(m, 3));
    end
end
end

function inspect_obj(o, tag)
props = {'elec_ch_pos_mni', 'bip_ch_pos_mni', 'anatomy', 'filt_ops'};
for i = 1:numel(props)
    p = props{i};
    if isprop(o, p)
        v = o.(p);
        if iscell(v)
            fprintf('%s.%s: cell %d\n', tag, p, numel(v));
        elseif isnumeric(v)
            fprintf('%s.%s: numeric %s\n', tag, p, mat2str(size(v)));
        elseif isstruct(v)
            fprintf('%s.%s: struct fields %s\n', tag, p, strjoin(fieldnames(v), ','));
        else
            fprintf('%s.%s: %s nonempty=%d\n', tag, p, class(v), ~isempty(v));
        end
    else
        fprintf('%s.%s: MISSING\n', tag, p);
    end
end
end
