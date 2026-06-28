function gammaPower = build_gamma_power_from_sn_obj(sn_obj, taskConfig, opts)
% BUILD_GAMMA_POWER_FROM_SN_OBJ  All-channel bipolar gamma profiles for group analysis.
%
%   Output matches encoding_group.py expectations:
%     gammaPower.bip.labels, S/N/W/J_meanPower, S/N/W/J_word_ave (+ optional *_sem)

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'useOddForInference')
    opts.useOddForInference = false;
end

if isempty(sn_obj.bip_elec_data)
    error('sn_obj has no bipolar data.');
end

waBase = {'words', taskConfig.testWords, ...
    'allElecs', true, ...
    'use_odd_for_inference', opts.useOddForInference};

[S_ave, S_sem] = get_bipolar_word_averages(sn_obj, taskConfig.S_condition, waBase);
[N_ave, N_sem] = get_bipolar_word_averages(sn_obj, taskConfig.N_condition, waBase);

bip = struct();
bip.labels = sn_obj.bip_ch_label(:);
bip.S_word_ave = S_ave;
bip.N_word_ave = N_ave;
bip.S_word_ave_sem = S_sem;
bip.N_word_ave_sem = N_sem;
bip.S_meanPower = mean(S_ave, 2);
bip.N_meanPower = mean(N_ave, 2);

if isfield(taskConfig, 'W_condition') && ~isempty(taskConfig.W_condition)
    [W_ave, W_sem] = get_bipolar_word_averages(sn_obj, taskConfig.W_condition, waBase);
    bip.W_word_ave = W_ave;
    bip.W_word_ave_sem = W_sem;
    bip.W_meanPower = mean(W_ave, 2);
else
    bip.W_word_ave = nan(size(S_ave));
    bip.W_word_ave_sem = nan(size(S_sem));
    bip.W_meanPower = nan(size(bip.S_meanPower));
end

if isfield(taskConfig, 'J_condition') && ~isempty(taskConfig.J_condition)
    [J_ave, J_sem] = get_bipolar_word_averages(sn_obj, taskConfig.J_condition, waBase);
    bip.J_word_ave = J_ave;
    bip.J_word_ave_sem = J_sem;
    bip.J_meanPower = mean(J_ave, 2);
else
    bip.J_word_ave = nan(size(S_ave));
    bip.J_word_ave_sem = nan(size(S_sem));
    bip.J_meanPower = nan(size(bip.S_meanPower));
end

gammaPower = struct('bip', bip);
end


function [values_ave, values_ave_sem] = get_bipolar_word_averages(sn_obj, condition, waBase)
[va, va_sem, ~] = sn_obj.get_word_averages(waBase{:}, 'condition', condition);
bipIdx = find_bipolar_word_average_index(va, sn_obj);
values_ave = va{bipIdx};
values_ave_sem = va_sem{bipIdx};
end


function idx = find_bipolar_word_average_index(va, sn_obj)
nBip = numel(sn_obj.bip_ch_label);
for k = 1:numel(va)
    if size(va{k}, 1) == nBip
        idx = k;
        return;
    end
end
if numel(va) >= 2
    idx = 2;
else
    idx = 1;
end
end
