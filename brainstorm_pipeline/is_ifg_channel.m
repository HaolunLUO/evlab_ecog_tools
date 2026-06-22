function [tf, info] = is_ifg_channel(sn_obj, chanIdx, signalType)
% IS_IFG_CHANNEL  True when a channel lies in inferior frontal gyrus (IFG).
%
%   [tf, info] = is_ifg_channel(sn_obj, chanIdx, signalType)
%
%   Detection order:
%     1) HCP-MMP / anatomy text labels (any contact for bipolar pairs)
%     2) MNI coordinate bounding box (hemisphere-specific IFG envelope)
%
%   info.method is 'hcp', 'mni', or '' when not IFG.
%   info.hcp_labels lists matched anatomy strings.

if nargin < 3 || isempty(signalType)
    signalType = 'bipolar';
end

info = struct('method', '', 'hcp_labels', {{}}, 'mni', [NaN NaN NaN]);

if isempty(chanIdx) || chanIdx < 1
    tf = false;
    return;
end

hcpLabels = get_channel_hcp_labels(sn_obj, chanIdx, signalType);
info.hcp_labels = hcpLabels;
for i = 1:numel(hcpLabels)
    if is_ifg_hcp_label(hcpLabels{i})
        tf = true;
        info.method = 'hcp';
        [info.mni, ~] = get_lang_channel_mni_coords(sn_obj, chanIdx, signalType);
        return;
    end
end

[mni, ok] = get_lang_channel_mni_coords(sn_obj, chanIdx, signalType);
info.mni = mni;
if ok && is_mni_in_ifg(mni)
    tf = true;
    info.method = 'mni';
    return;
end

tf = false;
end


function tf = is_mni_in_ifg(mni)
% Generous MNI152 IFG envelope (pars opercularis / triangularis / orbitalis).
x = mni(1); y = mni(2); z = mni(3);
if x < 0
    tf = x >= -68 && x <= -28 && y >= -8 && y <= 48 && z >= -18 && z <= 38;
else
    tf = x >= 28 && x <= 68 && y >= -8 && y <= 48 && z >= -18 && z <= 38;
end
end
