function pairs = match_lang_channels_by_mni(profilesA, profilesB, maxDistMm)
% MATCH_LANG_CHANNELS_BY_MNI  Pair lang channels by nearest MNI (within radius).
%
%   pairs = match_lang_channels_by_mni(profilesA, profilesB, maxDistMm)
%
%   Greedy mutual-nearest matching among channels with valid MNI.
%   profilesA / profilesB are struct arrays from build_lang_channel_profiles.

if nargin < 3 || isempty(maxDistMm)
    maxDistMm = 12;
end

pairs = struct([]);
if isempty(profilesA) || isempty(profilesB)
    return;
end

nA = numel(profilesA);
nB = numel(profilesB);
D = inf(nA, nB);
for i = 1:nA
    for j = 1:nB
        D(i, j) = norm(profilesA(i).mni - profilesB(j).mni);
    end
end

usedA = false(nA, 1);
usedB = false(nB, 1);
pairList = {};

while true
    Dmask = D;
    Dmask(usedA, :) = inf;
    Dmask(:, usedB) = inf;
    [dmin, idx] = min(Dmask(:));
    if ~isfinite(dmin) || dmin > maxDistMm
        break;
    end
    [ia, ib] = ind2sub(size(D), idx);

    pairList{end+1} = struct( ... %#ok<AGROW>
        'idxA', ia, 'idxB', ib, ...
        'distMm', dmin, ...
        'profileA', profilesA(ia), ...
        'profileB', profilesB(ib));
    usedA(ia) = true;
    usedB(ib) = true;
end

if ~isempty(pairList)
    pairs = [pairList{:}]';
end

end
