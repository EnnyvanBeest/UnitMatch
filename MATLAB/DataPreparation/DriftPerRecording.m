function DriftZ = DriftPerRecording(UMparam,MaxDriftUm)
% Depth (Z) offset of every recording of a UnitMatch run relative to its first recording, in um
% (add to a unit's depth). UMparam.drift is (recording pair k -> k+1, [X Y Z], pass):
%  pass 1: correction applied to recording k+1, estimated against recording k after k's own
%          correction, i.e. relative to the first recording (if the chain is not broken)
%  pass 2: what remains between recordings k and k+1 after pass 1 (measured, not applied)
% Offset of recording r: its pass-1 correction. Where pass 1 had too few matches for a pair (NaN:
% recording k+1 not corrected, and later pass-1 corrections are relative to it), the pass-2 estimate
% of that pair bridges the gap and is added to all later recordings. Pass-2 remainders of other
% pairs are not added (they are mostly noise and would accumulate over the chain).
% Pairs without any estimate count as 0 (a warning is given).
% Used by LoadPreparedClusInfo (clusinfo.AdjustedDepth) and ReplaceUniqueIDs.
if nargin < 2, MaxDriftUm = 300; end
nRec = numel(UMparam.KSDir);
DriftZ = zeros(nRec,1);
if ~isfield(UMparam,'drift') || isempty(UMparam.drift) || nRec < 2
    return
end
d = UMparam.drift;
nPair = min(nRec-1,size(d,1));
Applied = d(1:nPair,3,1);
if size(d,3) >= 2, Rest = d(1:nPair,3,end); else, Rest = zeros(nPair,1); end
Broken = isnan(Applied) & isnan(Rest);
if any(Broken)
    warning('LoadPreparedClusInfo:DriftChainBroken', ...
        'No drift estimate between recordings %s and the next (too few matches): later recordings may not be aligned to the first', ...
        mat2str(find(Broken)'))
end
Bridge = Rest .* isnan(Applied); % pass-2 estimate only for pairs that pass 1 could not align
Bridge(isnan(Bridge)) = 0; Applied(isnan(Applied)) = 0;
DriftZ(2:nPair+1) = Applied + cumsum(Bridge);
TooLarge = abs(DriftZ) > MaxDriftUm;
if any(TooLarge)
    warning('LoadPreparedClusInfo:DriftTooLarge', ...
        'Drift estimate of recording(s) %s larger than %d um (%s um): not used, depth left uncorrected', ...
        mat2str(find(TooLarge)'),MaxDriftUm,mat2str(round(DriftZ(TooLarge)')))
    DriftZ(TooLarge) = 0;
end
end
