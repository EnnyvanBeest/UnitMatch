function [Bad, Stats] = BadLFPChannels(Y, Shank, Depth, Opt)
% [Bad, Stats] = BadLFPChannels(Y, Shank, Depth, Opt)
% Bad LFP channels (noisy, dead or with a different gain), to leave out of CSDs and depth features:
% one deviating site becomes a narrow sink-source-sink stripe in every event-triggered CSD.
%  Y: time x channel LFP, or a cell array of such segments (e.g. the 1-40 Hz LFP of LFPDepthFeatures)
%  Shank, Depth (um) per channel
% Per channel, averaged over the segments (common median reference removed): log10 variance, and the
% correlation with the median of the channels within Opt.NbUm on the same shank. Both are compared
% with the channels above and with those below separately (the other column at the same depth counts
% for both), and the deviation from the closer side is used, so a real change along the probe (brain
% surface, area boundary) is not flagged. Bad: log variance more than Opt.SD robust SDs (over the
% shank; floor Opt.MinPowSD log10 units) away from both sides, or correlation more than Opt.SD robust
% SDs (floor Opt.MinCorrSD) below both sides.
% Bad: nChan x 1 logical; Stats: LogVar, LocCorr, zPow, zCorr per channel
if nargin < 4, Opt = struct; end
Def = struct('NbUm',100,'SD',5,'MinPowSD',0.05,'MinCorrSD',0.01);
for fn = fieldnames(Def)', if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end, end
if ~iscell(Y), Y = {Y}; end
Shank = Shank(:); Depth = Depth(:); nRec = numel(Shank);
% neighbours above (incl. same depth) and below (incl. same depth), within Opt.NbUm on the same shank
dd = Depth' - Depth; same = Shank == Shank';
nbUp = arrayfun(@(c) setdiff(find(same(c,:) & dd(c,:) > -5 & dd(c,:) <= Opt.NbUm),c),(1:nRec)','Uni',0);
nbDn = arrayfun(@(c) setdiff(find(same(c,:) & dd(c,:) < 5 & dd(c,:) >= -Opt.NbUm),c),(1:nRec)','Uni',0);
LV = 0; LC = 0; nS = 0;
for s = 1:numel(Y)
    X = double(Y{s}(:,1:nRec));
    if size(X,1) < 10, continue; end
    X = X - median(X,2);
    lc = nan(nRec,2);
    for c = 1:nRec
        if ~isempty(nbUp{c}), lc(c,1) = corr(X(:,c),median(X(:,nbUp{c}),2)); end
        if ~isempty(nbDn{c}), lc(c,2) = corr(X(:,c),median(X(:,nbDn{c}),2)); end
    end
    LC = LC + max(lc,[],2); % the better-matching side
    LV = LV + log10(max(var(X,0,1),realmin))';
    nS = nS + 1;
end
Bad = false(nRec,1);
Stats = struct('LogVar',nan(nRec,1),'LocCorr',nan(nRec,1),'zPow',nan(nRec,1),'zCorr',nan(nRec,1));
if nS == 0, return; end
LV = LV/nS; LC = LC/nS;
% deviation from the closer side
rP = nan(nRec,1); rC = nan(nRec,1);
for c = 1:nRec
    pv = [median(LV(nbUp{c})) median(LV(nbDn{c}))]; pv = pv(~isnan(pv));
    if ~isempty(pv), [~,k] = min(abs(LV(c) - pv)); rP(c) = LV(c) - pv(k); end
    rC(c) = LC(c) - min([median(LC(nbUp{c}),'omitnan') median(LC(nbDn{c}),'omitnan') Inf]); % below both sides
end
rC(isinf(rC)) = NaN;
zP = nan(nRec,1); zC = nan(nRec,1);
for sh = unique(Shank)'
    in = Shank == sh;
    okP = in & isfinite(rP); okC = in & isfinite(rC);
    zP(in) = (rP(in) - median(rP(okP)))/max(1.4826*mad(rP(okP),1),Opt.MinPowSD);
    zC(in) = (rC(in) - median(rC(okC)))/max(1.4826*mad(rC(okC),1),Opt.MinCorrSD);
end
Bad = abs(zP) > Opt.SD | zC < -Opt.SD;
Stats = struct('LogVar',LV,'LocCorr',LC,'zPow',zP,'zCorr',zC);
end
