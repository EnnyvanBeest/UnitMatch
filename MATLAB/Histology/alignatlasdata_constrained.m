function [Depth2Area, AlignmentFeatures, Fit] = alignatlasdata_constrained(histinfo, trackcoordinates, sp, clusinfo, AllenCCFPath, Opt)
% [Depth2Area, AlignmentFeatures, Fit] = alignatlasdata_constrained(histinfo, trackcoordinates, sp, clusinfo, AllenCCFPath, Opt)
% Aligns the histology track(s) of a probe to its recording depths with a constrained model:
%  1. Ephys profile per shank along the probe (depth = site position, um): log spike rate (all
%     spikes), density of good units and mean amplitude, smoothed. Activity = mean of the z-scored
%     rate and unit density (high in grey matter; low in white matter, ventricles and outside the brain).
%  2. Predicted profile from the histology: grey matter 1 (+ a little for the atlas' expected firing
%     rate), fibre tracts, ventricles and outside the brain 0.
%  3. Global fit, shared by all shanks (inserted together): depth = Scale * (distance of the track
%     point above the track end) + Offset, by maximum correlation between the two profiles. The scale
%     captures shrinkage/registration (Opt.ScaleRange), the offset where the tip is relative to the
%     track end. Each shank's offset is then refined within Opt.ShankOffsetUm.
%  4. Local warp per shank: a smooth displacement of the histology along the probe (knots every
%     Opt.KnotUm, at most Opt.MaxWarpUm, local stretch within 1 +- Opt.MaxStretchDev), fitted by
%     dynamic programming, with a pull towards the global fit (Opt.WarpPenalty).
% Inputs as alignatlasdata_automated: histinfo {shank} tables from the BrainGlobe track csv files
% (distance from the first (surface) point, region acronym), trackcoordinates {shank} (npy, atlas
% coordinates per track point), sp/clusinfo (LoadPreparedClusInfo; clusinfo.Shank 0-based).
% Output:
%  Depth2Area: one row per histology track point in the recorded range: Depth (aligned, um),
%    Shank (position among the shanks used, 1-based, as before), ShankAbs (shank number, 0-based),
%    Area, Color, Coordinates (track point), BaselineDepth (global fit only), WarpUm (Depth - BaselineDepth)
%    TrackPosUm (position of the track point above the track end, atlas um)
%  AlignmentFeatures: {shank} struct with the ephys profile at its own probe depths: Depth, Values
%    (nDepth x nFeatures), Names, Activity, the predicted profile with the global fit / aligned, and
%    the geometry: TipDepth, SurfaceTrackPos, TrackSpacingUm, Scale, Offset; LandmarkNames,
%    Landmarks (z-scored profiles), LandmarkAreasGlobal/Aligned (area indicators at the probe depths)
%  Fit: Scale, Offset, per-shank offsets, correlations of the global fit and after the warp;
%    Alternatives (Scale, Offset, Corr, Score, Components: the correlation per part of the score,
%    named in Fit.ComponentNames: activity and each landmark); LandmarkWarnings

% Prior on the offset: the deepest histology point is taken as the probe tip, which lies
% Opt.TipToSiteUm below the first recording site (depth 0). The track end can be far off (traced past
% the tip, tissue damage; up to ~1 mm seen), so the prior is weak and capped: penalty
% Opt.PriorWeight * min((offset error/Opt.OffsetPriorSdUm)^2, Opt.PriorCapSD^2) in correlation units.
% Landmarks: profiles that should be high in specific areas, each fitted next to the activity profile
% with its weight (score = (r_activity + sum w_k r_k)/(1 + sum w_k)). The activity profile dominates:
% the landmark weights together are at most Opt.MaxLandmarkWeight (relative to activity 1).
%  - LFP (optional, Opt.LFP from LFPDepthFeatures): the ripple event rate -> Ammon's horn (CA1-3),
%    weight Opt.LFPWeight. Only used if its peak is a narrow band (<= Opt.RippleMaxWidthUm at half
%    height): a broad band is not the CA1 pyramidal layer (e.g. thalamic bursts share the band)
%  - Opt.Landmarks (optional, e.g. TuningLandmarks): per unit Shank, Depth and Value (0/1), with
%    AreaPaths (Allen structure id path tokens) and Weight; the fraction of positive units per depth
%    (smoothed over Opt.LandmarkSmoothUm) is fitted to those areas, e.g. units with a receptive field
%    to visual areas. A landmark whose areas are not in the recorded range gives no evidence.
% Brain surface from the data: units cannot lie outside the brain, so the aligned surface may not be
% below the highest good units (hard-ish penalty). If the top >= Opt.SurfaceGapUm of the probe has no
% good units, the probe was partly outside the brain: the surface is expected Opt.SurfaceAboveUnitsUm
% above the highest units (+- Opt.SurfaceSdUm; weight Opt.SurfaceWeight per SD^2, capped at 3 SD).
if nargin < 6, Opt = struct; end
Def = struct('BinUm',10,'SmoothUm',30,'ScaleRange',[0.8 1.3],'ScaleStep',0.01,'OffsetRangeUm',[-1000 1000], ...
    'OffsetStepUm',10,'ShankOffsetUm',100,'KnotUm',50,'MaxWarpUm',150,'MaxStretchDev',0.3,'WarpPenalty',0.5,'MarginUm',50, ...
    'TipToSiteUm',175,'OffsetPriorSdUm',300,'PriorWeight',0.05,'PriorCapSD',2,'LFP',[],'LFPWeight',0.3, ...
    'MaxLandmarkWeight',0.5,'RippleMaxWidthUm',700, ...
    'SurfaceGapUm',300,'SurfaceAboveUnitsUm',100,'SurfaceSdUm',150,'SurfaceWeight',0.05, ...
    'Landmarks',struct('Name',{},'Shank',{},'Depth',{},'Value',{},'AreaPaths',{},'Weight',{}),'LandmarkSmoothUm',100);
Def.OffsetRangeUm = [-2000 1000]; % the track end can be >1 mm deeper than the tip
for fn = fieldnames(Def)'
    if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end
end

%% Atlas information per area
st = readtable(fullfile(AllenCCFPath,'structure_tree_safe_2017.csv'));
stAcr = lower(strrep(st.acronym,'/',''));
stCol = cellfun(@(c) pad(char(string(c)),6,'left','0'),st.color_hex_triplet,'Uni',0);
if ismember('expected_firing_rate',st.Properties.VariableNames), stEFR = st.expected_firing_rate; else, stEFR = ones(height(st),1); end

if ~iscell(histinfo), histinfo = {histinfo}; end
if ~iscell(trackcoordinates), trackcoordinates = {trackcoordinates}; end
ShanksUsed = unique(clusinfo.Shank(:))';
histinfo = histinfo(ShanksUsed+1);
trackcoordinates = trackcoordinates(ShanksUsed+1);
nSh = numel(ShanksUsed);

%% 1. Ephys profiles
% Shank of every spike via its cluster (cluster ids repeat across recordings: use recording too)
[isClu,cluRow] = ismember([double(sp.clu(:)) double(sp.RecSes(:))],[double(clusinfo.cluster_id(:)) double(clusinfo.RecSesID(:))],'rows');
spkShank = nan(numel(sp.clu),1); spkShank(isClu) = clusinfo.Shank(cluRow(isClu));
Dur = 0;
for r = unique(sp.RecSes(:))', t = sp.st(sp.RecSes==r); Dur = Dur + (max(t)-min(t)); end
Edges = floor(min(sp.spikeDepths)/Opt.BinUm)*Opt.BinUm : Opt.BinUm : ceil(max(sp.spikeDepths)/Opt.BinUm)*Opt.BinUm + Opt.BinUm;
dGrid = (Edges(1:end-1) + Edges(2:end))'/2;
Sm = @(x) smoothdata(x,1,'gaussian',max(3,round(5*Opt.SmoothUm/Opt.BinUm)));
Z = @(x) (x - mean(x,'omitnan'))./max(std(x,'omitnan'),eps);
E = cell(1,nSh); Feat = cell(1,nSh);
% Good units: noise clusters (e.g. on channels outside the brain) would make those channels look active
spkGood = false(numel(sp.clu),1); spkGood(isClu) = logical(clusinfo.Good_ID(cluRow(isClu)));
for si = 1:nSh
    sel = spkShank == ShanksUsed(si) & spkGood;
    rate = histcounts(sp.spikeDepths(sel),Edges)'/max(Dur,1);               % spikes/s per bin (good units)
    good = logical(clusinfo.Good_ID(:)) & clusinfo.Shank(:) == ShanksUsed(si);
    dens = histcounts(clusinfo.depth(good),Edges)';                          % good units per bin
    b = discretize(sp.spikeDepths(sel),Edges); amps = double(sp.spikeAmps(sel));
    amp = accumarray(b(~isnan(b)),amps(~isnan(b)),[numel(dGrid) 1],@mean,NaN);
    LogRate = Sm(log10(rate + 0.01)); Dens = Sm(dens); Amp = Sm(fillmissing(amp,'linear','EndValues','nearest'));
    E{si} = (Z(LogRate) + Z(Dens))/2;
    Feat{si} = struct('Depth',dGrid,'Values',[LogRate Dens Amp],'Names',{{'log10 spike rate (good units)','good units per bin','mean amplitude'}}, ...
        'Activity',E{si},'Ripple',[]);
end
% Landmarks: profiles along the probe that should be high in specific areas (AreaPaths, Allen
% structure id path tokens), each fitted with its weight next to the activity profile.
% LM{si}: nDepth x nLandmark (z-scored, NaN where unknown)
Fit = struct; Fit.LandmarkWarnings = {}; Fit.LandmarkNotes = {};
LMName = {}; LMAreas = {}; LMW = [];
LM = repmat({zeros(numel(dGrid),0)},1,nSh);
% LFP features on the same depth bins (per shank); the ripple event rate marks Ammon's horn
UseLFP = ~isempty(Opt.LFP);
if UseLFP
    ripCol = find(contains(Opt.LFP.Names,'ripple events'),1); % event rate: not dominated by spikes, unlike ripple-band power
    if isempty(ripCol), ripCol = find(contains(Opt.LFP.Names,'ripple'),1); end
    LMName{end+1} = 'LFP ripple events'; LMAreas{end+1} = {'/375/'}; LMW(end+1) = Opt.LFPWeight;
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        ch = Opt.LFP.Shank == ShanksUsed(si);
        if ~any(ch), continue; end
        b = discretize(Opt.LFP.Depth(ch),Edges); V = Opt.LFP.Values(ch,:); okb = ~isnan(b);
        Vb = nan(numel(dGrid),size(V,2));
        for k = 1:size(V,2), Vb(:,k) = accumarray(b(okb),V(okb,k),[numel(dGrid) 1],@mean,NaN); end
        Vb = Sm(fillmissing(Vb,'linear','EndValues','nearest'));
        Feat{si}.Values = [Feat{si}.Values Vb];
        Feat{si}.Names = [Feat{si}.Names strcat('LFP: ',Opt.LFP.Names)];
        if ~isempty(ripCol)
            if contains(Opt.LFP.Names{ripCol},'events'), r = Z(log10(Vb(:,ripCol) + 0.1)); else, r = Z(Vb(:,ripCol)); end
            % Width at half height of the band around the peak (above the median)
            [pk,ipk] = max(r); hm = median(r) + (pk - median(r))/2;
            lo = find(r(1:ipk) < hm,1,'last'); if isempty(lo), lo = 0; end
            hi = find(r(ipk:end) < hm,1,'first'); if isempty(hi), hi = numel(r) - ipk + 2; end
            widthUm = ((ipk + hi - 1) - lo - 1)*Opt.BinUm;
            if widthUm > Opt.RippleMaxWidthUm
                Fit.LandmarkNotes{end+1} = sprintf('Shank %d: ripple events in a %.0f um wide band (peak %.0f um): not used (not specific to CA1)', ...
                    ShanksUsed(si),widthUm,dGrid(ipk));
                fprintf('%s\n',Fit.LandmarkNotes{end})
            else
                LM{si}(:,end) = r;
            end
        end
    end
end
% Unit landmarks (e.g. TuningLandmarks: fraction of units with a receptive field / locked to a
% stimulus): fraction per depth, smoothed over Opt.LandmarkSmoothUm; unknown where < 2 units nearby
SmL = @(x) smoothdata(x,1,'gaussian',max(3,round(5*Opt.LandmarkSmoothUm/Opt.BinUm)));
MinDens = 2/(2.5*Opt.LandmarkSmoothUm/Opt.BinUm); % ~2 units within +- 1 sigma
for k = 1:numel(Opt.Landmarks)
    L = Opt.Landmarks(k);
    LMName{end+1} = L.Name; LMAreas{end+1} = cellstr(L.AreaPaths); LMW(end+1) = L.Weight; %#ok<AGROW>
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        on = L.Shank(:) == ShanksUsed(si);
        if sum(on) < 10 || ~any(L.Value(on)), continue; end
        cnt = SmL(histcounts(L.Depth(on),Edges)');
        pos = SmL(histcounts(L.Depth(on & L.Value(:) > 0),Edges)');
        frac = pos./cnt; frac(cnt < MinDens) = NaN;
        if mean(~isnan(frac)) < 0.2, continue; end
        LM{si}(:,end) = Z(frac);
        Feat{si}.Values = [Feat{si}.Values frac];
        Feat{si}.Names = [Feat{si}.Names {['tuning: ' L.Name]}];
    end
end
nLM = numel(LMName);
for si = 1:nSh
    Feat{si}.LandmarkNames = LMName; Feat{si}.Landmarks = LM{si};
    rc = find(strcmp(LMName,'LFP ripple events'),1);
    if ~isempty(rc), Feat{si}.Ripple = LM{si}(:,rc); end
end

%% 2. Predicted profile per histology track point
Hist = cell(1,nSh);
for si = 1:nSh
    H = histinfo{si}; vn = H.Properties.VariableNames;
    distV = vn{find(contains(vn,'Distance','IgnoreCase',true),1)};
    acrV  = vn{find(contains(vn,'Acronym','IgnoreCase',true),1)};
    tTrack = double(H.(distV)); acr = H.(acrV);
    nPt = min(numel(tTrack),size(trackcoordinates{si},1));
    tTrack = tTrack(1:nPt); acr = acr(1:nPt); coords = trackcoordinates{si}(1:nPt,:);
    [~,row] = ismember(lower(strrep(acr,'/','')),stAcr);
    pred = zeros(nPt,1); lm = zeros(nPt,nLM); col = repmat({'000000'},nPt,1);
    for k = 1:nPt
        if row(k) == 0, continue; end % 'Not found in brain' etc.
        col{k} = stCol{row(k)};
        pth = st.structure_id_path{row(k)};
        for j = 1:nLM, lm(k,j) = any(cellfun(@(a) contains(pth,a),LMAreas{j})); end % e.g. Ammon's horn for ripples
        if contains(pth,'/1009/') || contains(pth,'/73/') || strcmpi(stAcr{row(k)},'root'), continue; end % fibre tract, ventricle
        efr = stEFR(row(k)); if isnan(efr) || efr <= 0, efr = 1; end
        pred(k) = 1 + 0.25*log10(efr);
    end
    sTrack = max(tTrack) - tTrack; % distance above the track end (um, atlas)
    [sTrack,o] = sort(sTrack,'ascend'); % from the track end upwards
    Hist{si} = struct('s',sTrack,'pred',pred(o),'lm',lm(o,:),'acr',{acr(o)},'col',{col(o)},'coords',coords(o,:));
    Hist{si}.sSurface = max(sTrack(~strcmpi(acr(o),'Not found in brain')));
end
% A landmark that is clearly present while its areas are nowhere on the track: the track (or the
% landmark's area list) may be wrong; it cannot help the fit then
for k = 1:numel(Opt.Landmarks)
    j = find(strcmp(LMName,Opt.Landmarks(k).Name),1);
    for si = 1:nSh
        on = Opt.Landmarks(k).Shank(:) == ShanksUsed(si);
        if sum(on) >= 10 && mean(Opt.Landmarks(k).Value(on) > 0) >= 0.2 && ~any(Hist{si}.lm(:,j))
            Fit.LandmarkWarnings{end+1} = sprintf('Shank %d: %.0f%% of units "%s", but none of its areas on the histology track', ...
                ShanksUsed(si),100*mean(Opt.Landmarks(k).Value(on) > 0),LMName{j});
            warning('alignatlasdata_constrained:Landmark','%s',Fit.LandmarkWarnings{end})
        end
    end
end
% Highest good units per shank (and whether the top of the probe was without units: outside the brain)
TopUnit = nan(1,nSh); SurfaceObs = nan(1,nSh);
for si = 1:nSh
    gd = clusinfo.depth(logical(clusinfo.Good_ID(:)) & clusinfo.Shank(:) == ShanksUsed(si));
    if numel(gd) < 5, continue; end
    TopUnit(si) = prctile(gd,99);
    if dGrid(end) - TopUnit(si) >= Opt.SurfaceGapUm
        SurfaceObs(si) = TopUnit(si) + Opt.SurfaceAboveUnitsUm;
    end
end
SurfacePen = @(a,b) local_surfpen(a,b,Hist,TopUnit,SurfaceObs,Opt);
% Predicted value at probe depths d for a mapping s(d); NaN below the track end, 0 above the surface
PredAt = @(H,s) local_pred(H,s,'pred');
LmAt   = @(H,s) local_pred(H,s,'lm');
% Landmark weights per shank (0 for landmarks without data on that shank)
Wsh = cellfun(@(L) LMW .* any(~isnan(L),1),LM,'Uni',0);
Wsh = cellfun(@(w) w*min(1,Opt.MaxLandmarkWeight/max(sum(w),eps)),Wsh,'Uni',0); % activity dominates
% Score of a mapping s(d) for one shank: correlation of the activity with the predicted profile,
% combined with the correlation of each landmark with its areas
ScoreAt = @(si,s) local_score(E{si},LM{si},PredAt(Hist{si},s),LmAt(Hist{si},s),Wsh{si});
OffsetPrior = @(b) -Opt.PriorWeight*min(((b + Opt.TipToSiteUm)/Opt.OffsetPriorSdUm).^2, Opt.PriorCapSD^2);

%% 3. Global fit (scale and offset shared by the shanks), then per-shank offset
Scales = Opt.ScaleRange(1):Opt.ScaleStep:Opt.ScaleRange(2);
Offsets = Opt.OffsetRangeUm(1):Opt.OffsetStepUm:Opt.OffsetRangeUm(2);
C = nan(numel(Scales),numel(Offsets));
for ia = 1:numel(Scales)
    for ib = 1:numel(Offsets)
        r = arrayfun(@(si) ScoreAt(si,(dGrid - Offsets(ib))/Scales(ia)),1:nSh);
        C(ia,ib) = mean(r,'omitnan');
    end
end
Cp = C + OffsetPrior(Offsets); % with the prior on the offset
for ia = 1:numel(Scales) % and the brain surface constraint
    Cp(ia,:) = Cp(ia,:) - arrayfun(@(b) SurfacePen(Scales(ia),b),Offsets);
end
[~,best] = max(Cp(:)); [ia,ib] = ind2sub(size(Cp),best);
Fit.Scale = Scales(ia); Fit.Offset = Offsets(ib); Fit.CorrGlobal = C(ia,ib);
Fit.Landscape = struct('Scales',Scales,'Offsets',Offsets,'Corr',C,'CorrWithPrior',Cp);
% Alternative solutions (for the curator): best per offset over scales, local maxima >= 200 um apart
[cOff,iaOff] = max(Cp,[],1);
[~,loc] = findpeaks(fillmissing(cOff,'constant',-10),'MinPeakDistance',round(200/Opt.OffsetStepUm),'SortStr','descend');
loc = loc(1:min(4,end));
Fit.Alternatives = struct('Scale',num2cell(Scales(iaOff(loc))),'Offset',num2cell(Offsets(loc)), ...
    'Corr',num2cell(C(sub2ind(size(C),iaOff(loc),loc))),'Score',num2cell(cOff(loc)));
% Parts of the score per alternative (mean over shanks), so the curator can see what drives it
Fit.ComponentNames = [{'activity'} LMName];
for k = 1:numel(Fit.Alternatives)
    cs = nan(nSh,1 + nLM);
    for si = 1:nSh
        [~,cs(si,:)] = local_score(E{si},LM{si},PredAt(Hist{si},(dGrid - Fit.Alternatives(k).Offset)/Fit.Alternatives(k).Scale), ...
            LmAt(Hist{si},(dGrid - Fit.Alternatives(k).Offset)/Fit.Alternatives(k).Scale),Wsh{si});
    end
    Fit.Alternatives(k).Components = mean(cs,1,'omitnan');
end
fprintf('Global fit alternatives (scale, offset um, correlation, score with prior | %s):\n',strjoin(Fit.ComponentNames,', '))
for k = 1:numel(Fit.Alternatives)
    fprintf('   %.2f  %5.0f  %.2f  %.2f | %s\n',Fit.Alternatives(k).Scale,Fit.Alternatives(k).Offset,Fit.Alternatives(k).Corr,Fit.Alternatives(k).Score, ...
        strjoin(arrayfun(@(c) sprintf('%.2f',c),Fit.Alternatives(k).Components,'Uni',0),'  '))
end
Fit.ShankOffset = nan(1,nSh); Fit.CorrShankGlobal = nan(1,nSh); Fit.CorrAligned = nan(1,nSh);
for si = 1:nSh
    bs = Fit.Offset + (-Opt.ShankOffsetUm:Opt.OffsetStepUm:Opt.ShankOffsetUm);
    rs = arrayfun(@(b) ScoreAt(si,(dGrid - b)/Fit.Scale),bs);
    sp_ = arrayfun(@(b) local_surfpen(Fit.Scale,b,Hist(si),TopUnit(si),SurfaceObs(si),Opt),bs);
    [~,k] = max(rs + OffsetPrior(bs) - sp_); Fit.ShankOffset(si) = bs(k); Fit.CorrShankGlobal(si) = rs(k);
end

%% 4. Local warp per shank (dynamic programming over displacements at knots)
Depth2Area = cell(nSh,1);
for si = 1:nSh
    a = Fit.Scale; b0 = Fit.ShankOffset(si);
    sGlobal = (dGrid - b0)/a;
    pG = PredAt(Hist{si},sGlobal); okG = ~isnan(pG);
    beta = [pG(okG) ones(sum(okG),1)] \ E{si}(okG); % scale the prediction to the activity
    % scale each landmark's area indicator to its profile
    W = Wsh{si}; hG = LmAt(Hist{si},sGlobal); betaL = zeros(2,nLM);
    for j = find(W > 0)
        okj = okG & ~isnan(LM{si}(:,j));
        if sum(okj) < 10 || std(hG(okj,j)) == 0, W(j) = 0; continue; end
        betaL(:,j) = [hG(okj,j) ones(sum(okj),1)] \ LM{si}(okj,j);
    end
    knots = (dGrid(1):Opt.KnotUm:dGrid(end))'; if knots(end) < dGrid(end), knots(end+1) = dGrid(end); end %#ok<AGROW>
    U = (-Opt.MaxWarpUm:Opt.BinUm:Opt.MaxWarpUm)'; nU = numel(U);
    MaxStep = Opt.MaxStretchDev*Opt.KnotUm/a; % max change in displacement between knots (um)
    % cost of displacement u around each knot
    Cost = zeros(numel(knots),nU);
    for kk = 1:numel(knots)
        w = abs(dGrid - knots(kk)) <= Opt.KnotUm/2;
        for iu = 1:nU
            p = PredAt(Hist{si},sGlobal(w) + U(iu));
            res = E{si}(w) - (beta(1)*p + beta(2)); res(isnan(res)) = 1; % no histology there: penalised
            c = mean(res.^2);
            if any(W > 0)
                h = LmAt(Hist{si},sGlobal(w) + U(iu)); cw = 0;
                for j = find(W > 0)
                    L = LM{si}(w,j); known = ~isnan(L); % unknown landmark: no evidence
                    if ~any(known), continue; end
                    resL = L(known) - (betaL(1,j)*h(known,j) + betaL(2,j)); resL(isnan(resL)) = 1;
                    c = c + W(j)*mean(resL.^2); cw = cw + W(j);
                end
                c = c/(1 + cw);
            end
            Cost(kk,iu) = c + Opt.WarpPenalty*(U(iu)/Opt.MaxWarpUm)^2;
        end
    end
    % Viterbi
    Acc = Cost(1,:); Back = zeros(numel(knots),nU);
    for kk = 2:numel(knots)
        newAcc = inf(1,nU);
        for iu = 1:nU
            allowed = abs(U - U(iu)) <= MaxStep + eps;
            [v,j] = min(Acc(allowed)); idx = find(allowed);
            newAcc(iu) = v + Cost(kk,iu); Back(kk,iu) = idx(j);
        end
        Acc = newAcc;
    end
    path = zeros(numel(knots),1); [~,path(end)] = min(Acc);
    for kk = numel(knots):-1:2, path(kk-1) = Back(kk,path(kk)); end
    uKnot = U(path);
    u = interp1(knots,uKnot,dGrid,'linear','extrap');
    sAligned = sGlobal + u;
    pA = PredAt(Hist{si},sAligned);
    Fit.CorrAligned(si) = ScoreAt(si,sAligned);
    Feat{si}.PredictedGlobal = pG; Feat{si}.PredictedAligned = pA;
    Feat{si}.LandmarkAreasGlobal = hG; Feat{si}.LandmarkAreasAligned = LmAt(Hist{si},sAligned);
    Feat{si}.Warp = u; Feat{si}.ShankAbs = ShanksUsed(si);

    % Depth of every track point: invert s(d) (monotonic), beyond the recorded range with the global fit
    H = Hist{si};
    dPt = interp1(sAligned,dGrid,H.s,'linear',NaN);
    below = H.s < sAligned(1); above = H.s > sAligned(end);
    dPt(below) = dGrid(1) + (H.s(below) - sAligned(1))*a;
    dPt(above) = dGrid(end) + (H.s(above) - sAligned(end))*a;
    dBase = H.s*a + b0;
    % All in-brain track points (also outside the recorded range), so a curator who shifts the
    % alignment brings the neighbouring anatomy into the recorded range
    keep = ~strcmpi(H.acr,'Not found in brain');
    n = sum(keep);
    Depth2Area{si} = table(dPt(keep), repmat(si,n,1), repmat(ShanksUsed(si),n,1), H.acr(keep), H.col(keep), H.coords(keep,:), ...
        dBase(keep), dPt(keep) - dBase(keep), H.s(keep), 'VariableNames',{'Depth','Shank','ShankAbs','Area','Color','Coordinates','BaselineDepth','WarpUm','TrackPosUm'});
    % Geometry for the curation GUI: the tip (TipToSiteUm below the first site) should coincide with
    % the track end (TrackPosUm 0); the brain surface is the first in-brain track point
    inBrain = ~strcmpi(H.acr,'Not found in brain');
    Feat{si}.TipDepth = -Opt.TipToSiteUm;
    Feat{si}.SurfaceTrackPos = max(H.s(inBrain));
    Feat{si}.TrackSpacingUm = median(diff(H.s));
    Feat{si}.Scale = a; Feat{si}.Offset = b0;
    Feat{si}.SiteRange = [min(dGrid) max(dGrid)]; % recorded range on the probe (um)
end
Depth2Area = vertcat(Depth2Area{:});
AlignmentFeatures = Feat(:);
fprintf('Histology alignment: scale %.2f, offset %.0f um, correlation global %.2f -> aligned %s\n', ...
    Fit.Scale, Fit.Offset, Fit.CorrGlobal, mat2str(round(Fit.CorrAligned,2)))

%% Figure: activity vs predicted profile per shank, areas with global fit and aligned
figure('Name','Histology alignment','Position',[50 50 300*nSh 800]);
for si = 1:nSh
    subplot(1,nSh,si); hold on
    F = Feat{si}; pz = @(x) (x - min(x))./max(max(x)-min(x),eps);
    plot(pz(F.Activity),F.Depth,'k','LineWidth',1.5)
    plot(pz(F.PredictedGlobal)*0.9+0.05,F.Depth,'--','Color',[0.5 0.5 0.5])
    plot(pz(F.PredictedAligned)*0.9+0.05,F.Depth,'r')
    lmCol = lines(nLM + 1); lmCol = lmCol(2:end,:);
    for j = 1:nLM
        if all(isnan(F.Landmarks(:,j))), continue; end
        plot(pz(F.Landmarks(:,j)),F.Depth,'-','Color',lmCol(j,:))
        plot(F.LandmarkAreasAligned(:,j)*(0.95 - 0.03*j),F.Depth,':','Color',lmCol(j,:),'LineWidth',1.5)
    end
    T = Depth2Area(Depth2Area.Shank==si,:);
    for k = 1:height(T)
        c = sscanf(T.Color{k},'%2x%2x%2x')'/255;
        plot([1.1 1.3],[T.Depth(k) T.Depth(k)],'-','Color',c,'LineWidth',3)
        plot([1.35 1.55],[T.BaselineDepth(k) T.BaselineDepth(k)],'-','Color',c,'LineWidth',3)
    end
    xlim([0 1.6]); ylim([min(F.Depth) max(F.Depth)])
    set(gca,'XTick',[0.5 1.2 1.45],'XTickLabel',{'profiles','aligned','global'})
    title(sprintf('Shank %d: r %.2f -> %.2f',ShanksUsed(si),Fit.CorrShankGlobal(si),Fit.CorrAligned(si)))
    if si == 1
        ylabel('Depth on probe (um)')
        lg = {'activity','predicted (global)','predicted (aligned)'};
        for j = 1:nLM
            if all(isnan(F.Landmarks(:,j))), continue; end
            lg = [lg {LMName{j}, [LMName{j} ': areas (aligned)']}]; %#ok<AGROW>
        end
        legend(lg,'Location','southoutside','Interpreter','none')
    end
end
end

function p = local_pred(H,s,field)
% Value of H.(field) at histology positions s (um above the track end): NaN below the track end,
% 0 above the surface (outside the brain), else the nearest track point's value
% (one column per column of H.(field))
V = H.(field); s = s(:);
p = nan(numel(s),size(V,2));
inside = s >= 0 & s <= max(H.s);
if any(inside) && ~isempty(V)
    p(inside,:) = interp1(H.s,V,s(inside),'nearest');
end
p(s > max(H.s),:) = 0;
end

function pen = local_surfpen(a,b,Hist,TopUnit,SurfaceObs,Opt)
% Brain-surface penalty, averaged over shanks: the aligned surface (a*sSurface + b) may not be below
% the highest units; if the probe was partly outside the brain, it is expected near SurfaceObs
pen = 0; n = 0;
for si = 1:numel(Hist)
    if isnan(TopUnit(si)), continue; end
    dSurf = a*Hist{si}.sSurface + b;
    p = 0;
    if dSurf < TopUnit(si) % units outside the brain: not possible
        p = p + Opt.SurfaceWeight*min(((TopUnit(si) - dSurf)/100)^2, 9);
    end
    if ~isnan(SurfaceObs(si))
        p = p + Opt.SurfaceWeight*min(((dSurf - SurfaceObs(si))/Opt.SurfaceSdUm)^2, 9);
    end
    pen = pen + p; n = n + 1;
end
pen = pen/max(n,1);
end

function [r,comp] = local_score(E,L,p,h,w)
% Correlation of activity E with predicted profile p, combined with the correlation of each landmark
% profile L(:,j) with its area indicator h(:,j) (weights w). NaN if less than 60% of the depths have
% histology. comp: [r_activity r_landmark1 ...] (NaN for unused landmarks)
r = NaN; comp = nan(1,1 + numel(w));
ok = ~isnan(p) & ~isnan(E);
if mean(ok) <= 0.6 || std(p(ok)) == 0, return; end
r = corr(E(ok),p(ok)); comp(1) = r;
for j = find(w > 0)
    okj = ok & ~isnan(L(:,j));
    rr = 0; % areas not in range (or landmark unknown there): no evidence
    if sum(okj) >= 10 && std(h(okj,j)) > 0, rr = corr(L(okj,j),h(okj,j)); end
    r = r + w(j)*rr; comp(1 + j) = rr;
end
r = r/(1 + sum(w));
end
