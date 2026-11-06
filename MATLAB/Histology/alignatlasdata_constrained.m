function [Depth2Area, AlignmentFeatures, Fit, AlignmentState] = alignatlasdata_constrained(histinfo, trackcoordinates, sp, clusinfo, AllenCCFPath, Opt)
% [Depth2Area, AlignmentFeatures, Fit, AlignmentState] = alignatlasdata_constrained(histinfo, trackcoordinates, sp, clusinfo, AllenCCFPath, Opt)
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
% The track's area labels are cleaned of registration flicker first (local_cleanlabels: runs shorter
% than Opt.MinAreaRunUm, and A > B > A with B shorter than Opt.MaxSandwichUm, are merged).
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
%  AlignmentState: the prepared profiles (steps 1-2), to fit again with AlignmentFitStage (steps 3-4),
%    e.g. at a given scale/offset (JointChronicAlignment). Opt.FixedFit ([Scale Offset], default
%    empty) fixes the global fit here too.

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
%    Type 'profile' (e.g. EvokedCSDLandmarks, the visual early CSD sink -> VIS layer 4/superficial SC
%    via AreaRegex on the acronym): a value per depth instead of per unit.
%  - Ripple-triggered CSD (Opt.LFP.RippleCSD): the ripple landmark is used only where the sharp-wave
%    sink/source of CA1 is present, with the sink Opt.RippleSinkDistUm below the source (or, without
%    a CSD, where the ripple band is narrow). Confirmed
%    strongly by the CSD (source and sink >= Opt.RippleCSDStrongSD) it gets weight
%    Opt.RippleCSDWeight (as the activity profile), outside the cap, and
%    the profile becomes a Gaussian (sd Opt.RipplePyrSdUm) at the CSD source (pyramidal layer),
%    fitted to Opt.RippleAreas (default CA, prosubiculum and subiculum). Ripples also reach the
%    pre-, post- and parasubiculum (Opt.RippleTolerantAreas): no evidence there, but no penalty
%    for a missing landmark either (below; also for the fissure, which is found from the ripple CSD).
%  - Theta CSD (Opt.LFP.ThetaCSD): theta sinks/sources are confined to the hippocampus -> HIP,
%    weight Opt.ThetaCSDWeight.
%  - Hippocampal fissure: the theta-trough sink below a confirmed CA1 ripple CSD (local_fissure) ->
%    the lower (molecular-layer) edge of CA1/ProS/SUB, or the top of DG-mo below it, weight
%    Opt.FissureWeight outside the cap (as dentate spikes)
%  - Opt.CSD: CSD images (EvokedCSDLandmarks) passed on to the curation GUI (AlignmentFeatures.CSD).
% Brain surface from the data: units cannot lie outside the brain, so the aligned surface may not be
% below the highest good units (hard-ish penalty). If the top >= Opt.SurfaceGapUm of the probe has no
% good units, the probe was partly outside the brain: the surface is expected Opt.SurfaceAboveUnitsUm
% above the highest units (+- Opt.SurfaceSdUm; weight Opt.SurfaceWeight per SD^2, capped at 3 SD).
if nargin < 6, Opt = struct; end
Def = struct('BinUm',10,'SmoothUm',30,'ScaleRange',[0.8 1.3],'ScaleStep',0.01,'OffsetRangeUm',[-1000 1000], ...
    'OffsetStepUm',10,'ShankOffsetUm',100,'KnotUm',50,'MaxWarpUm',150,'MaxStretchDev',0.3,'WarpPenalty',0.5,'MarginUm',50, ...
    'TipToSiteUm',175,'OffsetPriorSdUm',300,'PriorWeight',0.05,'PriorCapSD',2,'LFP',[],'LFPWeight',0.3, ...
    'MaxLandmarkWeight',0.5,'RippleMaxWidthUm',700,'ThetaCSDWeight',0.3,'ThetaCSDMinSD',4,'RippleCSDWeight',1,'RipplePyrSdUm',75,'RippleCSDStrongSD',20,'DSWeight',0.5,'SlowCSDWeight',0, ...
    'FissureWeight',0.5,'FissureNearUm',75,'FissureMinSD',6,'AbsentLandmarkScore',-0.5,'RippleSinkDistUm',[60 450], ...
    'MinAreaRunUm',20,'MaxSandwichUm',60, ...
    'RippleAreas',{{'/375/','/484682470/','/502/'}}, ... % CA, ProS, SUB ({'/375/','/484682470/'}: without SUB)
    'RippleTolerantAreas',{{'/1037/','/1084/','/843/'}}, ... % POST, PRE, PAR
    'FixedFit',[], ...
    'SurfaceGapUm',300,'SurfaceAboveUnitsUm',100,'SurfaceSdUm',150,'SurfaceWeight',0.05, ...
    'Landmarks',struct('Name',{},'Type',{},'Shank',{},'Depth',{},'Value',{},'AreaPaths',{},'AreaRegex',{},'Weight',{}), ...
    'LandmarkSmoothUm',100,'CSD',[]);
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
LMName = {}; LMAreas = {}; LMRegex = {}; LMW = [];
LM = repmat({zeros(numel(dGrid),0)},1,nSh);
% LFP features on the same depth bins (per shank); the ripple event rate marks Ammon's horn
UseLFP = ~isempty(Opt.LFP); RipCSD = []; RipValid = false(1,nSh);
if UseLFP
    ripCol = find(contains(Opt.LFP.Names,'ripple events'),1); % event rate: not dominated by spikes, unlike ripple-band power
    if isempty(ripCol), ripCol = find(contains(Opt.LFP.Names,'ripple'),1); end
    % areas (Opt.RippleAreas): Ammon's horn and prosubiculum (continuous with the CA1 pyramidal layer)
    % and, by default, subiculum (CA1 ripples propagate into it; subicular bursts ride on them)
    LMName{end+1} = 'LFP ripple events'; LMAreas{end+1} = Opt.RippleAreas; LMRegex{end+1} = ''; LMW(end+1) = Opt.LFPWeight;
    if isfield(Opt.LFP,'RippleCSD'), RipCSD = Opt.LFP.RippleCSD; end
    % A CA1 sharp wave has its radiatum sink Opt.RippleSinkDistUm below the pyramidal-layer source;
    % a source/sink pair closer together (or further apart) is not one (e.g. a noisy channel)
    for k = 1:numel(RipCSD)
        dSS = RipCSD(k).SourceDepth - RipCSD(k).SinkDepth;
        if RipCSD(k).Valid && ~(dSS >= Opt.RippleSinkDistUm(1) && dSS <= Opt.RippleSinkDistUm(2))
            RipCSD(k).Valid = false;
            Fit.LandmarkNotes{end+1} = sprintf('Shank %d: ripple CSD sink %.0f um below its source (source %.1f SD at %.0f um, sink %.1f SD): not a CA1 sharp wave, not used', ...
                RipCSD(k).Shank,dSS,RipCSD(k).SourceZ,RipCSD(k).SourceDepth,RipCSD(k).SinkZ);
            fprintf('%s\n',Fit.LandmarkNotes{end})
        end
    end
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        ch = Opt.LFP.Shank == ShanksUsed(si);
        if ~any(ch), continue; end
        b = discretize(Opt.LFP.Depth(ch),Edges); V = Opt.LFP.Values(ch,:); okb = ~isnan(b);
        Vb = nan(numel(dGrid),size(V,2));
        for k = 1:size(V,2), Vb(:,k) = accumarray(b(okb),V(okb,k),[numel(dGrid) 1],@(x) mean(x,'omitnan'),NaN); end % bad channels: NaN
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
            % Ripple-triggered CSD (LFPDepthFeatures): a sharp-wave sink below a source at the
            % ripples identifies CA1; without it the events are not CA1. It overrides the width rule.
            rc = []; if ~isempty(RipCSD), rc = RipCSD([RipCSD.Shank] == ShanksUsed(si)); end
            if ~isempty(rc)
                [~,best] = max([rc.nEvents]); rc = rc(best);
                if rc.Valid && min(rc.SourceZ,-rc.SinkZ) >= Opt.RippleCSDStrongSD
                    % strong sharp-wave CSD: its source marks the CA1 pyramidal layer: a narrow profile there
                    LM{si}(:,end) = Z(exp(-0.5*((dGrid - rc.SourceDepth)/Opt.RipplePyrSdUm).^2)); RipValid(si) = true;
                    Fit.LandmarkNotes{end+1} = sprintf('Shank %d: strong ripple-triggered CSD: CA1 pyramidal layer at %.0f um (source %.1f SD, sink %.1f SD at %.0f um; %d events)', ...
                        ShanksUsed(si),rc.SourceDepth,rc.SourceZ,rc.SinkZ,rc.SinkDepth,rc.nEvents);
                elseif rc.Valid
                    LM{si}(:,end) = r; % weak: the ripple rate profile, normal weight
                    Fit.LandmarkNotes{end+1} = sprintf('Shank %d: ripple-triggered CSD confirms CA1 (source %.1f SD at %.0f um, sink %.1f SD at %.0f um; %d events)', ...
                        ShanksUsed(si),rc.SourceZ,rc.SourceDepth,rc.SinkZ,rc.SinkDepth,rc.nEvents);
                else
                    Fit.LandmarkNotes{end+1} = sprintf('Shank %d: ripple events without a CA1 sharp-wave CSD (source %.1f SD, sink %.1f SD; %d events): not used', ...
                        ShanksUsed(si),rc.SourceZ,rc.SinkZ,rc.nEvents);
                end
                fprintf('%s\n',Fit.LandmarkNotes{end})
            elseif widthUm > Opt.RippleMaxWidthUm
                Fit.LandmarkNotes{end+1} = sprintf('Shank %d: ripple events in a %.0f um wide band (peak %.0f um): not used (not specific to CA1)', ...
                    ShanksUsed(si),widthUm,dGrid(ipk));
                fprintf('%s\n',Fit.LandmarkNotes{end})
            else
                LM{si}(:,end) = r;
            end
        end
    end
end
% Theta CSD (LFPDepthFeatures): the CSD removes volume conduction, so theta sinks/sources are
% confined to the hippocampus (CA and DG, Allen 'HIP' /1080/). Landmark: the theta CSD strength per
% depth (RMS over the cycle, in noise SDs); used where it is clear (max >= Opt.ThetaCSDMinSD)
if UseLFP && isfield(Opt.LFP,'ThetaCSD') && ~isempty(Opt.LFP.ThetaCSD)
    LMName{end+1} = 'theta CSD'; LMAreas{end+1} = {'/1080/'}; LMRegex{end+1} = ''; LMW(end+1) = Opt.ThetaCSDWeight;
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        tc_ = Opt.LFP.ThetaCSD([Opt.LFP.ThetaCSD.Shank] == ShanksUsed(si));
        if isempty(tc_), continue; end
        [~,best] = max([tc_.nEvents]); tc_ = tc_(best);
        rmsZ = sqrt(mean((tc_.CSD/max(tc_.Noise,eps)).^2,2,'omitnan'));
        if max(rmsZ) < Opt.ThetaCSDMinSD
            Fit.LandmarkNotes{end+1} = sprintf('Shank %d: no clear theta CSD (max %.1f SD): not used',ShanksUsed(si),max(rmsZ));
            continue
        end
        b = discretize(tc_.Depth,Edges); okb = ~isnan(b) & ~isnan(rmsZ);
        prof = accumarray(b(okb),rmsZ(okb),[numel(dGrid) 1],@mean,NaN);
        LM{si}(:,end) = Z(prof);
        Feat{si}.Values = [Feat{si}.Values prof];
        Feat{si}.Names = [Feat{si}.Names {'theta CSD strength'}];
        Fit.LandmarkNotes{end+1} = sprintf('Shank %d: theta CSD up to %.1f SD, strongest at %.0f um',ShanksUsed(si),max(rmsZ),tc_.Depth(find(rmsZ == max(rmsZ),1)));
        fprintf('%s\n',Fit.LandmarkNotes{end})
    end
end
% Point landmarks from the LFP CSDs (LFPDepthFeatures), where confirmed (Valid): a Gaussian (sd
% given) at the depth, fitted to the areas:
%  - dentate spikes: the reference channel (largest dentate spikes) -> DG hilus/granule layer
%  - slow waves: the strongest deep sink -> isocortex layer 5 (Opt.SlowCSDWeight 0: shown, not
%    fitted; on EB001 the strongest slow-wave generator was hippocampal despite the exclusions)
%  - hippocampal fissure: the theta-trough sink below a confirmed CA1 ripple source (local_fissure).
%    It is the perforant-path input layer (CA1 lacunosum-moleculare, or the subicular molecular
%    layer where the track crosses subiculum; DG need not be below) -> the lower edge of the
%    CA1/ProS/SUB band and the top of DG-mo just below it (LMNear, within Opt.FissureNearUm).
%    Evidence for DG itself comes from dentate spikes.
PointLM = {'DentateCSD','dentate spikes (DG hilus)','RefDepth',{},'^DG-(po|sg)$',Opt.DSWeight,75; ...
           'SlowCSD','slow-wave sink (cortex L5)','SinkDepth',{'/315/'},'5$',Opt.SlowCSDWeight,100; ...
           'FissureCSD','hippocampal fissure (theta CSD)','SinkDepth',{},'^(CA1|CA2|CA3|ProS|SUB|DG-mo)$',Opt.FissureWeight,50};
LMNear = {}; % per landmark (optional): {partner regex A, partner regex B, um}: only points of A near B or B near A
if UseLFP && isfield(Opt.LFP,'ThetaCSD') && ~isempty(Opt.LFP.ThetaCSD) && ~isempty(RipCSD)
    Opt.LFP.FissureCSD = local_fissure(Opt.LFP.ThetaCSD,RipCSD,Opt);
end
for pl = 1:size(PointLM,1)
    if ~UseLFP || ~isfield(Opt.LFP,PointLM{pl,1}) || isempty(Opt.LFP.(PointLM{pl,1})), continue; end
    PC = Opt.LFP.(PointLM{pl,1});
    LMName{end+1} = PointLM{pl,2}; LMAreas{end+1} = PointLM{pl,4}; LMRegex{end+1} = PointLM{pl,5}; LMW(end+1) = PointLM{pl,6}; %#ok<AGROW>
    if strcmp(PointLM{pl,1},'FissureCSD'), LMNear{numel(LMName)} = {'^(CA1|CA2|CA3|ProS|SUB)$','^DG-mo$',Opt.FissureNearUm}; end %#ok<AGROW>
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        pc = PC([PC.Shank] == ShanksUsed(si));
        if isempty(pc), continue; end
        [~,best] = max([pc.nEvents]); pc = pc(best);
        if ~pc.Valid
            Fit.LandmarkNotes{end+1} = sprintf('Shank %d: %s not confirmed (%d events)',ShanksUsed(si),PointLM{pl,2},pc.nEvents);
            continue
        end
        d0 = pc.(PointLM{pl,3});
        LM{si}(:,end) = Z(exp(-0.5*((dGrid - d0)/PointLM{pl,7}).^2));
        Fit.LandmarkNotes{end+1} = sprintf('Shank %d: %s at %.0f um (%d events)',ShanksUsed(si),PointLM{pl,2},d0,pc.nEvents);
        fprintf('%s\n',Fit.LandmarkNotes{end})
    end
end
% Other landmarks (Opt.Landmarks), by Type:
%  'fraction' (default; e.g. TuningLandmarks): Value 0/1 per unit; fraction per depth, smoothed over
%    Opt.LandmarkSmoothUm; unknown where < 2 units nearby
%  'profile' (e.g. EvokedCSDLandmarks): Value per depth sample; mean per depth bin (smoothed over
%    Opt.SmoothUm), unknown where there are no samples
% Areas: Allen structure id path tokens (AreaPaths) and/or a regular expression on the acronym (AreaRegex)
SmL = @(x) smoothdata(x,1,'gaussian',max(3,round(5*Opt.LandmarkSmoothUm/Opt.BinUm)));
MinDens = 2/(2.5*Opt.LandmarkSmoothUm/Opt.BinUm); % ~2 units within +- 1 sigma
for k = 1:numel(Opt.Landmarks)
    L = Opt.Landmarks(k);
    typ = 'fraction'; if isfield(L,'Type') && ~isempty(L.Type), typ = L.Type; end
    rx = ''; if isfield(L,'AreaRegex'), rx = L.AreaRegex; end
    LMName{end+1} = L.Name; LMAreas{end+1} = cellstr(L.AreaPaths); LMRegex{end+1} = rx; LMW(end+1) = L.Weight; %#ok<AGROW>
    for si = 1:nSh
        LM{si}(:,end+1) = NaN;
        on = L.Shank(:) == ShanksUsed(si);
        if sum(on) < 10 || ~any(L.Value(on)), continue; end
        if strcmp(typ,'profile')
            b = discretize(L.Depth(on),Edges); v = L.Value(on); okb = ~isnan(b);
            prof = accumarray(b(okb),v(okb),[numel(dGrid) 1],@mean,NaN);
            known = ~isnan(prof);
            prof = Sm(fillmissing(prof,'linear','EndValues','nearest')); prof(~known) = NaN;
        else
            cnt = SmL(histcounts(L.Depth(on),Edges)');
            pos = SmL(histcounts(L.Depth(on & L.Value(:) > 0),Edges)');
            prof = pos./cnt; prof(cnt < MinDens) = NaN;
        end
        if mean(~isnan(prof)) < 0.2, continue; end
        LM{si}(:,end) = Z(prof);
        Feat{si}.Values = [Feat{si}.Values prof];
        if strcmp(typ,'profile'), Feat{si}.Names = [Feat{si}.Names {L.Name}]; else, Feat{si}.Names = [Feat{si}.Names {['tuning: ' L.Name]}]; end
    end
end
% CSD images for the curation GUI (ripple-triggered from Opt.LFP, stimulus-evoked from Opt.CSD)
for si = 1:nSh
    Feat{si}.CSD = struct('Name',{},'Depth',{},'t',{},'CSD',{});
    if UseLFP && ~isempty(RipCSD)
        for c = RipCSD([RipCSD.Shank] == ShanksUsed(si))
            Feat{si}.CSD(end+1) = struct('Name',sprintf('ripple CSD (%d ev)',c.nEvents),'Depth',c.Depth,'t',c.t,'CSD',c.CSD/max(c.Noise,eps));
        end
    end
    lfpCSD = {'ThetaCSD','theta CSD (%d cyc)'; 'DentateCSD','dentate spike CSD (%d ev)'; 'SlowCSD','slow-wave CSD (%d waves)'};
    for kc = 1:size(lfpCSD,1)
        if ~UseLFP || ~isfield(Opt.LFP,lfpCSD{kc,1}) || isempty(Opt.LFP.(lfpCSD{kc,1})), continue; end
        CC = Opt.LFP.(lfpCSD{kc,1});
        for c = CC([CC.Shank] == ShanksUsed(si))
            Feat{si}.CSD(end+1) = struct('Name',sprintf(lfpCSD{kc,2},c.nEvents),'Depth',c.Depth,'t',c.t,'CSD',c.CSD/max(c.Noise,eps));
        end
    end
    if isstruct(Opt.CSD) && ~isempty(Opt.CSD)
        for c = Opt.CSD([Opt.CSD.Shank] == ShanksUsed(si))
            Feat{si}.CSD(end+1) = struct('Name',c.Name,'Depth',c.Depth,'t',c.t,'CSD',c.CSD);
        end
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
    distV = vn(find(contains(vn,'Distance','IgnoreCase',true),1));
    acrV  = vn{find(contains(vn,'Acronym','IgnoreCase',true),1)};
    acr = H.(acrV);
    if ~isempty(distV)
        tTrack = double(H.(distV{1}));
    else
        % older brainreg-segment tracks (Position, Region ID, ...): no distance column. It is the
        % path length along the track points (um; as the 'Distance from first position' of newer files)
        nC = min(height(H),size(trackcoordinates{si},1));
        tTrack = [0; cumsum(sqrt(sum(diff(double(trackcoordinates{si}(1:nC,:))).^2,2)))];
    end
    nPt = min(numel(tTrack),size(trackcoordinates{si},1));
    tTrack = tTrack(1:nPt); acr = acr(1:nPt); coords = trackcoordinates{si}(1:nPt,:);
    acr = local_cleanlabels(acr,tTrack,Opt);
    [~,row] = ismember(lower(strrep(acr,'/','')),stAcr);
    pred = zeros(nPt,1); lm = zeros(nPt,nLM); col = repmat({'000000'},nPt,1);
    for k = 1:nPt
        if row(k) == 0, continue; end % 'Not found in brain' etc.
        col{k} = stCol{row(k)};
        pth = st.structure_id_path{row(k)};
        for j = 1:nLM % e.g. Ammon's horn for ripples, VIS layer 4 for the visual CSD sink
            % path tokens and/or acronym pattern; with both, both must match (e.g. isocortex AND layer 5)
            inPath = any(cellfun(@(a) contains(pth,a),LMAreas{j}));
            inRx = ~isempty(LMRegex{j}) && ~isempty(regexp(acr{k},LMRegex{j},'once'));
            if ~isempty(LMAreas{j}) && ~isempty(LMRegex{j}), lm(k,j) = inPath && inRx; else, lm(k,j) = inPath || inRx; end
        end
        if contains(pth,'/1009/') || contains(pth,'/73/') || strcmpi(stAcr{row(k)},'root'), continue; end % fibre tract, ventricle
        efr = stEFR(row(k)); if isnan(efr) || efr <= 0, efr = 1; end
        pred(k) = 1 + 0.25*log10(efr);
    end
    % Edge landmarks (LMNear = {A, B, um}): the lower edge of area A along the track (points of A with
    % a point outside A in the brain at most um deeper), plus points of B at most um below A, e.g. the
    % molecular-layer side of CA1/ProS/SUB, and the top of DG-mo where DG lies below
    inBrain = ~strcmpi(acr,'Not found in brain') & row > 0;
    for j = find(~cellfun(@isempty,LMNear(1:min(end,nLM))))
        nr = LMNear{j};
        inA = ~cellfun(@isempty,regexp(acr,nr{1},'once')); inB = ~cellfun(@isempty,regexp(acr,nr{2},'once'));
        edge = false(nPt,1);
        dt_ = tTrack(:)' - tTrack(:); % (k,m): how much deeper point m is than point k
        if any(inA)
            other = inBrain(:)' & ~inA(:)';
            edge = edge | (inA(:) & any(dt_ > 0 & dt_ <= nr{3} & other,2));
            edge = edge | (inB(:) & any(-dt_ > 0 & -dt_ <= nr{3} & inA(:)',2)); % B just below A
        end
        lm(:,j) = lm(:,j) & edge;
    end
    % Areas where the ripples (and the fissure found from them) may be although not listed: no penalty
    tol = zeros(nPt,nLM);
    jt = find(ismember(LMName,{'LFP ripple events','hippocampal fissure (theta CSD)'}));
    okr = row > 0;
    if ~isempty(jt) && any(okr)
        tol(okr,jt) = repmat(cellfun(@(p) any(cellfun(@(a) contains(p,a),Opt.RippleTolerantAreas)),st.structure_id_path(row(okr))),1,numel(jt));
    end
    sTrack = max(tTrack) - tTrack; % distance above the track end (um, atlas)
    [sTrack,o] = sort(sTrack,'ascend'); % from the track end upwards
    Hist{si} = struct('s',sTrack,'pred',pred(o),'lm',lm(o,:),'tol',tol(o,:),'acr',{acr(o)},'col',{col(o)},'coords',coords(o,:));
    Hist{si}.sSurface = max(sTrack(~strcmpi(acr(o),'Not found in brain')));
end
% A landmark that is clearly present while its areas are nowhere on the track: the track (or the
% landmark's area list) may be wrong; it cannot help the fit then
for k = 1:numel(Opt.Landmarks)
    j = find(strcmp(LMName,Opt.Landmarks(k).Name),1);
    if Opt.Landmarks(k).Weight <= 0, continue; end % display only
    isProfile = isfield(Opt.Landmarks,'Type') && strcmp(Opt.Landmarks(k).Type,'profile'); % only passed on when clear (e.g. CSD sink > threshold)
    for si = 1:nSh
        on = Opt.Landmarks(k).Shank(:) == ShanksUsed(si);
        if sum(on) < 10 || any(Hist{si}.lm(:,j)), continue; end
        if isProfile
            Fit.LandmarkWarnings{end+1} = sprintf('Shank %d: clear "%s", but none of its areas on the histology track',ShanksUsed(si),LMName{j});
        elseif mean(Opt.Landmarks(k).Value(on) > 0) >= 0.2
            Fit.LandmarkWarnings{end+1} = sprintf('Shank %d: %.0f%% of units "%s", but none of its areas on the histology track', ...
                ShanksUsed(si),100*mean(Opt.Landmarks(k).Value(on) > 0),LMName{j});
        else
            continue
        end
        warning('alignatlasdata_constrained:Landmark','%s',Fit.LandmarkWarnings{end})
    end
end
% Likewise for the confirmed hippocampal CSD landmarks (CA1 ripples, fissure, dentate spikes): their
% areas missing on the track means the track does not pass where the probe recorded (e.g. a track of
% another insertion). The score gives such a landmark no evidence (0), so it would not show otherwise.
for jn = {'LFP ripple events','hippocampal fissure (theta CSD)','dentate spikes (DG hilus)'}
    j = find(strcmp(LMName,jn{1}),1); if isempty(j), continue; end
    for si = 1:nSh
        confirmed = any(~isnan(LM{si}(:,j)));
        if strcmp(jn{1},'LFP ripple events'), confirmed = RipValid(si); end
        if confirmed && ~any(Hist{si}.lm(:,j)) && any(Hist{si}.tol(:,j))
            Fit.LandmarkNotes{end+1} = sprintf('Shank %d: confirmed "%s"; none of its areas on the histology track, but pre/post/parasubiculum (no penalty)',ShanksUsed(si),LMName{j});
            fprintf('%s\n',Fit.LandmarkNotes{end})
        elseif confirmed && ~any(Hist{si}.lm(:,j))
            Fit.LandmarkWarnings{end+1} = sprintf('Shank %d: confirmed "%s", but none of its areas on the histology track (wrong track?)',ShanksUsed(si),LMName{j});
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
% Landmark weights per shank (0 for landmarks without data on that shank)
Wsh = cellfun(@(L) LMW .* any(~isnan(L),1),LM,'Uni',0);
Wsh = cellfun(@(w) w*min(1,Opt.MaxLandmarkWeight/max(sum(w),eps)),Wsh,'Uni',0); % activity dominates
% Ripples confirmed by the ripple-triggered CSD (CA1 sharp-wave sink below a source at the ripples)
% are a precise anchor for the CA1 pyramidal layer: weight Opt.RippleCSDWeight, outside the cap
rc = find(strcmp(LMName,'LFP ripple events'),1);
if ~isempty(rc)
    for si = find(RipValid), Wsh{si}(rc) = Opt.RippleCSDWeight; end
end
% confirmed dentate spikes likewise mark one layer precisely: their weight outside the cap too
jd = find(strcmp(LMName,'dentate spikes (DG hilus)'),1);
if ~isempty(jd)
    for si = 1:nSh, if any(~isnan(LM{si}(:,jd))), Wsh{si}(jd) = Opt.DSWeight; end, end
end
% and the hippocampal fissure (theta sink below the confirmed CA1 pyramidal layer)
jf = find(strcmp(LMName,'hippocampal fissure (theta CSD)'),1);
if ~isempty(jf)
    for si = 1:nSh, if any(~isnan(LM{si}(:,jf))), Wsh{si}(jf) = Opt.FissureWeight; end, end
end
% These confirmed hippocampal CSD landmarks are tied to specific layers: a mapping that puts none of
% their areas in the recorded range is a mismatch (Opt.AbsentLandmarkScore), not just no evidence.
% (Tuning and visual-CSD landmarks also occur outside their listed areas: 0 there.)
Pen = repmat({zeros(1,nLM)},1,nSh);
for si = 1:nSh
    if ~isempty(rc) && RipValid(si), Pen{si}(rc) = Opt.AbsentLandmarkScore; end
    for jj = [jd jf], if any(~isnan(LM{si}(:,jj))), Pen{si}(jj) = Opt.AbsentLandmarkScore; end, end
end
% Steps 3-4 (AlignmentFitStage): global fit, per-shank offsets, local warp, figure. The prepared
% profiles are returned (AlignmentState) to fit again, e.g. at a given scale/offset (chronic probes)
Feat0 = cellfun(@(F) rmfield(F,'CSD'),Feat,'Uni',0); % CSD images: large, not used in the fit
AlignmentState = struct('Opt',rmfield(Opt,intersect(fieldnames(Opt),{'LFP','Landmarks','CSD'})),'E',{E},'LM',{LM}, ...
    'Hist',{Hist},'Wsh',{Wsh},'Pen',{Pen},'dGrid',dGrid,'TopUnit',TopUnit,'SurfaceObs',SurfaceObs, ...
    'Feat',{Feat0},'ShanksUsed',ShanksUsed,'LMName',{LMName},'Fit',Fit);
[Depth2Area, AlignmentFeatures, Fit] = AlignmentFitStage(AlignmentState, Opt.FixedFit);
for si = 1:nSh, AlignmentFeatures{si}.CSD = Feat{si}.CSD; end
end

function acr = local_cleanlabels(acr, t, Opt)
% Area labels along the track without the flicker of an imperfect registration: a track running
% close to a boundary jumps between the areas on either side (e.g. VISp4 > VISp6a > VISp4, or
% or > ar > or), which is not anatomy. Merged into a neighbour, shortest first:
%  - runs shorter than Opt.MinAreaRunUm (into the longer neighbour)
%  - runs shorter than Opt.MaxSandwichUm with the same area on both sides (A > B > A -> A)
% acr: labels in track order, t: distance along the track (um). 0 for both keeps the labels.
if numel(acr) < 3 || (Opt.MinAreaRunUm <= 0 && Opt.MaxSandwichUm <= 0), return; end
t = t(:); sp = median(diff(t)); if ~(sp > 0), return; end
while true
    st = [true; ~strcmp(acr(2:end),acr(1:end-1))]; s0 = find(st); s1 = [s0(2:end)-1; numel(acr)];
    len = (s1 - s0 + 1)*sp; nR = numel(s0);
    if nR < 2, return; end
    sand = false(nR,1); sand(2:end-1) = strcmp(acr(s0(1:end-2)),acr(s0(3:end)));
    cand = find(len < Opt.MinAreaRunUm | (sand & len < Opt.MaxSandwichUm));
    if isempty(cand), return; end
    [~,k] = min(len(cand)); k = cand(k);
    if k == 1, nbr = 2; elseif k == nR, nbr = nR - 1; elseif sand(k) || len(k-1) >= len(k+1), nbr = k - 1; else, nbr = k + 1; end
    acr(s0(k):s1(k)) = acr(s0(nbr));
end
end

function FC = local_fissure(ThetaCSD, RipCSD, Opt)
% Hippocampal fissure per shank: the strongest theta-trough sink below a confirmed CA1 ripple CSD.
% The overall strongest theta sink is not reliable (often above the pyramidal layer), but across
% recordings with a strong ripple CSD the theta-trough CSD has a sink 210-430 um below the ripple
% source (lacunosum-moleculare/fissure; median 300-450 um) with a source in stratum radiatum. The
% sign of the trough-locked CSD depends on where the theta reference channel is, so it is oriented
% with stratum radiatum (the ripple-CSD sink) as a source. Searched from 20 um below the ripple sink
% down to 650 um below the ripple source. Valid if the sink is >= Opt.FissureMinSD.
FC = struct('Shank',{},'SinkDepth',{},'SinkZ',{},'RadiatumZ',{},'nEvents',{},'Valid',{});
for sh = unique([ThetaCSD.Shank])
    rc = RipCSD([RipCSD.Shank] == sh & [RipCSD.Valid]);
    tc = ThetaCSD([ThetaCSD.Shank] == sh);
    if isempty(rc) || isempty(tc), continue; end
    [~,b] = max([rc.nEvents]); rc = rc(b); [~,b] = max([tc.nEvents]); tc = tc(b);
    m = mean(tc.CSD(:,abs(tc.t) <= 0.01),2)/max(tc.Noise,eps);
    rad = mean(m(abs(tc.Depth - rc.SinkDepth) <= 40),'omitnan');
    if rad < 0, m = -m; end
    win = tc.Depth <= rc.SinkDepth - 20 & tc.Depth >= rc.SourceDepth - 650 & ~isnan(m);
    if ~any(win), continue; end
    mw = m; mw(~win) = Inf; [z,i] = min(mw);
    FC(end+1) = struct('Shank',sh,'SinkDepth',tc.Depth(i),'SinkZ',z,'RadiatumZ',abs(rad),'nEvents',tc.nEvents, ...
        'Valid',z <= -Opt.FissureMinSD); %#ok<AGROW>
end
end
