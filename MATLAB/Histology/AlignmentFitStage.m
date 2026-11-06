function [Depth2Area, AlignmentFeatures, Fit] = AlignmentFitStage(S, Fixed, Landscape)
% [Depth2Area, AlignmentFeatures, Fit] = AlignmentFitStage(S, Fixed, Landscape)
% Fitting steps of alignatlasdata_constrained (3. global fit, per-shank offset; 4. local warp; and
% the figure) on its prepared profiles S (its 4th output, AlignmentState). See there for the model.
% Fixed (optional): [Scale Offset] of the global fit, e.g. from JointChronicAlignment; the per-shank
%   offsets and the warp are then fitted around it as usual
% Landscape (optional): Fit.Landscape of an earlier fit of the same S (then not computed again)
% AlignmentFeatures come without the CSD images (alignatlasdata_constrained adds them).
if nargin < 2, Fixed = []; end
if nargin < 3, Landscape = []; end
Opt = S.Opt; E = S.E; LM = S.LM; Hist = S.Hist; Wsh = S.Wsh; Pen = S.Pen; dGrid = S.dGrid;
TopUnit = S.TopUnit; SurfaceObs = S.SurfaceObs; Feat = S.Feat; ShanksUsed = S.ShanksUsed;
LMName = S.LMName; Fit = S.Fit; nSh = numel(ShanksUsed); nLM = numel(LMName);

SurfacePen = @(a,b) local_surfpen(a,b,Hist,TopUnit,SurfaceObs,Opt);
% Predicted value at probe depths d for a mapping s(d); NaN below the track end, 0 above the surface
PredAt = @(H,s) local_pred(H,s,'pred');
LmAt   = @(H,s) local_pred(H,s,'lm');
TolAt  = @(H,s) local_pred(H,s,'tol');
% Score of a mapping s(d) for one shank: correlation of the activity with the predicted profile,
% combined with the correlation of each landmark with its areas
ScoreAt = @(si,s) local_score(E{si},LM{si},PredAt(Hist{si},s),LmAt(Hist{si},s),Wsh{si},Pen{si},TolAt(Hist{si},s));
OffsetPrior = @(b) -Opt.PriorWeight*min(((b + Opt.TipToSiteUm)/Opt.OffsetPriorSdUm).^2, Opt.PriorCapSD^2);

%% 3. Global fit (scale and offset shared by the shanks), then per-shank offset
if isempty(Landscape)
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
    Landscape = struct('Scales',Scales,'Offsets',Offsets,'Corr',C,'CorrWithPrior',Cp);
end
Scales = Landscape.Scales; Offsets = Landscape.Offsets; C = Landscape.Corr; Cp = Landscape.CorrWithPrior;
if isempty(Fixed)
    [~,best] = max(Cp(:)); [ia,ib] = ind2sub(size(Cp),best);
    Fit.Scale = Scales(ia); Fit.Offset = Offsets(ib); Fit.CorrGlobal = C(ia,ib);
else
    Fit.Scale = Fixed(1); Fit.Offset = Fixed(2);
    Fit.CorrGlobal = mean(arrayfun(@(si) ScoreAt(si,(dGrid - Fit.Offset)/Fit.Scale),1:nSh),'omitnan');
end
Fit.Landscape = Landscape;
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
        s = (dGrid - Fit.Alternatives(k).Offset)/Fit.Alternatives(k).Scale;
        [~,cs(si,:)] = local_score(E{si},LM{si},PredAt(Hist{si},s),LmAt(Hist{si},s),Wsh{si},Pen{si},TolAt(Hist{si},s));
    end
    Fit.Alternatives(k).Components = mean(cs,1,'omitnan');
end
fprintf('Global fit alternatives (scale, offset um, correlation, score with prior | %s):\n',strjoin(Fit.ComponentNames,', '))
for k = 1:numel(Fit.Alternatives)
    fprintf('   %.2f  %5.0f  %.2f  %.2f | %s\n',Fit.Alternatives(k).Scale,Fit.Alternatives(k).Offset,Fit.Alternatives(k).Corr,Fit.Alternatives(k).Score, ...
        strjoin(arrayfun(@(c) sprintf('%.2f',c),Fit.Alternatives(k).Components,'Uni',0),'  '))
end
if ~isempty(Fixed), fprintf('Global fit fixed at scale %.2f, offset %.0f um (correlation %.2f)\n',Fit.Scale,Fit.Offset,Fit.CorrGlobal); end
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
ax = gobjects(1,nSh);
yl = [min(cellfun(@(F) min(F.Depth),Feat)) max(cellfun(@(F) max(F.Depth),Feat))]; % same depth axis on all shanks
for si = 1:nSh
    ax(si) = subplot(1,nSh,si); hold on
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
    xlim([0 1.6]); ylim(yl)
    set(gca,'XTick',[0.5 1.2 1.45],'XTickLabel',{'profiles','aligned','global'})
    title(sprintf('Shank %d: r %.2f -> %.2f',ShanksUsed(si),Fit.CorrShankGlobal(si),Fit.CorrAligned(si)))
    if si == 1
        ylabel('Depth on probe (um)')
        lg = {'activity','predicted (global)','predicted (aligned)'};
        for j = 1:nLM
            if all(isnan(F.Landmarks(:,j))), continue; end
            lg = [lg {LMName{j}, [LMName{j} ': areas (aligned)']}]; %#ok<AGROW>
        end
        hLg = legend(lg,'Location','southoutside','Interpreter','none');
    end
end
% The legend shrinks the first axes: fix the legend where it is and give all shanks that height
% so the depths line up (drawnow: the legend layout is applied lazily)
drawnow
p1 = ax(1).Position; pLg = hLg.Position;
hLg.Location = 'none'; hLg.Position = pLg;
for si = 1:nSh
    ax(si).Position([2 4]) = p1([2 4]);
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

function [r,comp] = local_score(E,L,p,h,w,pen,tol)
% Correlation of activity E with predicted profile p, combined with the correlation of each landmark
% profile L(:,j) with its area indicator h(:,j) (weights w). NaN if less than 60% of the depths have
% histology. comp: [r_activity r_landmark1 ...] (NaN for unused landmarks)
% pen(j) (optional, default 0): the value for landmark j when none of its areas is in the recorded
% range (confirmed hippocampal CSD landmarks: a mismatch, not just no evidence), unless areas where
% it may also be (tol(:,j), e.g. ripples in presubiculum) are in range: then 0
if nargin < 6 || isempty(pen), pen = zeros(size(w)); end
if nargin < 7 || isempty(tol), tol = zeros(size(h)); end
r = NaN; comp = nan(1,1 + numel(w));
ok = ~isnan(p) & ~isnan(E);
if mean(ok) <= 0.6 || std(p(ok)) == 0, return; end
r = corr(E(ok),p(ok)); comp(1) = r;
for j = find(w > 0)
    okj = ok & ~isnan(L(:,j));
    rr = 0; % areas not in range (or landmark unknown there): no evidence
    if sum(okj) >= 10 && std(h(okj,j)) > 0, rr = corr(L(okj,j),h(okj,j));
    elseif sum(okj) >= 10 && ~any(h(okj,j)) && ~any(tol(okj,j)), rr = pen(j); % areas not in range at all
    end
    r = r + w(j)*rr; comp(1 + j) = rr;
end
r = r/(1 + sum(w));
end
