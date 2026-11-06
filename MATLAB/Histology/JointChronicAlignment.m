function JointFit = JointChronicAlignment(AlignFiles, Opt)
% JointFit = JointChronicAlignment(AlignFiles, Opt)
% Histology alignment of a chronically implanted probe, fitted jointly over its recording days.
% AlignFiles: the days' *_HistoEphysAlignment_Auto.mat (GetHistologyOutput, constrained method), each
% fitted on its own recordings' site depths, so each already contains that day's drift.
% Joint model:
%  - one scale for all days (shrinkage/registration of the one histology track)
%  - offset of day d: b_d = B - DriftZ_d + e_d, with B common, DriftZ_d the day's drift from
%    UnitMatch (ChronicInfo.DriftZ: depth + DriftZ = depth in the reference recording; only
%    differences between days matter) and a deviation e_d (|e_d| <= Opt.MaxDevUm) with a penalty
%    Opt.DriftWeight*(e_d/Opt.DriftSdUm)^2. A day without landmarks so follows the days with strong
%    ones, while a clear misfit of the drift can still show.
%  - days without a drift estimate in the most common UnitMatch run share only the scale
% Joint score at (scale, B): the mean over days of each day's best score in its fit landscape
% (AlignmentFit.Landscape.CorrWithPrior, with its offset prior and brain-surface constraint).
% Each day is then fitted again at its scale/offset (per-shank offsets and local warp as usual,
% AlignmentFitStage) and its file updated: Depth2Area (still in that day's own probe depth; units are
% assigned by their own depth, see GetHistologyOutput), AlignmentFeatures, AlignmentFit, the figure,
% and JointFit (joint parameters; the day's offset; its independent fit). Manually curated files, and
% files flagged in the curation GUI (with a CurationDate; not the automatic landmark flag), are used for
% the estimate but their alignment is not changed: only their JointFit is saved (for the drift shift
% when a curation is propagated to the other days, and so the joint fit is not redone every run).
% Output JointFit: struct per file (empty for files that could not be used).
if nargin < 2, Opt = struct; end
Def = struct('DriftSdUm',50,'DriftWeight',0.05,'MaxDevUm',150);
for fn = fieldnames(Def)'
    if ~isfield(Opt,fn{1}), Opt.(fn{1}) = Def.(fn{1}); end
end
AlignFiles = cellstr(AlignFiles);
nF = numel(AlignFiles);
D = cell(1,nF); use = false(1,nF);
for k = 1:nF
    D{k} = load(AlignFiles{k});
    if ~isfield(D{k},'AlignmentState') || isempty(D{k}.AlignmentState) || ~isfield(D{k},'AlignmentFit') || ...
            ~isstruct(D{k}.AlignmentFit) || ~isfield(D{k}.AlignmentFit,'Landscape')
        warning('JointChronicAlignment:NoState','%s: no constrained alignment state (redo the alignment): not used',AlignFiles{k})
        continue
    end
    use(k) = true;
end
JointFit = cell(1,nF);
if sum(use) < 2
    fprintf('JointChronicAlignment: fewer than 2 days with an alignment state, nothing to fit jointly\n')
    return
end
idx = find(use);
L1 = D{idx(1)}.AlignmentFit.Landscape; Sc = L1.Scales; Of = L1.Offsets; step = Of(2) - Of(1);
for k = idx
    Lk = D{k}.AlignmentFit.Landscape;
    if ~isequal(Lk.Scales,Sc) || ~isequal(Lk.Offsets,Of)
        warning('JointChronicAlignment:Grid','%s: other fit grid: not used',AlignFiles{k}); use(k) = false;
    end
end
idx = find(use);

% Drift per day; coupled: days of the most common UnitMatch run with a drift estimate
DriftZ = nan(1,nF); UMDir = repmat({''},1,nF);
for k = idx
    if isfield(D{k},'ChronicInfo') && isstruct(D{k}.ChronicInfo)
        DriftZ(k) = D{k}.ChronicInfo.DriftZ; UMDir{k} = D{k}.ChronicInfo.UMDir;
    end
end
has = use & ~isnan(DriftZ) & ~cellfun(@isempty,UMDir);
coupled = false(1,nF);
if any(has)
    [run,~,ir] = unique(UMDir(has)); [~,mc] = max(accumarray(ir(:),1));
    coupled(has) = strcmp(UMDir(has),run{mc});
end
if sum(coupled) < 2, coupled(:) = false; end % one day: nothing to couple

% Joint score over (scale, B)
eGrid = (-round(Opt.MaxDevUm/step):round(Opt.MaxDevUm/step))*step;
ePen = Opt.DriftWeight*(eGrid/Opt.DriftSdUm).^2;
J = zeros(numel(Sc),numel(Of)); Edev = cell(1,nF);
for k = idx
    Cp = D{k}.AlignmentFit.Landscape.CorrWithPrior; Cp(isnan(Cp)) = -1; % < 60% of the depths with histology
    if coupled(k)
        best = -inf(numel(Sc),numel(Of)); bestE = zeros(numel(Sc),numel(Of));
        for ie = 1:numel(eGrid)
            ib = round((Of - DriftZ(k) + eGrid(ie) - Of(1))/step) + 1; ok = ib >= 1 & ib <= numel(Of);
            v = -inf(numel(Sc),numel(Of)); v(:,ok) = Cp(:,ib(ok)) - ePen(ie);
            bestE(v > best) = eGrid(ie); best = max(best,v);
        end
        Edev{k} = bestE;
    else
        best = repmat(max(Cp,[],2),1,numel(Of));
    end
    J = J + best;
end
J = J/numel(idx);
[jBest,i] = max(J(:)); [ia,iB] = ind2sub(size(J),i);
Scale = Sc(ia); B = Of(iB);

fprintf('Joint chronic alignment (%d days, %d coupled by drift): scale %.2f, common offset %.0f um\n',numel(idx),sum(coupled),Scale,B)
fprintf('   %-60s %8s %10s %10s %8s %8s %8s\n','file','drift','indep.off','joint off','dev','score','indep.')
for k = idx
    Fk = D{k}.AlignmentFit; Cp = Fk.Landscape.CorrWithPrior;
    if coupled(k)
        e = Edev{k}(ia,iB); b = B - DriftZ(k) + e;
    else
        [~,ib] = max(Cp(ia,:)); b = Of(ib); e = NaN;
    end
    [~,ibn] = min(abs(Of - b));
    JF = struct('Files',{AlignFiles(:)},'Coupled',coupled(k),'Scale',Scale,'CommonOffset',B,'DriftZ',DriftZ(k), ...
        'Deviation',e,'Offset',b,'Score',Cp(ia,ibn),'JointScore',jBest, ...
        'IndependentScale',Fk.Scale,'IndependentOffset',Fk.Offset,'IndependentScore',max(Cp(:)),'Opt',Opt);
    JointFit{k} = JF;
    fprintf('   %-60s %8.0f %10.0f %10.0f %8.0f %8.2f %8.2f\n',local_short(AlignFiles{k}),DriftZ(k),Fk.Offset,b,e,JF.Score,JF.IndependentScore)
    curated = isfield(D{k},'ManuallyCurated') && D{k}.ManuallyCurated;
    touched = isfield(D{k},'CurationDate') && ~isempty(D{k}.CurationDate); % accepted or flagged in the GUI
    if curated || touched
        fprintf('      manually curated/flagged: alignment not changed (JointFit saved)\n')
        S = struct('JointFit',JF);
        save(AlignFiles{k},'-struct','S','-append')
        continue
    end
    % Refit this day at the joint scale and its offset; keep its CSD images
    [Depth2Area, AlignmentFeatures, AlignmentFit] = AlignmentFitStage(D{k}.AlignmentState,[Scale b],Fk.Landscape);
    for si = 1:min(numel(AlignmentFeatures),numel(D{k}.AlignmentFeatures))
        if isfield(D{k}.AlignmentFeatures{si},'CSD'), AlignmentFeatures{si}.CSD = D{k}.AlignmentFeatures{si}.CSD; end
    end
    S = struct('Depth2Area',Depth2Area,'AlignmentFit',AlignmentFit,'DepthIsRaw',true,'JointFit',JF);
    S.AlignmentFeatures = AlignmentFeatures; % a cell: not through struct()
    save(AlignFiles{k},'-struct','S','-append')
    [fp,fn] = fileparts(AlignFiles{k});
    set(gcf,'Name',sprintf('%s (joint chronic fit)',fn))
    saveas(gcf,fullfile(fp,[fn '.fig'])); saveas(gcf,fullfile(fp,[fn '.bmp'])); close(gcf)
end
end

function s = local_short(f)
% last folders and file name, for the table
p = strsplit(f,filesep); s = strjoin(p(max(1,end-3):end),filesep);
end
