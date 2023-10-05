function QualityMetricsROCs(SaveDir)

TmpFile = load(fullfile(SaveDir, 'UnitMatch.mat')); % Access saved file
UMparam = TmpFile.UMparam; % Extract parameters
UMparam.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows

MatchTable = TmpFile.MatchTable; % Load Matchtable

% Extract cluster information
UniqueIDConversion = TmpFile.UniqueIDConversion;
if UMparam.GoodUnitsOnly
    GoodId = logical(UniqueIDConversion.GoodID);
else
    GoodId = true(1, length(UniqueIDConversion.GoodID));
end
UniqueID = UniqueIDConversion.UniqueID(GoodId);
OriID = UniqueIDConversion.OriginalClusID(GoodId);
OriIDAll = UniqueIDConversion.OriginalClusID;
recses = UniqueIDConversion.recsesAll(GoodId);
recsesall = UniqueIDConversion.recsesAll;

AllKSDir = UMparam.KSDir; %original KS Dir
nclus = length(UniqueID);
RecOpt = unique(recses);
nRec = length(RecOpt);
% Load qparams
for recid = 1:nRec
    if UMparam.RunPyKSChronicStitched
        d = dir(fullfile(AllKSDir{1}, '**', 'templates._bc_qMetrics.parquet'));
    else
        d = dir(fullfile(AllKSDir{RecOpt(recid)}, '**', 'templates._bc_qMetrics.parquet'));
    end
    qMetricsPath = d.folder;
    [~, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(qMetricsPath);

    ThisIdx = find(GoodId(recsesall == RecOpt(recid))); % Only take good units of this round
    qMclusterID = qMetric.clusterID - 1; %0-index
    tmpORI = OriIDAll(recsesall == RecOpt(recid));
    [~,UMidx] = ismember(qMclusterID, tmpORI(ThisIdx)); % Find the same cluster IDs
    if recid == 1
        qMetricAllGoodUnits(UMidx(UMidx~=0),:) = qMetric(ThisIdx, :); % Add sessions together, only take good units
        VarNames = qMetricAllGoodUnits.Properties.VariableNames;
    else
        [takethese,puthere] = ismember(qMetric.Properties.VariableNames,qMetricAllGoodUnits.Properties.VariableNames);
        qtmp(UMidx(UMidx~=0),puthere(puthere~=0)) = qMetric(ThisIdx, takethese);
        qMetricAllGoodUnits = cat(1, qMetricAllGoodUnits, qtmp);
    end
end

%% Units that formed a match with any other unit (across recordings)
MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
[ridx, cidx] = find(MatchProb > UMparam.ProbabilityThreshold);
SameUnit = ridx == cidx;
Pairs = cat(2, ridx(~SameUnit), cidx(~SameUnit));
UnitsWithAMatch = unique(Pairs(:));
if isempty(UnitsWithAMatch)
    %% save
    AUCqParams = struct;
    AUCqParams.qParamNames = qMetricAllGoodUnits.Properties.VariableNames;
    AUCqParams.AUCMvNM = nan(1,length(AUCqParams.qParamNames));
    save(fullfile(SaveDir,'qMetricAUCs'),'AUCqParams')
    return

end
    
UnitsWithoutAMatch = 1:nclus;
UnitsWithoutAMatch(ismember(UnitsWithoutAMatch, UnitsWithAMatch)) = [];

%% For every qParams, make an ROC
AUCqParams = struct;
AUCqParams.qParamNames = qMetricAllGoodUnits.Properties.VariableNames;
AUCqParams.AUCMvNM = nan(1,length(AUCqParams.qParamNames));
figure('name', 'ROCs with Quality Metrics')
for qmid = 1:size(qMetricAllGoodUnits, 2)
    subplot(ceil(sqrt(size(qMetricAllGoodUnits, 2))), round(sqrt(size(qMetricAllGoodUnits, 2))), qmid)

    labels = [ones(1, numel(UnitsWithAMatch)), zeros(1, numel(UnitsWithoutAMatch))];
    scores = [table2array(qMetricAllGoodUnits(UnitsWithAMatch, qmid))', table2array(qMetricAllGoodUnits(UnitsWithoutAMatch, qmid))'];
    [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
    h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);
    hold on
    plot([0, 1], [0, 1], 'k--')
    xlabel('False positive rate')
    ylabel('True positive rate')
    title([qMetricAllGoodUnits.Properties.VariableNames{qmid}, ', AUC = ', num2str(AUC3)])
    makepretty
    axis square
    drawnow %Something to look at while ACG calculations are ongoing

    AUCqParams.AUCMvNM(qmid) = AUC3;
end


%% save
save(fullfile(SaveDir,'qMetricAUCs'),'AUCqParams')
set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
saveas(gcf, fullfile(SaveDir, 'QMetricsROCs.fig'))
saveas(gcf, fullfile(SaveDir, 'QMetricsROCs.png'))