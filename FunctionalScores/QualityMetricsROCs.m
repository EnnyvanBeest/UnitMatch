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

    ThisSes = find(recsesall == RecOpt(recid) & GoodId' == 1); % Only take good units of this round
    % Which of these units have a match?
    tmpORI = OriIDAll(ThisSes); % Original ID
    tmpUID = UniqueIDConversion.UniqueID(ThisSes); % UID
    [UIDind,id1] = ismember(UniqueIDConversion.UniqueID,tmpUID); % Figure out which units have a match
    countthem = histcounts(id1,[0.5:1:length(tmpUID)+0.5]);
    hasMatch = countthem>1;

    if min(qMetric.clusterID) == 1 %0-index
        qMetric.clusterID = qMetric.clusterID - 1;
    end
    %     qMclusterID = qMetric.clusterID - 1; %0-index
    [taketheseUnits,putUnithere] = ismember(qMetric.clusterID, tmpORI); % Find the same cluster IDs

    clear qtmp
    if recid == 1
        qtmp(putUnithere(putUnithere~=0),:) = qMetric(taketheseUnits, :); % Qmetrics table
        qtmp.hasMatch = hasMatch';
        qMetricAllGoodUnits = qtmp; % Add sessions together, only take good units
        VarNames = qMetricAllGoodUnits.Properties.VariableNames;
    else
        [taketheseMetrics,puthere] = ismember(qMetric.Properties.VariableNames,qMetricAllGoodUnits.Properties.VariableNames);
        qtmp(putUnithere(putUnithere~=0),puthere(puthere~=0)) = qMetric(taketheseUnits, taketheseMetrics);
        qtmp.hasMatch = hasMatch';
        qMetricAllGoodUnits = cat(1, qMetricAllGoodUnits, qtmp);
    end
end

%% Units that formed a match with any other unit (across recordings)
% rowidx = find(MatchTable.MatchProb > UMparam.ProbabilityThreshold); % Row in table
% ClusIDPairs = [MatchTable.ID1(rowidx) MatchTable.ID2(rowidx)];
% recSesIDPairs =  [MatchTable.RecSes1(rowidx) MatchTable.RecSes2(rowidx)];
% 
% SameUnit = ridx == cidx;
% Pairs = cat(2, ridx(~SameUnit), cidx(~SameUnit));
% UnitsWithAMatch = unique(Pairs(:));
    
UnitsWithAMatch = find(qMetricAllGoodUnits.hasMatch);
UnitsWithoutAMatch = find(~qMetricAllGoodUnits.hasMatch);
%% For every qParams, make an ROC
AUCqParams = struct;
AUCqParams.qParamNames = qMetricAllGoodUnits.Properties.VariableNames;
AUCqParams.AUCMvNM = nan(1,length(AUCqParams.qParamNames));
if isempty(UnitsWithAMatch)
  save(fullfile(SaveDir,'qMetricAUCs'),'AUCqParams')
    return
end

figure('name', 'ROCs with Quality Metrics')
for qmid = 1:size(qMetricAllGoodUnits, 2) - 1

    labels = [ones(1, numel(UnitsWithAMatch)), zeros(1, numel(UnitsWithoutAMatch))];
    scores = [table2array(qMetricAllGoodUnits(UnitsWithAMatch, qmid))', table2array(qMetricAllGoodUnits(UnitsWithoutAMatch, qmid))'];
    if length(unique(scores)) == 1
        continue
    end
    subplot(ceil(sqrt(size(qMetricAllGoodUnits, 2))), round(sqrt(size(qMetricAllGoodUnits, 2))), qmid)

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