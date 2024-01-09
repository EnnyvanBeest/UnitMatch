%% Load and format data

UMFile = {'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays\Non_Stitched\AL032\AllProbes\AllIMRO\UnitMatch\UnitMatch.mat'};
summaryFunctionalPlots(UMFile, 'Rank', 1)
load(UMFile{1})

% Extract cluster information
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
recSesUni = unique(recsesall);

AllKSDir = UMparam.KSDir; %original KS Dir
nclus = length(UniqueID);

nRec = length(unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID))));
RecOpt = unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID)));
EuclDist = reshape(MatchTable.EucledianDistance, nclus, nclus);
SessionSwitch = arrayfun(@(X) find(recses == X, 1, 'first'), unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID))), 'Uni', 0);
SessionSwitch(cellfun(@isempty, SessionSwitch)) = [];
SessionSwitch = [cell2mat(SessionSwitch); nclus + 1];

%% Get spikes

sp = getSpikesFromPrepData(AllKSDir);

%% Get natim resp

%%% could be made a function, reused in main UM?

% Param for processing
proc.window = [-0.3 0.5 ... % around onset
    0.0 0.5]; % around offset
proc.binSize = 0.005; % in ms
proc.smoothSize = 5; % PSTH smoothing filter
gw = gausswin(proc.smoothSize,3);
proc.smWin = gw./sum(gw);

binsOn = (proc.window(1)+proc.binSize/2):proc.binSize:proc.window(2);
binsOff =  proc.window(2)+((proc.binSize/2):proc.binSize:(proc.window(4)-proc.window(3)));
bins = [binsOn binsOff];

nRec = numel(UMparam.AllRawPaths);
nClu = nan(1,nRec); for ss = 1:nRec; nClu(ss) = numel(unique(MatchTable.ID1(MatchTable.RecSes1 == ss))); end
spikeData_cv = cell(1,2*nRec);
clusterIDs = cell(1,nRec);
spikeData = cell(1,nRec);
for ss = 1:nRec
    % Get the original binFile (also for stitched?)
    if iscell(UMparam.AllRawPaths{ss})
        binFileRef = fullfile(UMparam.AllRawPaths{ss});
    else
        binFileRef = fullfile(UMparam.AllRawPaths{ss}.folder,UMparam.AllRawPaths{ss}.name);
    end

    % Find the associated experiments
    if ispc
        exp2keep = getNatImExpRef(binFileRef);
    else
        exp2keep = [];
        % Julie Fabre: on linux the servers are mounted differently, and the above
        % function needs to be substantially changed to work. Since I
        % don't need to run it, just added an "ispc" flag.
    end

    if ~isempty(exp2keep)
        % Get the spikes
        st = sp.st(sp.RecSes == RecOpt(ss));
        clu = sp.spikeTemplates(sp.RecSes == RecOpt(ss));

        spikesAll.times = st;
        spikesAll.clusters = clu;
        spikesAll.clusterIDs = unique(clu); % follow the same units across days

        % Get the natim responses
        spikeData{ss} = getNatImResp(spikesAll,exp2keep,binFileRef,proc);
        clusterIDs{ss} = spikesAll.clusterIDs;

        % Split in two halves and subselect units
        cluIdx = ismember(clusterIDs{ss},unique(MatchTable.ID1(MatchTable.RecSes1 == ss)));
        currIdx = (ss-1)*2;
        spikeData_cv{currIdx+1} = spikeData{ss}(:,:,cluIdx,1:2:end); % odd
        spikeData_cv{currIdx+2} = spikeData{ss}(:,:,cluIdx,2:2:end); % even
    end
end

%% Get binned FR for xcorr

%%% UGLY -- SHOULD BE MADE CLEANER (+ a function, reused in main UM?)

% Get binned FR
sr = cell(1,nRec);
for rid = 1:nRec
    edges = floor(min(sp.st(sp.RecSes == RecOpt(rid)))) - UMparam.binsz / 2:UMparam.binsz:ceil(max(sp.st(sp.RecSes == RecOpt(rid)))) + UMparam.binsz / 2;
    Good_Idx = find(GoodId & recsesall' == RecOpt(rid)); % Only care about good units at this point

    % bin data to create PSTH
    sr{rid} = nan(numel(Good_Idx), numel(edges)-1);
    for uid = 1:numel(Good_Idx)
        sr{rid}(uid, :) = histcounts(sp.st(sp.spikeTemplates == OriIDAll(Good_Idx(uid)) & sp.RecSes == recsesall(Good_Idx(uid))), edges);
    end
end

% Get ref population of "matched neurons"
MatchProbability = reshape(MatchTable.MatchProb, nclus, nclus);
[r, c] = find(MatchProbability > UMparam.ProbabilityThreshold); %Find matches
Pairs = cat(2, r, c);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs, 'rows');
dayopt = [1 2];
Unit2Take = OriID;
refPopIdx = cell(1,nRec);
for did = 1:2
    % We need a group of units that is likely to be a pair across at least two days
    if did==1
        pairidx = recses(Pairs(:,1)) == RecOpt(dayopt(did)) & recses(Pairs(:,2))==RecOpt(dayopt(did+1));
        PairsTmp = Pairs(pairidx,:);
        % Only use every 'unit' once --> take the highest scoring matches
        [~,id1,~]=unique(PairsTmp(:,1),'stable');
        PairsTmp = PairsTmp(id1,:);
        [~,id1,~]=unique(PairsTmp(:,2),'stable');
        PairsTmp = PairsTmp(id1,:);
        Unit2TakeIdx = PairsTmp(:,1); % Only take each unit once
    else
        Unit2TakeIdx = [];
    end
    if did==2
        pairidx = recses(Pairs(:,2)) == RecOpt(dayopt(did)) & recses(Pairs(:,1))==RecOpt(dayopt(did-1));
        PairsTmp = Pairs(pairidx,:);
        % Only use every 'unit' once --> take the highest scoring matches
        [~,id1,~]=unique(PairsTmp(:,1),'stable');
        PairsTmp = PairsTmp(id1,:);
        [~,id1,~]=unique(PairsTmp(:,2),'stable');
        PairsTmp = PairsTmp(id1,:);
        Unit2TakeIdx = [Unit2TakeIdx; PairsTmp(:,2)];
    end
    Unit2TakeIdxAll = find(recses == RecOpt(dayopt(did)));

    % Extract the part of the correlation matrix with these
    % pairs
    sortIdx = cell2mat(arrayfun(@(x) find(Unit2Take(Unit2TakeIdxAll) == x), Unit2Take(Unit2TakeIdx), 'uni', 0));
    refPopIdx{did} = sortIdx;
end

%% Plot  

% Get pairs of sessions
sess1 = 1;
sess2 = 2;

% Find matches
MatchTable_matches = MatchTable(MatchTable.UID1 == MatchTable.UID2 & MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :);
MatchTable_bestMatches = MatchTable_matches(MatchTable_matches.ACGCorr > 0.6 & MatchTable_matches.refPopCorr > 0.8 & MatchTable_matches.natImRespCorr > 0.8,:);
matches2plt = 3; % 3
clu1 = MatchTable_bestMatches(matches2plt,:).ID1; % 98
clu2 = MatchTable_bestMatches(matches2plt,:).ID2; % 99
colMatches = [1 0 0]; 

% Find non-matches
MatchTable_nonMatches = MatchTable(MatchTable.ID1 == clu1 & MatchTable.UID1 ~= MatchTable.UID2 & MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess1, :);
MatchTable_worstNonMatches = MatchTable_nonMatches(MatchTable_nonMatches.CentroidDist > 0 & ...
    MatchTable_nonMatches.ACGCorr < 0.6 | MatchTable_nonMatches.refPopCorr < 0.7 | MatchTable_nonMatches.natImRespCorr < 0.7,:);
nonMatches2plt = 6; % 6?
clu3 = MatchTable_worstNonMatches(nonMatches2plt,:).ID2;
colNonMatches = [0 0 0.7];

% Get depth of clusters for each sessions (used for sorting)
cluList1 = unique(MatchTable(MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :).ID1);
cluDepthSess1 = nan(1,numel(cluList1));
for cluIdx = 1:numel(cluList1)
   cluDepthSess1(cluIdx) = nanmean(sp.spikeDepths(sp.RecSes == sess1 & sp.spikeTemplates == cluList1(cluIdx)));
end
cluList2 = unique(MatchTable(MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :).ID2);
cluDepthSess2 = nan(1,numel(cluList2));
for cluIdx = 1:numel(cluList2)
   cluDepthSess2(cluIdx) = nanmean(sp.spikeDepths(sp.RecSes == sess2 & sp.spikeTemplates == cluList2(cluIdx)));
end

% Sorting neurons by depths for matrices
[~,sortIdx1] = sort(cluDepthSess1,'descend');
[~,sortIdx2] = sort(cluDepthSess2,'descend');
sortCluIdx = [sortIdx1 numel(sortIdx1)+sortIdx2];

% Plot figure
figure('Position', [400 270 800 700]);

% FR
subplot(4,4,13) % matrix
FRDiffMat = reshape(MatchTable.FRDiff, nclus, nclus);
imagesc(FRDiffMat(sortCluIdx,sortCluIdx))
hold on
colormap(flipud(gray))
clim([0 30])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% ACG
p(1) = subplot(4,4,2); % non-match
hold all
idx1_1 = sp.spikeTemplates == clu1 & sp.RecSes == sess1 & sp.st < max(sp.st)/2;
[CCGClu1_1, tClu1_1] = CCGBz([double(sp.st(idx1_1)); double(sp.st(idx1_1))], [ones(size(sp.st(idx1_1), 1), 1); ...
    ones(size(sp.st(idx1_1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
idx1_2 = sp.spikeTemplates == clu1 & sp.RecSes == sess1 & sp.st > max(sp.st)/2;
[CCGClu1_2, tClu1_2] = CCGBz([double(sp.st(idx1_2)); double(sp.st(idx1_2))], [ones(size(sp.st(idx1_2), 1), 1); ...
    ones(size(sp.st(idx1_2), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
plot(tClu1_2(tClu1_2>0), CCGClu1_2(tClu1_2>0,1),'color',[0.5 0.5 0.5]);
plot(tClu1_1(tClu1_1>0), CCGClu1_1(tClu1_1>0,1),'k');
xticks([0 0.01 UMparam.ACGduration/2])
yticks([0 20])
xlabel('Time (s)')
ylabel('Firing rate (sp/s)')
set(gca,'XScale','log')
makepretty
offsetAxes

p(1) = subplot(4,4,6); % non-match
hold all
plot(tClu1_1(tClu1_1>0), CCGClu1_1(tClu1_1>0,1),'k');
[CCGClu3, tClu3] = CCGBz([double(sp.st(idx3)); double(sp.st(idx3))], [ones(size(sp.st(idx3), 1), 1); ...
    ones(size(sp.st(idx3), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
plot(tClu3(tClu3>0), CCGClu3(tClu3>0,1),'color',colNonMatches);
xticks([0 0.01 UMparam.ACGduration/2])
yticks([0 20])
xlabel('Time (s)')
ylabel('Firing rate (sp/s)')
set(gca,'XScale','log')
makepretty
offsetAxes

p(2) = subplot(4,4,10); % match
hold all
idx2 = sp.spikeTemplates == clu2 & sp.RecSes == sess2;
plot(tClu1_1(tClu1_1>0), CCGClu1_1(tClu1_1>0,1),'k');
[CCGClu2, tClu2] = CCGBz([double(sp.st(idx2)); double(sp.st(idx2))], [ones(size(sp.st(idx2), 1), 1); ...
    ones(size(sp.st(idx2), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
plot(tClu2(tClu2>0), CCGClu2(tClu2>0,1),'color',colMatches);
xticks([0 0.01 UMparam.ACGduration/2])
yticks([0 20])
xlabel('Time (s)')
ylabel('Firing rate (sp/s)')
set(gca,'XScale','log')
makepretty
offsetAxes
linkaxes(p,'y')

subplot(4,4,14) % matrix
ACGCorrMat = reshape(MatchTable.ACGCorr, nclus, nclus);
ACGCorrMat(ACGCorrMat<0.7) = 0.7;
imagesc(ACGCorrMat(sortCluIdx,sortCluIdx))
hold on
colormap('RedBlue')
clim([0.45 0.95])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% NatImCorr 
p(1) = subplot(4,4,3); % non-match
hold all
respBin = bins<proc.window(2)+0.2 & bins>0;
resp1_1 = nanmean(spikeData{recSesUni == sess1}(:,:,clusterIDs{recSesUni == sess1} == clu1,1:2:end),4);
resp1_2 = nanmean(spikeData{recSesUni == sess1}(:,:,clusterIDs{recSesUni == sess1} == clu1,2:2:end),4);
% plot timecourse
plot(bins, mean(resp1_2),'color',[0.5 0.5 0.5])
plot(bins, mean(resp1_1),'k')
% plot average response
[~,natimOrd] = sort(nanmean(resp1_1(:,respBin),2),'descend');
nIm = numel(natimOrd);
xVal = bins(end)-bins(1) + 5*proc.binSize + linspace(0,bins(end)-bins(1),nIm);
plot(xVal, nanmean(resp1_2(natimOrd,respBin),2),'color',[0.5 0.5 0.5])
plot(xVal, nanmean(resp1_1(natimOrd,respBin),2),'k')
vline(0)
vline(proc.window(2))
xticks([0 proc.window(2) xVal(1) xVal(end)])
xticklabels({'0' num2str(proc.window(2)) '1' num2str(nIm)})
xlabel('Time (s) / Images')
ylabel('Firing rate (sp/s)')
makepretty
offsetAxes

% NatImCorr 
p(1) = subplot(4,4,7); % non-match
hold all
resp3 = nanmean(spikeData{recSesUni == sess1}(:,:,clusterIDs{recSesUni == sess1} == clu3,:),4);
% plot timecourse
plot(bins, mean(resp1_1),'k')
plot(bins, mean(resp3),'color', colNonMatches)
% plot average response
plot(xVal, nanmean(resp1_1(natimOrd,respBin),2),'k')
plot(xVal, nanmean(resp3(natimOrd,respBin),2),'color', colNonMatches)
vline(0)
vline(proc.window(2))
xticks([0 proc.window(2) xVal(1) xVal(end)])
xticklabels({'0' num2str(proc.window(2)) '1' num2str(nIm)})
xlabel('Time (s) / Images')
ylabel('Firing rate (sp/s)')
makepretty
offsetAxes

p(2) = subplot(4,4,11); % non-match
hold all
resp2 = nanmean(spikeData{recSesUni == sess2}(:,:,clusterIDs{recSesUni == sess2} == clu2,:),4);
% plot timecourse
plot(bins, mean(resp1_1),'k')
plot(bins, mean(resp2),'color', colMatches)
% plot average response
plot(xVal, nanmean(resp1_1(natimOrd,respBin),2),'k')
plot(xVal, nanmean(resp2(natimOrd,respBin),2),'color', colMatches)
vline(0)
vline(proc.window(2))
xticks([0 proc.window(2) xVal(1) xVal(end)])
xticklabels({'0' num2str(proc.window(2)) '1' num2str(nIm)})
xlabel('Time (s) / Images')
ylabel('Firing rate (sp/s)')
makepretty
offsetAxes
linkaxes(p,'y')

subplot(4,4,15) % matrix
natImCorrMat = reshape(MatchTable.natImRespCorr, nclus, nclus);
natImCorrMat(natImCorrMat<0.6) = 0.6;
imagesc(natImCorrMat(sortCluIdx,sortCluIdx))
hold on
colormap('RedBlue')
clim([0.35 0.85])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% XCorr
p(1) = subplot(4,4,4); % match
hold all
clu1Idx = find(cluList1 == clu1);
corrVec1_1 = corr(sr{sess1}(clu1Idx,1:floor(size(sr{sess1},2)/2))', sr{sess1}(refPopIdx{sess1},1:floor(size(sr{sess1},2)/2))');
corrVec1_2 = corr(sr{sess1}(clu1Idx,floor(size(sr{sess1},2)/2):end)', sr{sess1}(refPopIdx{sess1},floor(size(sr{sess1},2)/2):end)');
[c, neurOrd] = sort(corrVec1_1,'descend');
neurOrd(c > 0.99999) = [];
plot(corrVec1_2(neurOrd),'color',[0.5 0.5 0.5])
plot(corrVec1_1(neurOrd),'k')
xlabel('Ref. unit')
ylabel('Correlation')
xticks([1 numel(neurOrd)])
makepretty
offsetAxes

p(1) = subplot(4,4,8); % match
hold all
clu3Idx = find(cluList1 == clu3);
corrVec3 = corr(sr{sess1}(clu3Idx,:)', sr{sess1}(refPopIdx{sess1},:)');
corrVec3(corrVec3 > 0.99999) = nan; % same unit
plot(corrVec1_1(neurOrd),'k')
plot(corrVec3(neurOrd),'color',colNonMatches)
xlabel('Ref. unit')
ylabel('Correlation')
xticks([1 numel(neurOrd)])
makepretty
offsetAxes

p(2) = subplot(4,4,12); % non-match
hold all
clu2Idx = find(cluList2 == clu2);
corrVec2 = corr(sr{sess2}(clu2Idx,:)', sr{sess2}(refPopIdx{sess2},:)');
corrVec2(corrVec2 > 0.99999) = nan; % same unit
plot(corrVec1_1(neurOrd),'k')
plot(corrVec2(neurOrd),'color',colMatches)
xlabel('Ref. unit')
ylabel('Correlation')
xticks([1 numel(neurOrd)])
makepretty
offsetAxes
linkaxes(p,'y')

subplot(4,4,16) % matrix
xcorrCorrMat = reshape(MatchTable.refPopCorr, nclus, nclus);
xcorrCorrMat(xcorrCorrMat<0.85) = 0.85;
imagesc(xcorrCorrMat(sortCluIdx,sortCluIdx))
hold on
colormap('RedBlue')
clim([0.75 0.95])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% %% Plot  
% 
% % Get pairs
% sess1 = 1;
% sess2 = 2;
% MatchTable_pair = MatchTable(MatchTable.UID1 == MatchTable.UID2 & MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :);
% MatchTable_bestPairs = MatchTable_pair(MatchTable_pair.ACGCorr > 0.6 & MatchTable_pair.refPopCorr > 0.9 & MatchTable_pair.natImRespCorr > 0.85,:);
% % pair2plt = 7;
% pair2plt = 1;
% clu1 = MatchTable_bestPairs(pair2plt,:).ID1;
% clu2 = MatchTable_bestPairs(pair2plt,:).ID2;
% 
% % Get depth for each pair (used for sorting)
% cluList1 = unique(MatchTable(MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :).ID1);
% cluDepthSess1 = nan(1,numel(cluList1));
% for cluIdx = 1:numel(cluList1)
%    cluDepthSess1(cluIdx) = nanmean(sp.spikeDepths(sp.RecSes == sess1 & sp.spikeTemplates == cluList1(cluIdx)));
% end
% cluList2 = unique(MatchTable(MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :).ID2);
% cluDepthSess2 = nan(1,numel(cluList2));
% for cluIdx = 1:numel(cluList2)
%    cluDepthSess2(cluIdx) = nanmean(sp.spikeDepths(sp.RecSes == sess2 & sp.spikeTemplates == cluList2(cluIdx)));
% end
% 
% % Sorting neurons by depths for matrices
% [~,sortIdx1] = sort(cluDepthSess1,'descend');
% [~,sortIdx2] = sort(cluDepthSess2,'descend');
% sortCluIdx = [sortIdx1 numel(sortIdx1)+sortIdx2];
% 
% % Plot figure
% figure('Position', [400 270 800 500]);
% 
% % FR
% subplot(4,4,9) % matrix
% FRDiffMat = reshape(MatchTable.FRDiff, nclus, nclus);
% imagesc(FRDiffMat(sortCluIdx,sortCluIdx))
% hold on
% colormap(flipud(gray))
% clim([0 30])
% xlabel('Unit_i')
% ylabel('Unit_j')
% hold on
% arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% axis square
% freezeColors
% makepretty
% 
% % ACG
% p(1) = subplot(4,4,2); % day 1
% idx1 = sp.spikeTemplates == clu1 & sp.RecSes == sess1;
% [CCGClu1, tClu1] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
%     ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
% plot(tClu1(tClu1>0), CCGClu1(tClu1>0,1),'k');
% xticks([0 UMparam.ACGduration/2])
% yticks([0 20])
% xlabel('Time (s)')
% ylabel('Firing rate (sp/s)')
% set(gca,'XScale','log')
% makepretty
% offsetAxes
% 
% p(2) = subplot(4,4,6); % day 2
% idx2 = sp.spikeTemplates == clu2 & sp.RecSes == sess2;
% [CCGClu2, tClu2] = CCGBz([double(sp.st(idx2)); double(sp.st(idx2))], [ones(size(sp.st(idx2), 1), 1); ...
%     ones(size(sp.st(idx2), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
% plot(tClu2(tClu2>0), CCGClu2(tClu2>0,1),'r');
% xticks([0 UMparam.ACGduration/2])
% yticks([0 20])
% xlabel('Time (s)')
% ylabel('Firing rate (sp/s)')
% set(gca,'XScale','log')
% linkaxes(p,'xy')
% makepretty
% offsetAxes
% 
% subplot(4,4,10) % matrix
% ACGCorrMat = reshape(MatchTable.ACGCorr, nclus, nclus);
% imagesc(ACGCorrMat(sortCluIdx,sortCluIdx))
% hold on
% colormap('RedBlue')
% clim([-1 1])
% xlabel('Unit_i')
% ylabel('Unit_j')
% hold on
% arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% axis square
% freezeColors
% makepretty
% 
% % XCorr
% subplot(4,4,3); % day 
% hold all
% clu1Idx = find(cluList1 == clu1);
% corrVec1 = corr(sr{sess1}(clu1Idx,:)', sr{sess1}(refPopIdx{sess1},:)');
% [c, neurOrd] = sort(corrVec1,'descend');
% neurOrd(c > 0.99999) = [];
% pltPeriod = 1:10e3;
% pltPeriodTime = pltPeriod*UMparam.binsz;
% mat2plt = cat(1,sr{sess1}(clu1Idx,pltPeriod),zeros(3,numel(pltPeriod)),sr{sess1}(refPopIdx{sess1}(neurOrd),pltPeriod));
% imagesc( 1:size(mat2plt,1), pltPeriodTime, mat2plt')
% set(gca,'YDir','reverse')
% xticks([1 5 size(mat2plt,1)])
% xticklabels({'unit', '1', num2str(numel(refPopIdx{sess1}))})
% colormap(flipud(gray))
% clim([0 2])
% plot(5+(1:numel(neurOrd)), pltPeriodTime(end)+(pltPeriodTime(end)-pltPeriodTime(1))*(max(corrVec1(neurOrd))-corrVec1(neurOrd))/max(corrVec1(neurOrd)),'k')
% axis tight
% xlabel('Ref. unit')
% ylabel('Time (s)')
% freezeColors
% makepretty
% offsetAxes
% 
% subplot(4,4,7) % day 2
% hold all
% clu2Idx = find(cluList2 == clu2);
% corrVec2 = corr(sr{sess2}(clu2Idx,:)', sr{sess2}(refPopIdx{sess2},:)');
% mat2plt = cat(1,sr{sess2}(clu2Idx,pltPeriod),zeros(3,numel(pltPeriod)),sr{sess2}(refPopIdx{sess2}(neurOrd),pltPeriod));
% imagesc( 1:size(mat2plt,1), pltPeriodTime, mat2plt')
% set(gca,'YDir','reverse')
% colormap(flipud(gray))
% clim([0 2])
% plot(5+(1:numel(neurOrd)), pltPeriodTime(end)+(pltPeriodTime(end)-pltPeriodTime(1))*(max(corrVec2(neurOrd))-corrVec2(neurOrd))/max(corrVec2(neurOrd)),'r')
% axis tight
% xlabel('Ref. unit')
% ylabel('Time (s)')
% freezeColors
% makepretty
% offsetAxes
% 
% subplot(4,4,11) % matrix
% xcorrCorrMat = reshape(MatchTable.refPopCorr, nclus, nclus);
% xcorrCorrMat(xcorrCorrMat<.8) = .8;
% imagesc(xcorrCorrMat(sortCluIdx,sortCluIdx))
% hold on
% colormap('RedBlue')
% clim([.6 1])
% xlabel('Unit_i')
% ylabel('Unit_j')
% hold on
% arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% axis square
% freezeColors
% makepretty
% 
% % NatImCorr 
% subplot(4,4,4) % day 1
% hold all
% respBin = bins > 0 & bins < 0.4;
% resp1 = nanmean(spikeData{recSesUni == sess1}(:,:,clusterIDs{recSesUni == sess1} == clu1,:),4);
% [~,natimOrd] = sort(nanmean(resp1(:,respBin),2),'descend');
% imagesc(1:size(resp1,1), bins, resp1(natimOrd,:)');
% set(gca,'YDir','reverse')
% hline(0)
% hline(1.0)
% colormap(flipud(gray))
% clim([0 40])
% xticks([1 112])
% ylabel('Images')
% ylabel('Time (s)')
% tmp = nanmean(resp1(natimOrd,respBin),2);
% plot(bins(end)+(bins(end)-bins(1))*(max(tmp)-tmp)/max(tmp),'k')
% axis tight
% freezeColors
% makepretty
% offsetAxes
% 
% subplot(4,4,8) % day 2
% hold all
% resp2 = nanmean(spikeData{recSesUni == sess2}(:,:,clusterIDs{recSesUni == sess2} == clu2,:),4);
% imagesc(1:size(resp2,1), bins, resp2(natimOrd,:)');
% set(gca,'YDir','reverse')
% hline(0)
% hline(1.0)
% colormap(flipud(gray))
% clim([0 40])
% xticks([1 112])
% ylabel('Images')
% ylabel('Time (s)')
% tmp = nanmean(resp2(natimOrd,respBin),2);
% plot(bins(end)+(bins(end)-bins(1))*(max(tmp)-tmp)/max(tmp),'r')
% axis tight
% freezeColors
% makepretty
% offsetAxes
% 
% subplot(4,4,12) % matrix
% natImCorrMat = reshape(MatchTable.natImRespCorr, nclus, nclus);
% imagesc(natImCorrMat(sortCluIdx,sortCluIdx))
% hold on
% colormap('RedBlue')
% clim([-1 1])
% xlabel('Unit_i')
% ylabel('Unit_j')
% hold on
% arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
% axis square
% freezeColors
% makepretty