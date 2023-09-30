%% Load and format data

UMFile = {'\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\2ConsecutiveDays\Non_Stitched\AL032\AllProbes\AllIMRO\UnitMatch\UnitMatch.mat'};
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

% Param for processing
proc.window = [-0.3 1.5 ... % around onset
    0.0 0.0]; % around offset
proc.binSize = 0.002; % in ms
proc.smoothSize = 5; % PSTH smoothing filter
gw = gausswin(proc.smoothSize,3);
proc.smWin = gw./sum(gw);

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

binsOn = (proc.window(1)+proc.binSize/2):proc.binSize:proc.window(2);
binsOff =  (proc.window(3)+proc.binSize/2):proc.binSize:proc.window(4);
bins = [binsOn binsOff];

%% Get binned FR for xcorr

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
    sr{did} = sr{did}(sortIdx,:);
end

%% Plot  

% Get pairs
sess1 = 1;
sess2 = 2;
MatchTable_pair = MatchTable(MatchTable.UID1 == MatchTable.UID2 & MatchTable.RecSes1 == sess1 & MatchTable.RecSes2 == sess2, :);
MatchTable_bestPairs = MatchTable_pair(MatchTable_pair.ACGRankScore < 2 & MatchTable_pair.RankScore < 2 & MatchTable_pair.NImgRankScore < 2,:);
pair2plt = 7;
clu1 = MatchTable_bestPairs(pair2plt,:).ID1;
clu2 = MatchTable_bestPairs(pair2plt,:).ID2;

% Plot figure
figure('Position', [400 270 800 500]);

% FR
subplot(3,4,1)
% Nothing?

subplot(3,4,9) % matrix
FRDiffMat = reshape(MatchTable.FRDiff, nclus, nclus);
imagesc(FRDiffMat)
hold on
colormap(flipud(gray))
clim([0 30])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
title('NatIm Fingerprint')
axis square
freezeColors
makepretty

% ACG
p(1) = subplot(3,4,2); % day 1
idx1 = sp.spikeTemplates == clu1 & sp.RecSes == sess1;
[CCGClu1, tClu1] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
    ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
plot(tClu1, CCGClu1(:,1),'k');
xticks([-UMparam.ACGduration/2 0 UMparam.ACGduration/2])
yticks([0 20])
xlabel('time (s)')
ylabel('Firing rate (sp/s)')
makepretty
offsetAxes

p(2) = subplot(3,4,6); % day 2
idx2 = sp.spikeTemplates == clu2 & sp.RecSes == sess2;
[CCGClu2, tClu2] = CCGBz([double(sp.st(idx2)); double(sp.st(idx2))], [ones(size(sp.st(idx2), 1), 1); ...
    ones(size(sp.st(idx2), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
plot(tClu2, CCGClu2(:,1),'r');
xticks([-UMparam.ACGduration/2 0 UMparam.ACGduration/2])
yticks([0 20])
xlabel('time (s)')
ylabel('Firing rate (sp/s)')
linkaxes(p,'xy')
makepretty
offsetAxes

subplot(3,4,10) % matrix
ACGCorrMat = reshape(MatchTable.ACGCorr, nclus, nclus);
imagesc(ACGCorrMat)
hold on
colormap('RedBlue')
clim([-1 1])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% XCorr
subplot(3,4,3) % day 1

subplot(3,4,7) % day 2

subplot(3,4,11) % matrix
xcorrCorrMat = reshape(MatchTable.FingerprintCor, nclus, nclus);
imagesc(xcorrCorrMat)
hold on
colormap('RedBlue')
clim([-1 1])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty

% NatImCorr 
subplot(3,4,4) % day 1
resp1 = nanmean(spikeData{recSesUni == sess1}(:,:,clusterIDs{recSesUni == sess1} == clu1,:),4);
[~,natimOrd] = sort(nanmean(resp1(:,bins > 0 & bins < 0.4),2),'descend');
imagesc(bins, 1:size(resp1,1), resp1(natimOrd,:));
hold all
vline(0)
vline(1.0)
colormap(flipud(gray))
clim([0 40])
yticks([1 112])
freezeColors
makepretty

subplot(3,4,8) % day 2
resp2 = nanmean(spikeData{recSesUni == sess2}(:,:,clusterIDs{recSesUni == sess2} == clu2,:),4);
imagesc(bins, 1:size(resp2,1), resp2(natimOrd,:));
hold all
vline(0)
vline(1.0)
colormap(flipud(gray))
clim([0 40])
yticks([1 112])
freezeColors
makepretty

subplot(3,4,12) % matrix
natImCorrMat = reshape(MatchTable.NatImCorr, nclus, nclus);
imagesc(natImCorrMat)
hold on
colormap('RedBlue')
clim([-1 1])
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
axis square
freezeColors
makepretty