function DrawPairsUnitMatch(SaveDir, DrawBlind, loadMATsToSave)
if nargin < 2
    DrawBlind = 0;
end
Redo = 0;
if nargin < 3 || isempty(loadMATsToSave)
    loadMATsToSave = 1;
end
if ~loadMATsToSave
    TmpFile = matfile(fullfile(SaveDir, 'UnitMatch.mat')); % Access saved file
else
    TmpFile = load(fullfile(SaveDir, 'UnitMatch.mat'));
end

UMparam = TmpFile.UMparam; % Extract parameters
MatchTable = TmpFile.MatchTable; % Load Matchtable
try
    AllSessionCorrelations = TmpFile.AllSessionCorrelations;
catch ME
    disp('Run ComputeFunctionalScores first')
    return
end
WaveformInfo = TmpFile.WaveformInfo;

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
% Pairs redefined:
uId = unique(UniqueID);
Pairs = arrayfun(@(X) find(UniqueID == X), uId, 'Uni', 0);
Pairs(cellfun(@length, Pairs) == 1) = [];

MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
lowselfscores = find(diag(MatchProb) < UMparam.ProbabilityThreshold);
for id = 1:length(lowselfscores) % Add these for plotting - inspection
    Pairs{end+1} = [lowselfscores(id), lowselfscores(id)];
end

RankSc = reshape(MatchTable.refPopSig,nclus,nclus);
Rank = reshape(MatchTable.refPopRank,nclus,nclus);
[HighRankLowProbr HighRankLowProbc] = find(RankSc==1 & Rank == 1 & MatchProb<0.5);
PairsRank = [HighRankLowProbr,HighRankLowProbc];
for id = 1:size(PairsRank,1)
    Pairs{end+1} = PairsRank(id,:);
end

PyKSLabel = false(nclus, nclus);
label = MatchProb > UMparam.ProbabilityThreshold;

if UMparam.RunPyKSChronicStitched
    PairsPyKS = [];
    for uid = 1:nclus
        pairstmp = find(OriID == OriID(uid))';
        if length(pairstmp) > 1
            PairsPyKS = cat(1, PairsPyKS, pairstmp);
        end
    end

    for pid = 1:size(PairsPyKS, 1)
        PyKSLabel(PairsPyKS(pid, 1), PairsPyKS(pid, 2)) = true;
        PyKSLabel(PairsPyKS(pid, 2), PairsPyKS(pid, 1)) = true;
    end
    PyKSLabel(logical(eye(size(PyKSLabel)))) = true;
    [r, c] = find(label == 0 & PyKSLabel == 1);
    PairsPKS = cat(2, r, c);
    PairsPKS(r == c, :) = [];
    PairsPKS = sort(PairsPKS, 2, 'ascend');
    PairsPKS = unique(PairsPKS, 'stable', 'rows');
    Pairs = cat(2, Pairs, arrayfun(@(X) PairsPKS(X, :), 1:length(PairsPKS), 'Uni', 0)); % Add these for plotting - inspection
end

if DrawBlind
    % Keep length 2 for each pair
    tmptbl = dir(fullfile(UMparam.SaveDir, 'BlindFigures', 'BlindTable.mat'));

    if ~Redo && ~isempty(tmptbl)
        tmptbl = load(fullfile(tmptbl(1).folder, tmptbl(1).name));
        Pair1 = cell2mat(arrayfun(@(X) find(ismember(OriID, tmptbl.tbl.ClusID1(X)) & ismember(recses, tmptbl.tbl.RecID1(X))), 1:height(tmptbl.tbl), 'Uni', 0));
        Pair2 = cell2mat(arrayfun(@(X) find(ismember(OriID, tmptbl.tbl.ClusID2(X)) & ismember(recses, tmptbl.tbl.RecID2(X))), 1:height(tmptbl.tbl), 'Uni', 0));
        Pairs = arrayfun(@(X) [Pair1(X), Pair2(X)], 1:length(Pair1), 'Uni', 0);
    else
        Pairs = cellfun(@(X) X(1:2), Pairs, 'Uni', 0);
        % Need to add some low probabilities, but around same location
        EucledianDistance = reshape(MatchTable.EucledianDistance, nclus, nclus);
        [r, c] = find(label == 0 & PyKSLabel == 0 & EucledianDistance < UMparam.NeighbourDist);
        LowProb = cat(2, r, c);
        LowProb(r == c, :) = [];
        LowProb = sort(LowProb, 2, 'ascend');
        LowProb = unique(LowProb, 'stable', 'rows');
        % Match number of pairs for equal set
        npairs = length(Pairs);
        LowProb = LowProb(randsample(length(LowProb), npairs, 0), :);
        Pairs = cat(2, Pairs, arrayfun(@(X) LowProb(X, :), 1:length(LowProb), 'Uni', 0));

        % add some diagonal (i=j within)
        addIisJ = randsample(nclus, round(0.1*npairs), 0);
        addIisJ = repmat(addIisJ, 1, 2);
        Pairs = cat(2, Pairs, arrayfun(@(X) addIisJ(X, :), 1:length(addIisJ), 'Uni', 0));
        % Randomly shuffle Pairs so we don't have matches grouped together and
        % different order for different runs
        Pairs = Pairs(randsample(length(Pairs), length(Pairs), 0));

    end

    if size(Pairs, 2) > UMparam.drawmax
        DrawPairs = randsample(1:size(Pairs, 2), UMparam.drawmax, 'false');
    else
        DrawPairs = 1:size(Pairs, 2);
    end
    PlotTheseUnits_UM_Blind(Pairs(DrawPairs), MatchTable, UniqueIDConversion, WaveformInfo, AllSessionCorrelations, UMparam)
else
    if size(Pairs, 2) > UMparam.drawmax
        DrawPairs = randsample(1:size(Pairs, 2), UMparam.drawmax, 'false');
    else
        DrawPairs = 1:size(Pairs, 2);
    end
    PlotTheseUnits_UM(Pairs(DrawPairs), MatchTable, UniqueIDConversion, WaveformInfo, AllSessionCorrelations, UMparam)
end


return