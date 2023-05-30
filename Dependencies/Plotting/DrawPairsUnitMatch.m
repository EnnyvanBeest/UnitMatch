function DrawPairsUnitMatch(SaveDir)

TmpFile = matfile(fullfile(SaveDir,'UnitMatch.mat')); % Access saved file
UMparam = TmpFile.UMparam; % Extract parameters
MatchTable = TmpFile.MatchTable; % Load Matchtable
AllSessionCorrelations = TmpFile.AllSessionCorrelations;
WaveformInfo = TmpFile.WaveformInfo;

% Extract cluster information
UniqueIDConversion = TmpFile.UniqueIDConversion;
if UMparam.GoodUnitsOnly
    GoodId = logical(UniqueIDConversion.GoodID);
else
    GoodId = true(1,length(UniqueIDConversion.GoodID));
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
Pairs = arrayfun(@(X) find(UniqueID==X),uId,'Uni',0);
Pairs(cellfun(@length,Pairs)==1) = [];

MatchProb = reshape(MatchTable.MatchProb,nclus,nclus);
lowselfscores = find(diag(MatchProb)<UMparam.ProbabilityThreshold);
for id =1:length(lowselfscores) % Add these for plotting - inspection
    Pairs{end+1} = [lowselfscores(id) lowselfscores(id)];
end

if UMparam.RunPyKSChronicStitched
    label = MatchProb>UMparam.ProbabilityThreshold;
    PairsPyKS = [];
    for uid = 1:nclus
        pairstmp = find(OriID==OriID(uid))';
        if length(pairstmp)>1
            PairsPyKS = cat(1,PairsPyKS,pairstmp);
        end
    end

    PyKSLabel = false(nclus,nclus);
    for pid = 1:size(PairsPyKS,1)
        PyKSLabel(PairsPyKS(pid,1),PairsPyKS(pid,2)) = true;
        PyKSLabel(PairsPyKS(pid,2),PairsPyKS(pid,1)) = true;
    end
    PyKSLabel(logical(eye(size(PyKSLabel)))) = true;
    [r,c] = find(label==0 & PyKSLabel==1);
    PairsPKS = cat(2,r,c);
    PairsPKS(r==c,:) = [];
    PairsPKS = sort(PairsPKS,2,'ascend');
    PairsPKS = unique(PairsPKS,'stable','rows');
    for id = 1:size(PairsPKS,1) % Add these for plotting - inspection
        Pairs{end+1} = [PairsPKS(id,1) PairsPKS(id,2)];
    end
end

if size(Pairs,2)>UMparam.drawmax
    DrawPairs = randsample(1:size(Pairs,2),UMparam.drawmax,'false');
else
    DrawPairs = 1:size(Pairs,2);
end

PlotTheseUnits_UM(Pairs(DrawPairs),MatchTable,UniqueIDConversion,WaveformInfo,AllSessionCorrelations,UMparam)



return