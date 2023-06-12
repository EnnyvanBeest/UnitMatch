for midx = 1:length(MiceOpt)
% compare
tmpdir = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
tmpMatchTbl = matfile(fullfile(tmpdir.folder,tmpdir.name));
MatchTable = tmpMatchTbl.MatchTable;


% Now load manual curation
tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','BlindTable.mat'));
tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','manualCuration*.mat'));
if isempty(manualscore)
    continue
end
for id = 1%:length(manualscore)

    tmpman = load(fullfile(manualscore(id).folder,manualscore(id).name));
    if id==1
        Manual = tmpman.match;
    else
        Manual = cat(1,Manual,tmpman.match);
    end

end

% Pairs of Interest
tbl = tmpmanual.tbl;


% UM results
RowIdx = cell2mat(arrayfun(@(X) find(MatchTable.ID1 == tbl.ClusID1(X) & MatchTable.ID2 == tbl.ClusID2(X) & ...
    MatchTable.RecSes1 == tbl.RecID1(X) & MatchTable.RecSes2 == tbl.RecID2(X)),1:height(tbl),'Uni',0))

MatchProb = MatchTable.MatchProb(RowIdx);
PyKS = MatchTable.ID1(RowIdx) == MatchTable.ID2(RowIdx);

% Now Compare
figure('name',['Comparing Methods ' MiceOpt{midx}])
subplot(1,3,1)
scatter(nanmean(Manual,1),MatchProb,10,[0 0 0],'filled')
ylabel('MatchProbability')
xlabel('ManualCuration')
xlim([-1.5 1.5])
makepretty


subplot(1,3,2)
scatter(PyKS,MatchProb,10,[0 0 0],'filled')
ylabel('MatchProb')
xlabel('PyKS')
xlim([-0.5 1.5])
makepretty

subplot(1,3,3)
scatter(PyKS,nanmean(Manual,1),10,[0 0 0],'filled')
xlabel('PyKS')
ylabel('Manual')
xlim([-0.5 1.5])
makepretty

AvgMan = nanmean(Manual,1);

% If PyKS said match
sum(AvgMan(PyKS==1)>0)./sum(PyKS==1)
sum(AvgMan(PyKS==0)>0)./sum(PyKS==0)
sum(AvgMan(PyKS==0)<0)./sum(PyKS==0)



% If UM said Match
sum(AvgMan(MatchProb>0.5)>0)./sum(MatchProb>0.5)
sum(AvgMan(MatchProb<0.5)>0)./sum(MatchProb<0.5)
sum(AvgMan(MatchProb<0.5)<0)./sum(MatchProb<0.5)

end