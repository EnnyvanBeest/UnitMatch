% UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_new\AL032\19011111882\2\UnitMatch\UnitMatch.mat"};
% UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_new\AL036\19011116882\3\UnitMatch\UnitMatch.mat"};
UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_new\CB016\Probe0\IMRO_1\UnitMatch\UnitMatch.mat"};
load(UMFile{1})

%%
splitUnitsUIDs = unique(MatchTable((MatchTable.UID1 == MatchTable.UID2) ...
    & (MatchTable.ID1 ~= MatchTable.ID2) ...
    & (MatchTable.RecSes1 == MatchTable.RecSes2),:).UID1);
splitUnitsIdx = ismember(MatchTable.UID1,splitUnitsUIDs);

%%
MatchTable.FlipMatchProb = mat2vec(reshape(MatchTable.MatchProb,[sqrt(size(MatchTable,1)), sqrt(size(MatchTable,1))])');
tmp = MatchTable(1:sqrt(size(MatchTable,1))+1:size(MatchTable),:).MatchProb;
MatchTable.MatchProbID1 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))])');
MatchTable.MatchProbID2 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))]));
tmp = MatchTable(1:sqrt(size(MatchTable,1))+1:size(MatchTable),:).ISICorr;
MatchTable.ISICorrID1 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))])');
MatchTable.ISICorrID2 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))]));

%%

FPNames = {'ISICorr','natImRespCorr','refPopCorr'};
probeBins = [0 10.^(-7:0.2:0)]; %[0:0.05:1]; %
probeBins(end) = 1.01; % to include p=1 in last bin
x = .5*probeBins(1:end-1)+.5*probeBins(2:end);
matchIdx = MatchTable.UID1Liberal == MatchTable.UID2Liberal;
figure;
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};
    y = nan(1, numel(probeBins)-1);
    err = nan(1, numel(probeBins)-1);
    n = nan(1, numel(probeBins)-1);
    y_match = nan(1, numel(probeBins)-1);
    err_match = nan(1, numel(probeBins)-1);
    n_match = nan(1, numel(probeBins)-1);
    y_nonmatch = nan(1, numel(probeBins)-1);
    err_nonmatch = nan(1, numel(probeBins)-1);
    n_nonmatch = nan(1, numel(probeBins)-1);
    subplot(1,numel(FPNames),fpIdx); hold all
    for pb = 1:numel(probeBins)-1
        idx = MatchTable.MatchProb >= probeBins(pb) & MatchTable.MatchProb < probeBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10 & ...
            ~splitUnitsIdx;
        % all
        n(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx)));
        y(pb) = nanmean(MatchTable.(FPNameCurr)(idx));
        err(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx))/sqrt(n(pb));
        % matches only
        n_match(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx & matchIdx)));
        y_match(pb) = nanmean(MatchTable.(FPNameCurr)(idx & matchIdx));
        err_match(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx & matchIdx))/sqrt(n(pb));
        % non-matches only
        n_nonmatch(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx & ~matchIdx)));
        y_nonmatch(pb) = nanmean(MatchTable.(FPNameCurr)(idx & ~matchIdx));
        err_nonmatch(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx & ~matchIdx))/sqrt(n(pb));
    end
%     errorbar(x,y,err,'k')
%     errorbar(x,y_match,err_match,'r')
%     errorbar(x,y_nonmatch,err_nonmatch,'b')
%     plot(x,y,'k.-')
    plot(x,y_match,'r.-')
    plot(x,y_nonmatch,'b.-')
    ylim([0 1])
    xlabel('Matching probability')
    ylabel('Correlation')
    title(FPNameCurr)
%     set(gca,'XScale','log')
end

%%
figure; hold all
% histogram(MatchTable(matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')
pb = 1;
idx = MatchTable.MatchProb >= probeBins(pb) & MatchTable.MatchProb < probeBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')
pb = 2;
idx = MatchTable.MatchProb >= probeBins(pb) & MatchTable.MatchProb < probeBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')

%%
figure; hold all
pb = 1;
idx = MatchTable.MatchProb >= probeBins(pb) & MatchTable.MatchProb < probeBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).FlipMatchProb,0:0.03:1,'Normalization','probability')
pb = 2;
idx = MatchTable.MatchProb >= probeBins(pb) & MatchTable.MatchProb < probeBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).FlipMatchProb,0:0.03:1,'Normalization','probability')