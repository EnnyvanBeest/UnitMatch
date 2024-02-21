if 1
%% Within day comparison
datapath = 'H:\MatchingUnits\Yuan\WithinDay';
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.5;

FPFNYuan = nan(length(miceopt),2); %FP/FN
FPFNUM = nan(length(miceopt),2); %FP/FN

% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));

    ResultTable = tmpYuan.output.results_wth(:,[2,3]); % This is with their 10 micron cut-off
    % ResultTable = tmpYuan.output.all_results_post(:,[2,3]);

    % Proportion - only for within day
    FPFNYuan(midx,2) = (tmpYuan.output.KSgood_f1 - sum(ResultTable(:,1)==ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were no found?
    FPFNYuan(midx,1) = (sum(ResultTable(:,1)~=ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were found extra?

    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapath,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    AssignUniqueID(UMOutputFile.folder)
    ComputeFunctionalScores(UMOutputFile.folder)

    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    if nclus ~= tmpYuan.output.KSgood_f1
        disp('Not same amount of units.. Unfair comparison??')
        keyboard
    end
   
    MatchProbability = reshape(tmpUM.MatchTable.MatchProb,nclus,nclus);
    AverageMP = nanmean(cat(3,MatchProbability,MatchProbability'),3);
    FPFNUM(midx,2) = (nclus - sum(AverageMP(logical(eye(nclus)))>ProbThresh))./nclus; % False negative

    tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.RecSes1 == tmpUM.MatchTable.RecSes2));
    Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)];
    FunctionalScoreUM = tmpUM.MatchTable.refPopCorr(tblidx);



    % proportion
    FPFNUM(midx,1) = sum(Pairs(:,1)~=Pairs(:,2))./nclus;

    % Evaluate the false positives
    scores = [WavformMSE(SameIdx(:))', WavformMSE(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'WavformMSE'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

end

% Histogram
figure; barwitherr(squeeze(nanstd(cat(3,FPFNYuan,FPFNUM),[],1))./sqrt(length(miceopt)-1),squeeze(nanmean(cat(3,FPFNYuan,FPFNUM),1)));
set(gca,'XTickLabel',{'False Positives','False Negatives'})
ylabel('Proportion (of total neurons)')
legend({'Yuan et al.',['UnitMatch, p=' num2str(ProbThresh)]})
title('Within-day comparison')

makepretty
offsetAxes
end
%% Across days
datapath = 'H:\MatchingUnits\Yuan\Across2Days';
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.5;

nMatches = nan(length(miceopt),2);
PercOverlapWithFunctional = nan(length(miceopt),3); % Shared, Unique to UM, Unique to Yuan

% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapath,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    AssignUniqueID(UMOutputFile.folder)
    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.RecSes1 < tmpUM.MatchTable.RecSes2));
    Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)];

    nMatches(midx,1) = size(Pairs,1);
    FunctionalRankUM = tmpUM.MatchTable.refPopRank(tblidx);

    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));
    Ntotal = tmpYuan.output.KSgood_f1 + tmpYuan.output.KSgood_f2;
    if nclus ~= Ntotal
        disp('Not same amount of units.. Unfair comparison??')
        keyboard
    end

    ResultTable = tmpYuan.output.results_wth(:,[3,2]); % This is an index
    nMatches(midx,2) = size(ResultTable,1);

    RealIDTable = nan(size(ResultTable));
    for sesid = 1:2
        ClusIDs = tmpUM.UniqueIDConversion.OriginalClusID(tmpUM.UniqueIDConversion.recsesAll'==sesid);
        RealIDTable(:,sesid) = ClusIDs(ResultTable(:,sesid));
    end

    OverlapIdx = sum(ismember(Pairs,RealIDTable),2)==2;
    NonOverlapIdx = sum(ismember(Pairs,RealIDTable),2)<2;
    nOverlap = sum(sum(ismember(Pairs,RealIDTable),2)==2);
    NonOverlapIdxYuan = sum(ismember(RealIDTable,Pairs),2)<2;
    nUniqueUM = size(Pairs,1)-nOverlap;
    nUniqueYuan = size(RealIDTable,1)-nOverlap;
    FunctionalRankYuan = nan(1,size(RealIDTable,1));
    for pairid = 1:size(ResultTable,1)
       Idx = find(ismember(tmpUM.MatchTable.ID1,RealIDTable(pairid,1)) & ismember(tmpUM.MatchTable.RecSes1,1) & ...
             ismember(tmpUM.MatchTable.ID2,RealIDTable(pairid,2)) & ismember(tmpUM.MatchTable.RecSes2,2));
       FunctionalRankYuan(pairid) = tmpUM.MatchTable.refPopRank(Idx);
    end

    %
     PercOverlapWithFunctional(midx,1) = sum(FunctionalRankUM(OverlapIdx) == 1)./sum(OverlapIdx)*100;
     PercOverlapWithFunctional(midx,2) = sum(FunctionalRankUM(NonOverlapIdx) == 1)./sum(NonOverlapIdx)*100;
     PercOverlapWithFunctional(midx,3) = sum(FunctionalRankYuan(NonOverlapIdxYuan)==1)./sum(NonOverlapIdxYuan)*100;
    
end
figure('name',['Overlap with functional scores'])
subplot(1,2,1)
barwitherr(nanstd(PercOverlapWithFunctional,[],1)./sqrt(length(miceopt)-1),nanmean(PercOverlapWithFunctional,1))
set(gca,'XTickLabel',{'Both UM & Yuan','UM Only','Yuan Only'})
ylabel('Percentage overlap')
title('Overlap with Rank == 1')
makepretty
offsetAxes

subplot(1,2,2)
bar(nMatches)
legend({'UnitMatch','Yuan et al'})
title('Total number of matches')
makepretty
offsetAxes
