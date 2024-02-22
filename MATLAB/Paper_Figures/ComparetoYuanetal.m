if 1
%% Within day comparison
datapath = 'H:\MatchingUnits\Yuan\WithinDay';
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.5;
FPFNYuan = nan(length(miceopt),2); %FP/FN
FPFNUM = nan(length(miceopt),2); %FP/FN
AUCs = nan(length(miceopt),2); % UM/Yuan
% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));

    ResultTable = tmpYuan.output.results_wth(:,[3,2]); % This is with their 10 micron cut-off
    % ResultTable = tmpYuan.output.all_results_post(:,[2,3]);

    % Proportion - only for within day
    FPFNYuan(midx,2) = (tmpYuan.output.KSgood_f1 - sum(ResultTable(:,1)==ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were no found?
    FPFNYuan(midx,1) = (sum(ResultTable(:,1)~=ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were found extra?
 
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapath,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    % AssignUniqueID(UMOutputFile.folder)
    % ComputeFunctionalScores(UMOutputFile.folder)

    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    if nclus ~= tmpYuan.output.KSgood_f1
        disp('Not same amount of units.. Unfair comparison??')
        keyboard
    end
   
    % False negatives
    MatchProbability = reshape(tmpUM.MatchTable.MatchProb,nclus,nclus);
    AverageMP = nanmean(cat(3,MatchProbability,MatchProbability'),3);
    FPFNUM(midx,2) = (nclus - sum(AverageMP(logical(eye(nclus)))>ProbThresh))./nclus; % False negative

    % False positives
    tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.ID1 <= tmpUM.MatchTable.ID2)); % Only 1 cv
    Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)]; 
    FPFNUM(midx,1) = sum(Pairs(:,1)~=Pairs(:,2))./nclus;

    % All refPopCorr
    CV1idx = find(tmpUM.MatchTable.ID1 < tmpUM.MatchTable.ID2 & tmpUM.MatchTable.EucledianDistance < tmpUM.UMparam.NeighbourDist); % Only neighbours
    % Extract functional scores
    FunctionalScoreUM = tmpUM.MatchTable.refPopCorr(tblidx);
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx)));

    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreUM'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreUM))];
    [x,y,~,AUCs(midx,1)] = perfcurve(labels,scores,1);
    
    % For Yuan et al.
    ClusIDs = tmpUM.UniqueIDConversion.OriginalClusID;
    RealIDTable = ClusIDs(ResultTable);

    % Numbers
    OverlapIdx = ismember(Pairs,RealIDTable,'rows');
    nOverlap = sum(OverlapIdx);
    NonOverlapIdx = ~ismember(Pairs,RealIDTable,'rows');
    NonOverlapIdxYuan = ~ismember(RealIDTable,Pairs,'rows');
    nUniqueUM = size(Pairs,1)-nOverlap;
    nUniqueYuan = size(RealIDTable,1)-nOverlap;
    
    % Functional scores Yuan
    YuanTblIdx = nan(1,size(RealIDTable,1));
    for pairid = 1:size(ResultTable,1)
        YuanTblIdx(pairid) = find(ismember(tmpUM.MatchTable.ID1,RealIDTable(pairid,1)) & ...
            ismember(tmpUM.MatchTable.ID2,RealIDTable(pairid,2)));
    end
    FunctionalScoreYuan = tmpUM.MatchTable.refPopCorr(YuanTblIdx);
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,YuanTblIdx)));
    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreYuan'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreYuan))];
    [x,y,~,AUCs(midx,2)] = perfcurve(labels,scores,1);

end

%% Histogram
cols = hsv(length(miceopt));
figure; 
subplot(3,2,1)
h = bar(squeeze(nanmean(cat(3,FPFNUM,FPFNYuan),1)));
hold on
for hid = 1:length(h)
    for modeid = 1:2
        if hid == 2
        scatter(repmat(h(hid).XEndPoints(modeid),1,length(miceopt))+(rand(1,length(miceopt))-0.5).*0.1,FPFNYuan(:,modeid)',20,cols,'filled')
        else
        scatter(repmat(h(hid).XEndPoints(modeid),1,length(miceopt))+(rand(1,length(miceopt))-0.5).*0.1,FPFNUM(:,modeid)',20,cols,'filled')
        end
    end
end
set(gca,'XTickLabel',{'Unexpected matches','Unexpected non-matches'})
ylabel('Proportion (of total neurons)')
legend({'UnitMatch','Yuan et al.'})
title('Within-day comparison')
makepretty
offsetAxes

subplot(3,2,2)
hold on
for modeid = 1:2
    line([modeid-0.5 modeid+0.5],[nanmean(AUCs(:,modeid),1) nanmean(AUCs(:,modeid),1)],'color',h(modeid).FaceColor);
    scatter(repmat(modeid,1,length(miceopt))+(rand(1,length(miceopt))-0.5).*0.1,AUCs(:,modeid),20,cols,'filled')
end
set(gca,'XTick',1:2,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('AUCvalues')
title('Within-day AUC values - cross-correlation')
makepretty
offsetAxes


end
%% Across days
datapath = 'H:\MatchingUnits\Yuan\Across2Days';
datapathUM = 'H:\MatchingUnits\Output'
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.5;

nMatches = nan(length(miceopt),3); % Both (overlap), UnitMatch only, Yuan only
PercOverlapWithFunctional = nan(length(miceopt),3); % Shared, Unique to UM, Unique to Yuan
AUCsAcross = nan(length(miceopt),3); % AUC values %Shared, Unique to UM, Unique to Yuan
DriftUM = nan(2,length(miceopt));
DriftYuan = nan(1,length(miceopt));
% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapathUM,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    % AssignUniqueID(UMOutputFile.folder)
    ComputeFunctionalScores(UMOutputFile.folder)

    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    DriftUM(:,midx) = nansum(tmpUM.UMparam.drift(:,2:3,1),1); % total Z-drift

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.RecSes1 < tmpUM.MatchTable.RecSes2));
    Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)];

    % Yuan et al?
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));
    DriftYuan(midx) = -tmpYuan.output.z_mode;
    Ntotal = tmpYuan.output.KSgood_f1 + tmpYuan.output.KSgood_f2;
    if nclus ~= Ntotal
        disp('Not same amount of units.. Unfair comparison??')
    end

    ResultTable = tmpYuan.output.results_wth(:,[3,2]); % This is an index, 1`st session in column 3, 2nd in column 2

    RealIDTable = nan(size(ResultTable));
    for sesid = 1:2
        ClusIDs = tmpUM.UniqueIDConversion.OriginalClusID(tmpUM.UniqueIDConversion.recsesAll'==sesid);
        RealIDTable(:,sesid) = ClusIDs(ResultTable(:,sesid));
    end

    OverlapIdx = ismember(Pairs,RealIDTable,'rows');
    NonOverlapIdx = ~ismember(Pairs,RealIDTable,'rows');
    nOverlap = sum(OverlapIdx);
    NonOverlapIdxYuan = ~ismember(RealIDTable,Pairs,'rows');
    nUniqueUM = size(Pairs,1)-nOverlap;
    nUniqueYuan = size(RealIDTable,1)-nOverlap;

    nMatches(midx,1) = nOverlap/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
    nMatches(midx,2) = nUniqueUM/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
    nMatches(midx,3) = nUniqueYuan/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);

    % All refPopCorr
    CV1idx = find(tmpUM.MatchTable.RecSes1 < tmpUM.MatchTable.RecSes2 & tmpUM.MatchTable.EucledianDistance < tmpUM.UMparam.maxdist); % Only neighbours
    % Extract functional scores
    FunctionalScoreOverlap = tmpUM.MatchTable.refPopCorr(tblidx(OverlapIdx));
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx(OverlapIdx))));

    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreOverlap'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreOverlap))];
    [x,y,~,AUCsAcross(midx,1)] = perfcurve(labels,scores,1);
    
    %  UM 
    FunctionalScoreUM = tmpUM.MatchTable.refPopCorr(tblidx);
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx)));

    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreUM'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreUM))];
    [x,y,~,AUCsAcross(midx,2)] = perfcurve(labels,scores,1);

    % Functional scores Yuan
    YuanTblIdx = nan(1,size(RealIDTable,1));
    for pairid = 1:size(ResultTable,1)
        YuanTblIdx(pairid) = find(ismember(tmpUM.MatchTable.ID1,RealIDTable(pairid,1)) & ismember(tmpUM.MatchTable.RecSes1,1) &...
            ismember(tmpUM.MatchTable.ID2,RealIDTable(pairid,2)) & ismember(tmpUM.MatchTable.RecSes2,2));
    end
    FunctionalScoreYuan = tmpUM.MatchTable.refPopCorr(YuanTblIdx);
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,YuanTblIdx)));
    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreYuan'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreYuan))];
    [x,y,~,AUCsAcross(midx,3)] = perfcurve(labels,scores,1);
      
end

%%
subplot(3,2,3)
CondCols = [0 0 0; 0 0 1; 1 0 0];
hold on
for modeid = 1:3
    line([modeid-0.5 modeid+0.5],[nanmean(nMatches(:,modeid),1) nanmean(nMatches(:,modeid),1)],'color',CondCols(modeid,:));
    scatter(repmat(modeid,1,length(miceopt))+(rand(1,length(miceopt))-0.5).*0.1,nMatches(:,modeid),20,cols,'filled')
end
set(gca,'XTick',1:3,'XTickLabel',{'Overlap','UnitMatch','Yuan et al.'})
ylabel('nMatches')
title('Across day matches')
makepretty
offsetAxes

subplot(3,2,4)
hold on
for modeid = 1:3
    line([modeid-0.5 modeid+0.5],[nanmean(AUCsAcross(:,modeid),1) nanmean(AUCsAcross(:,modeid),1)],'color',CondCols(modeid,:));
    h=scatter(repmat(modeid,1,length(miceopt))+(rand(1,length(miceopt))-0.5).*0.1,AUCsAcross(:,modeid),20,cols,'filled')
end
set(gca,'XTick',1:3,'XTickLabel',{'Overlap','UnitMatch','Yuan et al.'})
ylabel('AUC values')
title('Across day AUC values - cross-correlation')
makepretty
offsetAxes

subplot(3,2,5)
scatter(abs(DriftUM(2,:)),abs(DriftYuan),20,cols,'filled')
axis square
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
lims = [min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])];
set(gca,'xlim',lims,'ylim',lims)

xlabel('UnitMatch')
ylabel('Yuan et al.')
title('Drift Estimate')

%% Across many days
keyboard
datapath = 'H:\MatchingUnits\Yuan\AcrossMultipleDays';
datapathUM = 'H:\MatchingUnits\OutputMonthApart'
miceopt = {'AL032'};
ProbThresh = 0.5;

nMatches = nan(length(miceopt),3); % Both (overlap), UnitMatch only, Yuan only
PercOverlapWithFunctional = nan(length(miceopt),3); % Shared, Unique to UM, Unique to Yuan
AUCsAcross = nan(length(miceopt),3); % AUC values %Shared, Unique to UM, Unique to Yuan
% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapathUM,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    % AssignUniqueID(UMOutputFile.folder)
    % ComputeFunctionalScores(UMOutputFile.folder)

    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));
    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.RecSes1 < tmpUM.MatchTable.RecSes2));
    Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)];

    % Yuan et al?
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','chain_summary.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));
    DriftYuan(midx) = -tmpYuan.output.z_mode;
    Ntotal = tmpYuan.output.KSgood_f1 + tmpYuan.output.KSgood_f2;
    if nclus ~= Ntotal
        disp('Not same amount of units.. Unfair comparison??')
    end

    ResultTable = tmpYuan.output.results_wth(:,[3,2]); % This is an index, 1`st session in column 3, 2nd in column 2

    RealIDTable = nan(size(ResultTable));
    for sesid = 1:2
        ClusIDs = tmpUM.UniqueIDConversion.OriginalClusID(tmpUM.UniqueIDConversion.recsesAll'==sesid);
        RealIDTable(:,sesid) = ClusIDs(ResultTable(:,sesid));
    end

    OverlapIdx = ismember(Pairs,RealIDTable,'rows');
    NonOverlapIdx = ~ismember(Pairs,RealIDTable,'rows');
    nOverlap = sum(OverlapIdx);
    NonOverlapIdxYuan = ~ismember(RealIDTable,Pairs,'rows');
    nUniqueUM = size(Pairs,1)-nOverlap;
    nUniqueYuan = size(RealIDTable,1)-nOverlap;

    nMatches(midx,1) = nOverlap/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
    nMatches(midx,2) = nUniqueUM/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
    nMatches(midx,3) = nUniqueYuan/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);

    % All refPopCorr
    CV1idx = find(tmpUM.MatchTable.RecSes1 < tmpUM.MatchTable.RecSes2 & tmpUM.MatchTable.EucledianDistance < tmpUM.UMparam.maxdist); % Only neighbours
    % Extract functional scores
    FunctionalScoreOverlap = tmpUM.MatchTable.refPopCorr(tblidx(OverlapIdx));
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx(OverlapIdx))));

    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreOverlap'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreOverlap))];
    [x,y,~,AUCsAcross(midx,1)] = perfcurve(labels,scores,1);
    
    %  UM Only
    FunctionalScoreUM = tmpUM.MatchTable.refPopCorr(tblidx(NonOverlapIdx));
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx(NonOverlapIdx))));

    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreUM'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreUM))];
    [x,y,~,AUCsAcross(midx,2)] = perfcurve(labels,scores,1);

    % Functional scores Yuan
    YuanTblIdx = nan(1,size(RealIDTable,1));
    for pairid = 1:size(ResultTable,1)
        YuanTblIdx(pairid) = find(ismember(tmpUM.MatchTable.ID1,RealIDTable(pairid,1)) & ismember(tmpUM.MatchTable.RecSes1,1) &...
            ismember(tmpUM.MatchTable.ID2,RealIDTable(pairid,2)) & ismember(tmpUM.MatchTable.RecSes2,2));
    end
    FunctionalScoreYuan = tmpUM.MatchTable.refPopCorr(YuanTblIdx(NonOverlapIdxYuan));
    FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,YuanTblIdx(NonOverlapIdxYuan))));
    % AUC
    scores = [FunctionalScoreOther', FunctionalScoreYuan'];
    labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreYuan))];
    [x,y,~,AUCsAcross(midx,3)] = perfcurve(labels,scores,1);
      
end
