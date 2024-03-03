if 1
%% Within day comparison
datapath = 'H:\MatchingUnits\Yuan\WithinDay';
datapathUM = 'H:\FigShare_UnitMatch\OutputOneMouse\'
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
    UMOutputFile = dir(fullfile(datapathUM,miceopt{midx},'**','UnitMatch.mat')); % Find output file
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
    nUniqueUM = size(Pairs,1);
    nUniqueYuan = size(RealIDTable,1);
    
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
subplot(3,3,1)
h = bar(squeeze(nanmean(cat(3,FPFNUM,FPFNYuan),1)));
hold on
for hid = 1:length(h)
    for modeid = 1:2
        if hid == 2
        scatter(repmat(h(hid).XEndPoints(modeid),1,length(miceopt)),FPFNYuan(:,modeid)',20,cols,'filled')
        else
        scatter(repmat(h(hid).XEndPoints(modeid),1,length(miceopt)),FPFNUM(:,modeid)',20,cols,'filled')
        end

    end
end
for modeid = 1: 2
    for midx = 1:length(miceopt)
        plot([h(1).XEndPoints(modeid);h(2).XEndPoints(modeid)],cat(2,FPFNUM(midx,modeid),FPFNYuan(midx,modeid))','color',cols(midx,:))
    end
end

set(gca,'XTickLabel',{'Unexpected matches','Unexpected non-matches'})
ylabel('Within day')
legend({'UnitMatch','Yuan et al.'})
title('Proportion of matches')
makepretty
offsetAxes

subplot(3,3,2)
hold on
for modeid = 1:2
    scatter(repmat(modeid,1,length(miceopt)),AUCs(:,modeid),20,cols,'filled')
end
for midx = 1:length(miceopt)
    plot([1,2],AUCs(midx,:),'color',cols(midx,:))
end
set(gca,'XTick',1:2,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('AUCvalues')
title('Ref pop cross-correlation')
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

    DriftUM(:,midx) = nansum(tmpUM.UMparam.drift(1,2:3,:),3); % total Z-drift

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
    nUniqueUM = size(Pairs,1);
    nUniqueYuan = size(RealIDTable,1);

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
subplot(3,3,4)
CondCols = [0 0 0; 0 0 1; 1 0 0];
hold on
for modeid = 2:3
    scatter(repmat(modeid,1,length(miceopt)),nMatches(:,modeid),20,cols,'filled')
end
for midx = 1:length(miceopt)
    plot([2:3],nMatches(midx,2:3),'color',cols(midx,:))
end
set(gca,'XTick',2:3,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('Across two days')
makepretty
offsetAxes

subplot(3,3,5)
hold on
for modeid = 2:3
    h=scatter(repmat(modeid,1,length(miceopt)),AUCsAcross(:,modeid),20,cols,'filled')
end
for midx = 1:length(miceopt)
    plot([2,3],AUCsAcross(midx,2:3),'color',cols(midx,:))
end
set(gca,'XTick',2:3,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('AUC values')
makepretty
offsetAxes

subplot(3,3,6)
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
datapath = 'H:\MatchingUnits\Yuan\AcrossManyDays';
datapathUM = 'H:\UnitMatch\AL032\Probe0\IMRO_2\UnitMatch';% 'H:\MatchingUnits\OutputMonthApart'
miceopt = {'AL032'};
ProbThresh = 0.5;

nMatches = nan(3,3); % Both (overlap), UnitMatch only, Yuan only
AUCsAcrossManyDays = nan(3,3); % AUC values %Shared, Unique to UM, Unique to Yuan
% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapathUM,'**','UnitMatch.mat')); % Find output file
    % AssignUniqueID(UMOutputFile.folder)
    % ComputeFunctionalScores(UMOutputFile.folder)

    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));
    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
   
    % Yuan et al?
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','chain_summary.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));

    for did = 1:3 % Loop over days to find the matches Yuan found
        tmpYuan = load(fullfile(YuanOutputFile.folder,['result_' num2str(did) '_' num2str(did+1)],'Output.mat'));
        if did ==1
            Ntotal = tmpYuan.output.KSgood_f1 + tmpYuan.output.KSgood_f2;
        else
            Ntotal = Ntotal + tmpYuan.output.KSgood_f2;
        end
        ResultTable = tmpYuan.output.results_wth(:,[3,2]); % This is an index, 1`st session in column 3, 2nd in column 2

        RealIDTable = nan(size(ResultTable));
        for sesid = 1:2
            ClusIDs = tmpUM.UniqueIDConversion.OriginalClusID(tmpUM.UniqueIDConversion.recsesAll'==did+sesid-1);
            RealIDTable(:,sesid) = ClusIDs(ResultTable(:,sesid));
        end
        
        YuanTblIdx = nan(1,0); % Save out match table Idx of all identified pairs
        % Functional scores Yuan
        for pairid = 1:size(ResultTable,1)
            YuanTblIdx = cat(2,YuanTblIdx,find(ismember(tmpUM.MatchTable.ID1,RealIDTable(pairid,1)) & ismember(tmpUM.MatchTable.RecSes1,did) &...
                ismember(tmpUM.MatchTable.ID2,RealIDTable(pairid,2)) & ismember(tmpUM.MatchTable.RecSes2,did+1)));
        end

        % UnitMatch for this loop
        tblidx = find((tmpUM.MatchTable.UID1 == tmpUM.MatchTable.UID2) & (tmpUM.MatchTable.RecSes1 == did & tmpUM.MatchTable.RecSes2 == did+1));
        Pairs = [tmpUM.MatchTable.ID1(tblidx) tmpUM.MatchTable.ID2(tblidx)];

        OverlapIdx = ismember(Pairs,RealIDTable,'rows');
        NonOverlapIdx = ~ismember(Pairs,RealIDTable,'rows');
        nOverlap = sum(OverlapIdx);
        NonOverlapIdxYuan = ~ismember(RealIDTable,Pairs,'rows');
        nUniqueUM = size(Pairs,1);
        nUniqueYuan = size(RealIDTable,1);

        nMatches(did,1) = nOverlap./min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
        nMatches(did,2) = nUniqueUM/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);
        nMatches(did,3) = nUniqueYuan/min([tmpYuan.output.KSgood_f1,tmpYuan.output.KSgood_f2]);

        % All refPopCorr
        CV1idx = find(tmpUM.MatchTable.RecSes1 == did & tmpUM.MatchTable.RecSes2 == did+1 & tmpUM.MatchTable.EucledianDistance < tmpUM.UMparam.maxdist); % Only neighbours
        % Extract functional scores
        FunctionalScoreOverlap = tmpUM.MatchTable.refPopCorr(tblidx(OverlapIdx));
        FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx(OverlapIdx))));

        % AUC
        scores = [FunctionalScoreOther', FunctionalScoreOverlap'];
        labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreOverlap))];
        [x,y,~,AUCsAcrossManyDays(did,1)] = perfcurve(labels,scores,1);

        %  UM Only
        FunctionalScoreUM = tmpUM.MatchTable.refPopCorr(tblidx(NonOverlapIdx));
        FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,tblidx(NonOverlapIdx))));

        % AUC
        scores = [FunctionalScoreOther', FunctionalScoreUM'];
        labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreUM))];
        [x,y,~,AUCsAcrossManyDays(did,2)] = perfcurve(labels,scores,1);


        FunctionalScoreYuan = tmpUM.MatchTable.refPopCorr(YuanTblIdx(NonOverlapIdxYuan));
        FunctionalScoreOther = tmpUM.MatchTable.refPopCorr(CV1idx(~ismember(CV1idx,YuanTblIdx(NonOverlapIdxYuan))));
        % AUC
        scores = [FunctionalScoreOther', FunctionalScoreYuan'];
        labels = [zeros(1,length(FunctionalScoreOther)), ones(1,length(FunctionalScoreYuan))];
        [x,y,~,AUCsAcrossManyDays(did,3)] = perfcurve(labels,scores,1);


    end


end


subplot(3,3,7)
cols = [1 0 0; 0 1 0; 0 0 1];
hold on
for modeid = 2:3
    scatter(repmat(modeid,1,3),nMatches(:,modeid),20,cols,'filled')
end
for midx = 1:3
    plot([2,3],nMatches(midx,2:3),'color',cols(midx,:))
end
set(gca,'XTick',2:3,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('Across 4 recordings')
makepretty
offsetAxes

subplot(3,3,8)
hold on
for modeid = 2:3
    h=scatter(repmat(modeid,1,3),AUCsAcrossManyDays(:,modeid)',20,cols(1:3,:),'filled')
end
for midx = 1:3
    plot([2,3],AUCsAcrossManyDays(midx,2:3),'color',cols(midx,:))
end
set(gca,'XTick',2:3,'XTickLabel',{'UnitMatch','Yuan et al.'})
ylabel('AUC values')
makepretty
offsetAxes
%% 