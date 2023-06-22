for midx = 1:length(MiceOpt)
    % compare
    tmpdir = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
    tmpMatchTbl = matfile(fullfile(tmpdir.folder,tmpdir.name));
    MatchTable = tmpMatchTbl.MatchTable;
    nclus = sqrt(height(MatchTable));
    MatchProbability = reshape(MatchTable.MatchProb,nclus,nclus);
    RankScore = reshape(MatchTable.RankScore,nclus,nclus);
    UMparam = tmpMatchTbl.UMparam;
    UniqueIDConversion = tmpMatchTbl.UniqueIDConversion;
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

    for fid = 1:2
        % Now load manual curation
        if fid == 1
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end
            ScorerNames = cell(length(manualscore),2);
            AddName = 'WO';

            % Pairs of Interest
            tbl = tmpmanual.tbl;

            % UM results
            Pair1 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID1(X)) & ismember(recses,tbl.RecID1(X))),1:height(tbl),'Uni',0));
            Pair2 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID2(X)) & ismember(recses,tbl.RecID2(X))),1:height(tbl),'Uni',0));
            Pairs = arrayfun(@(X) [Pair1(X) Pair2(X)],1:length(Pair1),'Uni',0);
            MatchProb = cell2mat(cellfun(@(X) MatchProbability(X(1),X(2)),Pairs,'Uni',0))';
            Rank = cell2mat(cellfun(@(X) RankScore(X(1),X(2)),Pairs,'Uni',0))';

            PyKS = cell2mat(cellfun(@(X) OriID(X(1)) == OriID(X(2)),Pairs,'Uni',0));
            Order =  1:height(tbl);
        else
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            tbl2 = tmpmanual.tbl;
            % Order?
            Order = cell2mat(arrayfun(@(X) find(tbl2.ClusID1 == tbl.ClusID1(X) & tbl2.ClusID2 == tbl.ClusID2(X) & tbl2.RecID1 == tbl.RecID1(X) & tbl2.RecID2 == tbl.RecID2(X),1,'first'),1:height(tbl),'Uni',0));

            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end
            AddName = 'With';

        end
        for id = 1:length(manualscore)

            tmpman = load(fullfile(manualscore(id).folder,manualscore(id).name));
            if id==1 && fid==1
                Manual = tmpman.match(Order);
            else
                Manual = cat(1,Manual,tmpman.match(Order));
            end
            tmpname = strsplit(manualscore(id).name,'manualCuration_');
            tmpname = strsplit(tmpname{2},'.mat');

            ScorerNames{id,fid} = [tmpname{1} AddName];
        end
    end
    ScorerNames = ScorerNames(~cellfun(@isempty,ScorerNames));

    disp('User scoring correlations:')
    tmpcor = corr(Manual')

    %
    %     % Now Compare
    %     figure('name',['Comparing Methods ' MiceOpt{midx}])
    %     subplot(1,3,1)
    %     scatter(nanmean(Manual,1),MatchProb,10,[0 0 0],'filled')
    %     ylabel('MatchProbability')
    %     xlabel('ManualCuration')
    %     xlim([-1.5 1.5])
    %     makepretty
    %
    %
    %     subplot(1,3,2)
    %     scatter(PyKS,MatchProb,10,[0 0 0],'filled')
    %     ylabel('MatchProb')
    %     xlabel('PyKS')
    %     xlim([-0.5 1.5])
    %     makepretty
    %
    %     subplot(1,3,3)
    %     scatter(PyKS,nanmean(Manual,1),10,[0 0 0],'filled')
    %     xlabel('PyKS')
    %     ylabel('Manual')
    %     xlim([-0.5 1.5])
    %     makepretty

    AvgMan = nanmean(Manual,1);

    % If PyKS said match
    disp('Compared to Stitched Kilosort:')
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==1)>0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==1)<0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a non-match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==0)>0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(PyKS==0)<0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a non-match'])

    % If UM said Match
    disp('Compared to UnitMatch:')
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb>0.5)>0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb>0.5)<0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a non-match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb<0.5)>0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a match'])
    disp(['Manual Scorers found ' num2str(round(sum(AvgMan(MatchProb<0.5)<0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a non-match'])


    if midx == 1
        clasfig = figure('name','Classifier Comparisons');
    else
        figure(clasfig)
    end
    subplot(2,2,1)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,PyKS==1)==1)./sum(PyKS==1),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==1)>0.5)./sum(PyKS==1),30,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(Rank(PyKS==1)<=2)./sum(PyKS==1),30,[0.75 0.75 0.75],'filled');

    title('Detection Performance matches relative to Kilosort Stitched')
    ylim([0 1])
    makepretty

    subplot(2,2,2)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,PyKS==0)==-1)./sum(PyKS==0),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==0)<0.5)./sum(PyKS==0),30,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(Rank(PyKS==0)>2)./sum(PyKS==0),30,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to Kilosort Stitched')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','Rank<3'})
    end
    makepretty

    subplot(2,2,3)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,MatchProb>0.5)==1)./sum(MatchProb>0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(PyKS(MatchProb>0.5)==1)./sum(MatchProb>0.5),30,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(Rank(MatchProb>0.5)<=2)./sum(MatchProb>0.5),30,[0.75 0.75 0.75],'filled');

    title('Detection Performance matches relative to UnitMatch')
    ylim([0 1])
    makepretty

    subplot(2,2,4)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,MatchProb<0.5)==-1)./sum(MatchProb<0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(PyKS(MatchProb<0.5)==0)./sum(MatchProb<0.5),30,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(Rank(MatchProb<0.5)>2)./sum(MatchProb<0.5),30,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to UnitMatch')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'Kilosort','Rank<3'})
    end
    makepretty


    if midx == 1
        clasfigMan = figure('name','Classifier Comparisons to Average Manual Score');
    else
        figure(clasfigMan)
    end
    subplot(1,3,1)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,AvgMan>0.5)==1)./sum(AvgMan>0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan>0.5)>0.5)./sum(AvgMan>0.5),30,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan>0.5)==1)./sum(AvgMan>0.5),30,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(Rank(AvgMan>0.5)<=2)./sum(AvgMan>0.5),30,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')
    title(['Average manual score match'])
    ylim([0 1])

    xlabel('Dataset')

    makepretty

    subplot(1,3,2)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,AvgMan<-0.5)==1)./sum(AvgMan<-0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan<-0.5)>0.5)./sum(AvgMan<-0.5),30,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan<-0.5)==1)./sum(AvgMan<-0.5),30,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(Rank(AvgMan<-0.5)<=2)./sum(AvgMan<-0.5),30,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')

    title('Average manual score non-match')
    ylim([0 1])

    xlabel('Dataset')

    makepretty

    subplot(1,3,3)
    if midx == 1
        line([0 length(MiceOpt)+1],[1 1],'color',[0 0 0],'LineStyle','--')
    end
    hold on
    clear h
    cols = lines(size(Manual,1));
    for id = 1:size(Manual,1)
        h(id) = scatter(midx+0.05*id,sum(Manual(id,abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),30,cols(id,:),'filled');
    end
    h(id+1) = scatter(midx,sum(MatchProb(abs(AvgMan)<=0.5)>0.5)./sum(abs(AvgMan)<=0.5),30,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),30,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(Rank(abs(AvgMan)<=0.5)<=2)./sum(abs(AvgMan)<=0.5),30,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','Kilosort','Rank<3'})
    end

    title('Average manual score uncertain')
    ylim([0 1])

    xlabel('Dataset')
    makepretty


    if midx == 1
        figFract = figure('name','Classifier Comparisons to Average Manual Score');
    else
        figure(figFract)
    end
    hold on
    clear h
    h(1) = scatter(midx,sum(AvgMan>0.5)./length(AvgMan),60,[0 0.5 0],'filled');
    h(2) = scatter(midx,sum(AvgMan<-0.5)./length(AvgMan),60,[0.5 0 0],'filled');
    h(3) = scatter(midx,sum(abs(AvgMan)<=0.5)./length(AvgMan),60,[0 0 0.5],'filled');
    title('Fraction')
    if midx == length(MiceOpt)
        legend(h,{'Everyone says match','Everyone says no match','Uncertain'})
    end

end
