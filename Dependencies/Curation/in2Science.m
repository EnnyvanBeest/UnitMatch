%% In2Science
% 3 high-school students were with us for a week. They were given the
% assignment to decide if two different waveforms came from the same neuron
% yes or no. Here is how they judged that compared to lab members, and
% compared to UnitMatch and KiloSort
SaveDir = 'H:\MatchingUnits\Output\Concatenated1Day\' % This is the directory where we saved the data
MiceOpt = {'AL032'} % This is the mouse we have compared you to
% Initialize:
Check = cell(4,length(MiceOpt));


% We're going to loop over mice
for midx = 1:length(MiceOpt)
    % compare

    % Load the data
    tmpdir = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
    tmpMatchTbl = matfile(fullfile(tmpdir.folder,tmpdir.name));
    MatchTable = tmpMatchTbl.MatchTable;

    % Extract information: how did the model do?
    nclus = sqrt(height(MatchTable));
    MatchProbability = reshape(MatchTable.MatchProb,nclus,nclus);
    RankScore = reshape(MatchTable.RankScore,nclus,nclus);
    RankThreshold = reshape(MatchTable.SigFingerprintR,nclus,nclus);
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

    % Individual scores:
    Scores2Inclde = UMparam.Scores2Include;
    % Add functional scores
    Scores2Inclde = {Scores2Inclde{:},'FingerprintCor'};
    for scid = 1:length(Scores2Inclde)
        eval([Scores2Inclde{scid} ' = reshape(MatchTable.' Scores2Inclde{scid} ',nclus,nclus);'])
    end

    
    for fid = 1:2
        % Now load manual curation
        if fid == 1
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures_WOFunct','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end
            ScorerNames = cell(length(manualscore),1);
          
            % Pairs of Interest
            tbl = tmpmanual.tbl; % We need this information for later


            % UM results
            Pair1 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID1(X)) & ismember(recses,tbl.RecID1(X))),1:height(tbl),'Uni',0));
            Pair2 = cell2mat(arrayfun(@(X) find(ismember(OriID,tbl.ClusID2(X)) & ismember(recses,tbl.RecID2(X))),1:height(tbl),'Uni',0));
            Pairs = arrayfun(@(X) [Pair1(X) Pair2(X)],1:length(Pair1),'Uni',0);
            MatchProb = cell2mat(cellfun(@(X) MatchProbability(X(1),X(2)),Pairs,'Uni',0))';
            Rank = cell2mat(cellfun(@(X) RankScore(X(1),X(2)),Pairs,'Uni',0))';
            RankThreshold = cell2mat(cellfun(@(X) RankThreshold(X(1),X(2)),Pairs,'Uni',0))';
            WithinSameSession = cell2mat(cellfun(@(X) recses(X(1)) == recses(X(2)),Pairs,'Uni',0));

            % Py KS results
            PyKS = cell2mat(cellfun(@(X) OriID(X(1)) == OriID(X(2)),Pairs,'Uni',0));

            % Individual scores
            for scid = 1:length(Scores2Inclde)
                eval([Scores2Inclde{scid} ' = cell2mat(cellfun(@(X) ' Scores2Inclde{scid} '(X(1),X(2)),Pairs,''Uni'',0));'])
            end

        else
            % Extract randomization for the blind images
            tmpmanual = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','BlindTable.mat'));
            tmpmanual = load(fullfile(tmpmanual.folder,tmpmanual.name));
            tbl2 = tmpmanual.tbl;
            % Order?
            Order = cell2mat(arrayfun(@(X) find(tbl2.ClusID1 == tbl.ClusID1(X) & tbl2.ClusID2 == tbl.ClusID2(X) & tbl2.RecID1 == tbl.RecID1(X) & tbl2.RecID2 == tbl.RecID2(X),1,'first'),1:height(tbl),'Uni',0));
            ReOrderedTbl = tbl2(Order,:);
            manualscore = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','BlindFigures','manualCuration*.mat'));
            if isempty(manualscore)
                continue
            end

            % Actually extract scores from different human scorers
            for id = 1:length(manualscore)

                tmpman = load(fullfile(manualscore(id).folder,manualscore(id).name));
                if id==1
                    Manual = tmpman.match(Order);
                else
                    Manual = cat(1,Manual,tmpman.match(Order));
                end
                tmpname = strsplit(manualscore(id).name,'manualCuration_');
                tmpname = strsplit(tmpname{2},'.mat');

                ScorerNames{id,fid} = [tmpname{1}];
            end
        end
    end
    ScorerNames = ScorerNames(~cellfun(@isempty,ScorerNames));

    %% How well do the different scoring methods correlate?
    % You didn't score everything, so we need to only look at the ones you
    % scored
    In2Science = Manual(ismember(ScorerNames,'In2Science'),:);
    idx = find(Manual(ismember(ScorerNames,'In2Science'),:)~=0); 
    AvgMan = nanmean(Manual,1);
    disp('User scoring correlations:')
    AllScoringMethods = cat(2,Manual(:,idx)',AvgMan(idx)',MatchProb(idx),PyKS(idx)',Rank(idx)==1,RankThreshold(idx));
    % Normalize between 0 and 1
    AllScoringMethods = (AllScoringMethods - nanmin(AllScoringMethods,[],1))./(nanmax(AllScoringMethods,[],1)-nanmin(AllScoringMethods,[],1));
    
    %% Correlation
    tmpcor = corr(AllScoringMethods);
    AllScorerNames = {ScorerNames{:},'AvgScorer','UM','KS','Rank1','RankTr'};
    if midx == 1
        corrfig = figure('name','Scorer correlations');
    else
        figure(corrfig)
    end
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    h=imagesc(tmpcor,[0,1]);
    set(gca,'XTick',1:size(tmpcor,1),'XTickLabel',{ScorerNames{:},'AvgScorer','UM','KS','Rank1','RankTr'},'YTick',1:size(tmpcor,1),'YTickLabel',{ScorerNames{:},'AvgScorer', 'UM','KS','Rank1','RankTr'})
    colormap(flipud(gray))
    title(MiceOpt{midx})
    colorbar
    axis square
    makepretty


    %% Any particular parameter correlations?
    if midx == 1
        % Save out correlation of each method with each of the parameters:
        CorrParam = nan(length(Scores2Inclde),length(AllScorerNames),length(MiceOpt));
        corrperparamfig = figure('name','Corerlations with parameters');
    else
        figure(corrperparamfig);
    end
    for scid = 1:length(Scores2Inclde)
%         eval(['Idx = ~isnan(' Scores2Inclde{scid} ') & Manual(ismember(ScorerNames,''In2Science''),:)~=0;'])
        eval(['CorrParam(scid,:,midx) = corr(' Scores2Inclde{scid} '(idx)'',AllScoringMethods);'])
    end

    % Invert correlation for FR Diff (as less different = more likely
    % match)
    CorrParam(ismember(Scores2Inclde,'FRDiff'),:,midx) = -1*CorrParam(ismember(Scores2Inclde,'FRDiff'),:,midx);

    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    h=imagesc(CorrParam(:,:,midx)',[0,1]);
    set(gca,'YTick',1:length(AllScorerNames),'YTickLabel',AllScorerNames,'XTick',1:length(Scores2Inclde),'XTickLabel',Scores2Inclde)
    colormap(flipud(gray))
    title(MiceOpt{midx})
    colorbar
    axis square
    makepretty

    %% If PyKS said match
    In2Science = In2Science(idx);
    PyKS = PyKS(idx);
    MatchProb = MatchProb(idx);
    AvgMan = nanmean(Manual(~ismember(ScorerNames,'In2Science'),:),1);
    AvgMan = AvgMan(idx);
    Manual = Manual(:,idx);
    RankThreshold = RankThreshold(idx);
    disp('Compared to Stitched Kilosort:')
    disp(['In2Science students found ' num2str(round(sum(In2Science(PyKS==1)>0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(PyKS==1)<0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a non-match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(PyKS==0)>0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(PyKS==0)<0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a non-match'])

    disp('Compared to Stitched Kilosort:')
    disp(['Lab members found ' num2str(round(sum(AvgMan(PyKS==1)>0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(PyKS==1)<0)./sum(PyKS==1)*1000)/10) '% of PyKS Matches to be a non-match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(PyKS==0)>0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(PyKS==0)<0)./sum(PyKS==0)*1000)/10) '% of PyKS Non-matches to be a non-match'])

    % If UM said Match
    disp('Compared to UnitMatch:')
    disp(['In2Science students found ' num2str(round(sum(In2Science(MatchProb>0.5)>0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(MatchProb>0.5)<0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a non-match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(MatchProb<0.5)>0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a match'])
    disp(['In2Science students found ' num2str(round(sum(In2Science(MatchProb<0.5)<0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a non-match'])

    disp(['Lab members found ' num2str(round(sum(AvgMan(MatchProb>0.5)>0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(MatchProb>0.5)<0)./sum(MatchProb>0.5)*1000)/10) '% of UnitMatch Matches to be a non-match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(MatchProb<0.5)>0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a match'])
    disp(['Lab members found ' num2str(round(sum(AvgMan(MatchProb<0.5)<0)./sum(MatchProb<0.5)*1000)/10) '% of UnitMatch Non-matches to be a non-match'])


    %% Classifier comparisons
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
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==1)>0.5)./sum(PyKS==1),50,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(PyKS==1)==1)./sum(PyKS==1),50,[0.75 0.75 0.75],'filled');

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
    h(id+1) = scatter(midx,sum(MatchProb(PyKS==0)<0.5)./sum(PyKS==0),50,[0 0 0],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(PyKS==0)==0)./sum(PyKS==0),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to Kilosort Stitched')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','RankThrs'})
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
    h(id+1) = scatter(midx,sum(PyKS(MatchProb>0.5)==1)./sum(MatchProb>0.5),50,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(MatchProb>0.5)==1)./sum(MatchProb>0.5),50,[0.75 0.75 0.75],'filled');

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
    h(id+1) = scatter(midx,sum(PyKS(MatchProb<0.5)==0)./sum(MatchProb<0.5),50,[0.5 0.5 0.5],'filled');
    h(id+2) = scatter(midx,sum(RankThreshold(MatchProb<0.5)==0)./sum(MatchProb<0.5),50,[0.75 0.75 0.75],'filled');

    title('Detection Performance non-matches relative to UnitMatch')
    ylim([0 1])
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'Kilosort','RankThrs'})
    end
    makepretty


    %% Against average manual
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
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan>0.5)>0.5)./sum(AvgMan>0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan>0.5)==1)./sum(AvgMan>0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(AvgMan>0.5)==1)./sum(AvgMan>0.5),50,[0.75 0.75 0.75],'filled');

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
    h(id+1) = scatter(midx,sum(MatchProb(AvgMan<-0.5)>0.5)./sum(AvgMan<-0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(AvgMan<-0.5)==1)./sum(AvgMan<-0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(AvgMan<-0.5)==1)./sum(AvgMan<-0.5),50,[0.75 0.75 0.75],'filled');

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
    h(id+1) = scatter(midx,sum(MatchProb(abs(AvgMan)<=0.5)>0.5)./sum(abs(AvgMan)<=0.5),50,[0 0 0],'filled');
    h(id+2) = scatter(midx-0.1,sum(PyKS(abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),50,[0.5 0.5 0.5],'filled');
    h(id+3) = scatter(midx-0.05,sum(RankThreshold(abs(AvgMan)<=0.5)==1)./sum(abs(AvgMan)<=0.5),50,[0.75 0.75 0.75],'filled');

    ylabel('Defined as match (proportion)')
    if midx == length(MiceOpt)
        legend([h(:)],{ScorerNames{:} 'UnitMatch','Kilosort','RankThrs'})
    end

    title('Average manual score uncertain')
    ylim([0 1])

    xlabel('Dataset')
    makepretty


end
