function ComputeFunctionalScores(SaveDir)

TmpFile = matfile(fullfile(SaveDir,'UnitMatch','UnitMatch.mat')); % Access saved file
UMparam = TmpFile.UMparam; % Extract parameters
UMparam.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows

MatchTable = TmpFile.MatchTable; % Load Matchtable

% Extract groups
WithinIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
MatchIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
NonMatchIdx = find((MatchTable.UID1 ~= MatchTable.UID2)); % Not the same unit

% Extract cluster information
UniqueIDConversion = TmpFile.UniqueIDConversion;
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
OriID = UniqueIDConversion.OriginalClusID(GoodId);
OriIDAll = UniqueIDConversion.OriginalClusID;
recses = UniqueIDConversion.recsesAll(GoodId);
recsesall = UniqueIDConversion.recsesAll;

AllKSDir = UMparam.KSDir; %original KS Dir
nclus = length(UniqueID);

% Load SP
disp('Loading spike information...')
nKSFiles = length(AllKSDir);
sp = cell(1,nKSFiles);
for did = 1:nKSFiles
    tmp = matfile(fullfile(AllKSDir{did},'PreparedData.mat'));
    sptmp = tmp.sp;
    clear tmp

    % Only keep parameters used
    sp{did}.st = sptmp.st;
    sp{did}.spikeTemplates = sptmp.spikeTemplates;
    sp{did}.spikeAmps = sptmp.spikeAmps;
    % Replace recsesid with subsesid
    sp{did}.RecSes = repmat(did,size(sp{did}.st));
end
clear sptmp


% Add all spikedata in one spikes struct - can be used for further analysis
sp = [sp{:}];
spnew = struct;
fields = fieldnames(sp(1));
for fieldid=1:length(fields)
    try
        eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
    catch ME
        if strcmp(ME.message,'Out of memory.')
            eval(['spnew.' fields{fieldid} ' = sp(1).' fields{fieldid} ';'])
            for tmpid = 2:length(sp)
                eval(['spnew.' fields{fieldid} ' = cat(1,spnew.' fields{fieldid} ', sp(tmpid).' fields{fieldid} ');'])
            end
        else
            eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])
        end
    end
end
sp = spnew;
clear spnew

if ~any(ismember(MatchTable.Properties.VariableNames,'FingerprintCor')) % If it already exists in table, skip this entire thing

    %% Compute cross-correlation matrices for individual recordings
    disp('Computing cross-correlation fingerprint')
    nRec = length(unique(UniqueIDConversion.recsesAll));
    SessionCorrelations = cell(1,nRec);
    for rid = 1:nRec
        % Define edges for this dataset
        edges = floor(min(sp.st(sp.RecSes==rid)))-UMparam.binsz/2:UMparam.binsz:ceil(max(sp.st(sp.RecSes==rid)))+UMparam.binsz/2;
        Good_Idx = find(GoodId & recsesall'==rid); % Only care about good units at this point

        % bin data to create PSTH
        sr = nan(numel(Good_Idx),numel(edges)-1);
        for uid = 1:numel(Good_Idx)
            sr(uid,:) =  histcounts(sp.st(sp.spikeTemplates == OriIDAll(Good_Idx(uid)) & sp.RecSes == recsesall(Good_Idx(uid))),edges);
        end

        % Define folds (two halves)
        idx_fold1 = 1:floor(size(sr,2)./2);
        idx_fold2 = floor(size(sr,2)./2)+1:floor(size(sr,2)./2)*2;

        % Find cross-correlation in first and second half of session
        SessionCorrelations{rid}.fold1 = corr(sr(:,idx_fold1)',sr(:,idx_fold1)')';
        SessionCorrelations{rid}.fold2 = corr(sr(:,idx_fold2)',sr(:,idx_fold2)')';

        % Nan the diagonal
        SessionCorrelations{rid}.fold1(logical(eye(size(SessionCorrelations{rid}.fold1)))) = nan;
        SessionCorrelations{rid}.fold2(logical(eye(size(SessionCorrelations{rid}.fold2)))) = nan;
    end

    %% Get correlation matrics for fingerprint correlations
    sessionCorrelationsAll = cell(1,nRec);
    for did = 1:nRec
        % Load sp for correct day
        if length(AllKSDir)>1
            tmp = matfile(fullfile(AllKSDir{did},'PreparedData.mat'));
        else %Stitched
            tmp = matfile(fullfile(AllKSDir{1},'PreparedData.mat'));
        end
        if length(SessionCorrelations)==nRec
            sessionCorrelationsAll{did} = SessionCorrelations{did};
        elseif length(SessionCorrelations)==1  %Normal situation
            if iscell(SessionCorrelations)
                sessionCorrelationsAll{did} = SessionCorrelations{1};
            else
                sessionCorrelationsAll{did} = SessionCorrelations;
            end
        else
            disp('This is a weird situation...')
            keyboard
        end
    end

    %% Compute functional score (cross-correlation fingerprint)
    if nRec<5
        drawdrosscorr = 1;
    else
        drawdrosscorr = 0;
    end
    MatchProbability = reshape(MatchTable.MatchProb,nclus,nclus);
    [r, c] = find(MatchProbability>UMparam.ProbabilityThreshold); %Find matches
    Pairs = cat(2,r,c);
    Pairs = sortrows(Pairs);
    Pairs = unique(Pairs,'rows');
    [FingerprintR,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(sessionCorrelationsAll,Pairs,OriID,recses,drawdrosscorr);
        
    % Save in table
    MatchTable.FingerprintCor = FingerprintR(:);
    MatchTable.RankScoreAll = RankScoreAll(:);
    MatchTable.SigFingerprintR = SigMask(:);

    TmpFile.Properties.Writable = true;
    TmpFile.MatchTable = MatchTable; % Overwrite
    TmpFile.Properties.Writable = false;

    %% Compare to functional scores
    EuclDist = reshape(MatchTable.EucledianDistance,nclus,nclus);
    SessionSwitch = arrayfun(@(X) find(recses==X,1,'first'),1:nRec,'Uni',0);
    SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
    SessionSwitch = [cell2mat(SessionSwitch) nclus+1];

    % Plotting order (sort units based on distance)
    [~,SortingOrder] = arrayfun(@(X) sort(EuclDist(1,SessionSwitch(X):SessionSwitch(X+1)-1)),1:nRec,'Uni',0);
    SortingOrder = arrayfun(@(X) squeeze(SortingOrder{X}+SessionSwitch(X)-1),1:nRec,'Uni',0);
    if size(SortingOrder{1},1)==1
        SortingOrder = cat(2,SortingOrder{:});
    else
        SortingOrder = cat(1,SortingOrder{:});
    end
    figure;
    
    subplot(1,3,1)
    imagesc(RankScoreAll(SortingOrder,SortingOrder)==1 & SigMask(SortingOrder,SortingOrder)==1)
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    title('Rankscore == 1*')
    makepretty

    subplot(1,3,2)
    imagesc(MatchProbability(SortingOrder,SortingOrder)>UMparam.ProbabilityThreshold)
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    title(['Match Probability>' num2str(UMparam.ProbabilityThreshold)])
    makepretty

    subplot(1,3,3)
    imagesc(MatchProbability(SortingOrder,SortingOrder)>=UMparam.ProbabilityThreshold | (MatchProbability(SortingOrder,SortingOrder)>0.05 & RankScoreAll(SortingOrder,SortingOrder)==1 & SigMask(SortingOrder,SortingOrder)==1));
    % imagesc(MatchProbability>=0.99 | (MatchProbability>=0.05 & RankScoreAll==1 & SigMask==1))
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    title('Matching probability + rank')
    makepretty
    saveas(gcf,fullfile(SaveDir,'RankScoreVSProbability.fig'))
    saveas(gcf,fullfile(SaveDir,'RankScoreVSProbability.bmp'))

    tmpf = triu(FingerprintR);
    tmpm = triu(MatchProbability);
    tmpr = triu(RankScoreAll);
    tmpr = tmpr(tmpf~=0);
    tmpm = tmpm(tmpf~=0);
    tmpf = tmpf(tmpf~=0);
    figure;
    scatter(tmpm,tmpf,14,tmpr,'filled')
    colormap(cat(1,[0 0 0],winter))
    xlabel('Match Probability')
    ylabel('Cross-correlation fingerprint')
    makepretty
    saveas(gcf,fullfile(SaveDir,'RankScoreVSProbabilityScatter.fig'))
    saveas(gcf,fullfile(SaveDir,'RankScoreVSProbabilityScatter.bmp'))
end

%% Cross-correlation ROC?
figure('name','Functional score separatability')

subplot(3,3,1)
imagesc(reshape(MatchTable.FingerprintCor,nclus,nclus))
hold on
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
freezeColors 

FingerprintCor = reshape(MatchTable.FingerprintCor,nclus,nclus);
Subtr = repmat(diag(FingerprintCor),1,size(FingerprintCor,1));
FingerprintCor = FingerprintCor - Subtr; % Subtract diagonalcorrelations  

subplot(3,3,2)
bins = min(FingerprintCor(:)):0.1:max(FingerprintCor(:));
Vector = [bins(1)+0.1/2:0.1:bins(end)-0.1/2];
hw = histcounts(FingerprintCor(WithinIdx),bins)./length(WithinIdx);
hm = histcounts(FingerprintCor(MatchIdx),bins)./length(MatchIdx);
hn = histcounts(FingerprintCor(NonMatchIdx),bins)./length(NonMatchIdx);
plot(Vector,hw,'color',[0.5 0.5 0.5])
hold on
plot(Vector,hm,'color',[0 0.5 0])
plot(Vector,hn,'color',[0 0 0])
xlabel('Cross-correlation Fingerprint')
ylabel('Proportion|Group')
legend('i=j; within recording','matches','non-matches','Location','best')
makepretty

subplot(3,3,3)
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
scores = [FingerprintCor(MatchIdx)', FingerprintCor(NonMatchIdx)'];
[X,Y,~,AUC1] = perfcurve(labels,scores,1);
h(1) = plot(X,Y,'color',[0 0.25 0]);
hold all
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
scores = [FingerprintCor(MatchIdx)', FingerprintCor(WithinIdx)'];
[X,Y,~,AUC2] = perfcurve(labels,scores,1);
h(2) = plot(X,Y,'color',[0 0.5 0]);
labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
scores = [FingerprintCor(WithinIdx)', FingerprintCor(NonMatchIdx)'];
[X,Y,~,AUC3] = perfcurve(labels,scores,1);
h(3) = plot(X,Y,'color',[0.25 0.25 0.25]);

plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
title(sprintf('Cross-Correlation Fingerprint AUC: %.3f, %.3f, %.3f', AUC1,AUC2,AUC3))
makepretty
drawnow %Something to look at while ACG calculations are ongoing

if ~any(ismember(MatchTable.Properties.VariableNames,'ACGCorr')) % If it already exists in table, skip this entire thing

    %% Compute ACG and correlate them between units
    % This is very time consuming
    disp('Computing ACG, this will take some time...')
    tvec =  -UMparam.ACGduration/2:UMparam.ACGbinSize:UMparam.ACGduration/2;
    ACGMat = nan(length(tvec),2,nclus);
    FR = nan(2,nclus);
    parfor clusid = 1:nclus
        for cv = 1:2
            idx1=find(sp.spikeTemplates == OriID(clusid) & sp.RecSes == recses(clusid));
            if ~isempty(idx1)
                if cv==1
                    idx1 = idx1(1:floor(length(idx1)/2));
                else
                    idx1 = idx1(ceil(length(idx1)/2):end);
                end

                % Compute Firing rate
                nspkspersec = histcounts(sp.st(idx1),[min(sp.st(idx1)):1:max(sp.st(idx1))]);
                FR(cv,clusid) = nanmean(nspkspersec);

                % compute ACG
                [ccg, t] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
                    ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
                ACGMat(:,cv,clusid) = ccg(:, 1, 1);
            end
        end
    end

    %% Correlation between ACG
    ACGCorr = corr(squeeze(ACGMat(:,1,:)),squeeze(ACGMat(:,2,:)));
    MatchTable.ACGCorr = ACGCorr(:);

    %% FR difference
    FR = repmat(permute(FR,[2,1]),[1,1,nclus]);
    FRDiff = abs(squeeze(FR(:,2,:) - permute(FR(:,1,:),[3,2,1])));
    MatchTable.FRDiff = FRDiff(:);

    % Write to table
    TmpFile.Properties.Writable = true;
    TmpFile.MatchTable = MatchTable; % Overwrite
    TmpFile.Properties.Writable = false;

end

%% Plot ACG
subplot(3,3,4)
imagesc(reshape(MatchTable.ACGCorr,nclus,nclus))
hold on
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
freezeColors 

subplot(3,3,5)
ACGCor = reshape(MatchTable.ACGCorr,nclus,nclus);
Subtr = repmat(diag(ACGCor),1,size(ACGCor,1));
ACGCor = ACGCor - Subtr; % Subtract diagonalcorrelations  
bins = min(ACGCor(:)):0.1:max(ACGCor(:));
Vector = [bins(1)+0.1/2:0.1:bins(end)-0.1/2];
hw = histcounts(ACGCor(WithinIdx),bins)./length(WithinIdx);
hm = histcounts(ACGCor(MatchIdx),bins)./length(MatchIdx);
hn = histcounts(ACGCor(NonMatchIdx),bins)./length(NonMatchIdx);
plot(Vector,hw,'color',[0.5 0.5 0.5])
hold on
plot(Vector,hm,'color',[0 0.5 0])
plot(Vector,hn,'color',[0 0 0])
xlabel('Autocorrelogram Correlation')
ylabel('Proportion|Group')
legend('i=j; within recording','matches','non-matches','Location','best')
makepretty

subplot(3,3,6)
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
scores = [ACGCor(MatchIdx)', ACGCor(NonMatchIdx)'];
[X,Y,~,AUC1] = perfcurve(labels,scores,1);
h(1) = plot(X,Y,'color',[0 0.25 0]);
hold all
labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
scores = [ACGCor(MatchIdx)', ACGCor(WithinIdx)'];
[X,Y,~,AUC2] = perfcurve(labels,scores,1);
h(2) = plot(X,Y,'color',[0 0.5 0]);
labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
scores = [ACGCor(WithinIdx)', ACGCor(NonMatchIdx)'];
[X,Y,~,AUC3] = perfcurve(labels,scores,1);
h(3) = plot(X,Y,'color',[0.25 0.25 0.25]);

plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
title(sprintf('Autocorrelogram AUC: %.3f, %.3f, %.3f', AUC1,AUC2,AUC3))
makepretty
freezeColors 

%% Plot FR
subplot(3,3,7)
imagesc(reshape(MatchTable.FRDiff,nclus,nclus))
hold on
colormap(gray)
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)

FRDiff = reshape(MatchTable.FRDiff,nclus,nclus);
Subtr = repmat(diag(FRDiff),1,size(FRDiff,1));
FRDiff = FRDiff - Subtr; % Subtract diagonalcorrelations 

subplot(3,3,8)
bins = min(FRDiff(:)):0.1:5;
Vector = [bins(1)+0.1/2:0.1:bins(end)-0.1/2];
hw = histcounts(FRDiff(WithinIdx),bins)./length(WithinIdx);
hm = histcounts(FRDiff(MatchIdx),bins)./length(MatchIdx);
hn = histcounts(FRDiff(NonMatchIdx),bins)./length(NonMatchIdx);
plot(Vector,hw,'color',[0.5 0.5 0.5])
hold on
plot(Vector,hm,'color',[0 0.5 0])
plot(Vector,hn,'color',[0 0 0])
xlabel('Firing rate differences')
ylabel('Proportion|Group')
legend('i=j; within recording','matches','non-matches','Location','best')
makepretty

subplot(3,3,9)
labels = [zeros(1,numel(MatchIdx)), ones(1,numel(NonMatchIdx))];
scores = [FRDiff(MatchIdx)', FRDiff(NonMatchIdx)'];
[X,Y,~,AUC1] = perfcurve(labels,scores,1);
h(1) = plot(X,Y,'color',[0 0.25 0]);
hold all
labels = [zeros(1,numel(MatchIdx)), ones(1,numel(WithinIdx))];
scores = [FRDiff(MatchIdx)', FRDiff(WithinIdx)'];
[X,Y,~,AUC2] = perfcurve(labels,scores,1);
h(2) = plot(X,Y,'color',[0 0.5 0]);
labels = [zeros(1,numel(WithinIdx)), ones(1,numel(NonMatchIdx))];
scores = [FRDiff(WithinIdx)', FRDiff(NonMatchIdx)'];
[X,Y,~,AUC3] = perfcurve(labels,scores,1);
h(3) = plot(X,Y,'color',[0.25 0.25 0.25]);

plot([0 1],[0 1],'k--')
xlabel('False positive rate')
ylabel('True positive rate')
legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
title(sprintf('Firing rate differences AUC: %.3f, %.3f, %.3f', AUC1,AUC2,AUC3))
makepretty

% save
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,fullfile(SaveDir,'FunctionalScoreSeparability.fig'))
saveas(gcf,fullfile(SaveDir,'FunctionalScoreSeparability.png'))
