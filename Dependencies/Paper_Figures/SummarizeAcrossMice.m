% SummarizeAcrossMice
FSCoreFig = figure('name','Functional Scores')

% INitialize
bins = -1:0.1:1;
Vector = [bins(1)+0.1/2:0.1:bins(end)-0.1/2];
FingerPrintRAcrossMice = nan(length(Vector),3,length(MiceOpt)); % Vector / within,Match,Non-match / mice
FingerprintAUC = nan(3,length(MiceOpt));
ACGRAcrossMice = nan(length(Vector),3,length(MiceOpt)); % Vector / within,Match,Non-match / mice
ACGAUC = nan(3,length(MiceOpt));
RFDistAcrossMice = nan(length(Vector),3,length(MiceOpt)); % Vector / within,Match,Non-match / mice
RFAUC = nan(3,length(MiceOpt));

for midx = 1:length(MiceOpt)
    tmpfile = dir(fullfile(SaveDir,MiceOpt{midx},'UnitMatch','UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end
    tmpFile = matfile(fullfile(tmpfile.folder,tmpfile.name));
    MatchTable = tmpFile.MatchTable; %Extract matchtable
    UMparam = tmpFile.UMparam; % Extract parameters

    % Extract groups
    WithinIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
    MatchIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
    NonMatchIdx = find((MatchTable.UID1 ~= MatchTable.UID2)); % Not the same unit

    % Extract cluster information
    UniqueIDConversion = tmpFile.UniqueIDConversion;
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

    %% Fingerprint correlation
    FingerprintCor = reshape(MatchTable.FingerprintCor,nclus,nclus);


    subplot(3,3,1)
    hold on
   
    hw = histcounts(FingerprintCor(WithinIdx),bins)./length(WithinIdx);
    hm = histcounts(FingerprintCor(MatchIdx),bins)./length(MatchIdx);
    hn = histcounts(FingerprintCor(NonMatchIdx),bins)./length(NonMatchIdx);
    FingerPrintRAcrossMice(:,:,midx) = cat(1,hw,hm,hn)';
    plot(Vector,hw,'color',[0 0 0.5])
    plot(Vector,hm,'color',[0 0.5 0])
    plot(Vector,hn,'color',[0.5 0 0])
    xlabel('Cross-correlation Fingerprint')
    ylabel('Proportion|Group')
    legend('i=j; within recording','matches','non-matches','Location','best')
    axis square
    makepretty

    subplot(3,3,2)
    hold on
    labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
    scores = [FingerprintCor(MatchIdx)', FingerprintCor(NonMatchIdx)'];
    [X,Y,~,AUC1] = perfcurve(labels,scores,1);
    h(1) = plot(X,Y,'color',[0.25 0.25 0]);
    hold all
    labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
    scores = [FingerprintCor(MatchIdx)', FingerprintCor(WithinIdx)'];
    [X,Y,~,AUC2] = perfcurve(labels,scores,1);
    h(2) = plot(X,Y,'color',[0 0.5 0.5]);
    labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
    scores = [FingerprintCor(WithinIdx)', FingerprintCor(NonMatchIdx)'];
    [X,Y,~,AUC3] = perfcurve(labels,scores,1);
    h(3) = plot(X,Y,'color',[0.5 0 0.5]);
    axis square
    FingerprintAUC(:,midx) = [AUC1,AUC2,AUC3];

    plot([0 1],[0 1],'k--')
    xlabel('False positive rate')
    ylabel('True positive rate')
    legend([h(:)],'Match vs No Match','Match vs Within','Within vs No Match','Location','best')
    title('Cross-Correlation Fingerprint')
    makepretty
    drawnow %Something to look at while ACG calculations are ongoing


    %% Autocorrelogram
    ACGCor = reshape(MatchTable.ACGCorr,nclus,nclus);

    subplot(3,3,4)
    hold on
   
    hw = histcounts(ACGCor(WithinIdx),bins)./length(WithinIdx);
    hm = histcounts(ACGCor(MatchIdx),bins)./length(MatchIdx);
    hn = histcounts(ACGCor(NonMatchIdx),bins)./length(NonMatchIdx);
    ACGRAcrossMice(:,:,midx) = cat(1,hw,hm,hn)';
    plot(Vector,hw,'color',[0 0 0.5])
    plot(Vector,hm,'color',[0 0.5 0])
    plot(Vector,hn,'color',[0.5 0 0])
    xlabel('Autocorrelogram correlation')
    ylabel('Proportion|Group')
    axis square
    makepretty

    subplot(3,3,5)
    hold on
    labels = [ones(1,numel(MatchIdx)), zeros(1,numel(NonMatchIdx))];
    scores = [ACGCor(MatchIdx)', ACGCor(NonMatchIdx)'];
    [X,Y,~,AUC1] = perfcurve(labels,scores,1);
    h(1) = plot(X,Y,'color',[0.25 0.25 0]);
    hold all
    labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
    scores = [ACGCor(MatchIdx)', ACGCor(WithinIdx)'];
    [X,Y,~,AUC2] = perfcurve(labels,scores,1);
    h(2) = plot(X,Y,'color',[0 0.5 0.5]);
    labels = [ones(1,numel(WithinIdx)), zeros(1,numel(NonMatchIdx))];
    scores = [ACGCor(WithinIdx)', ACGCor(NonMatchIdx)'];
    [X,Y,~,AUC3] = perfcurve(labels,scores,1);
    h(3) = plot(X,Y,'color',[0.5 0 0.5]);
    axis square
    ACGAUC(:,midx) = [AUC1,AUC2,AUC3];

    plot([0 1],[0 1],'k--')
    xlabel('False positive rate')
    ylabel('True positive rate')
    title('Auto-correlogram correlations')
    makepretty
    drawnow %Something to look at while ACG calculations are ongoing

    %% Receptive Field (?)
    if any(ismember(MatchTable.Properties.VariableNames,'RFDist'))
        RFDist = reshape(MatchTable.RFDist,nclus,nclus);
        subplot(3,3,7)
        hold on

        bins = linspace(0,50,length(Vector)+1);
        stepsz = unique(diff(bins));
        hw = histcounts(RFDist(WithinIdx),bins)./length(WithinIdx);
        hm = histcounts(RFDist(MatchIdx),bins)./length(MatchIdx);
        hn = histcounts(RFDist(NonMatchIdx),bins)./length(NonMatchIdx);
        RFDistAcrossMice(:,:,midx) = cat(1,hw,hm,hn)';
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2,hw,'color',[0 0 0.5])
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2,hm,'color',[0 0.5 0])
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2,hn,'color',[0.5 0 0])
        xlabel('Receptive Field Distance')
        ylabel('Proportion|Group')
        axis square
        makepretty

        subplot(3,3,8)
        hold on
        labels = [zeros(1,numel(MatchIdx)), ones(1,numel(NonMatchIdx))];
        scores = [RFDist(MatchIdx)', RFDist(NonMatchIdx)'];
        [X,Y,~,AUC1] = perfcurve(labels,scores,1);
        h(1) = plot(X,Y,'color',[0.25 0.25 0]);
        hold all
        labels = [ones(1,numel(MatchIdx)), zeros(1,numel(WithinIdx))];
        scores = [RFDist(MatchIdx)', RFDist(WithinIdx)'];
        [X,Y,~,AUC2] = perfcurve(labels,scores,1);
        h(2) = plot(X,Y,'color',[0 0.5 0.5]);
        labels = [zeros(1,numel(WithinIdx)), ones(1,numel(NonMatchIdx))];
        scores = [RFDist(WithinIdx)', RFDist(NonMatchIdx)'];
        [X,Y,~,AUC3] = perfcurve(labels,scores,1);
        h(3) = plot(X,Y,'color',[0.5 0 0.5]);
        axis square
        RFAUC(:,midx) = [AUC1,AUC2,AUC3];

        plot([0 1],[0 1],'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        title('RF Distance')
        makepretty
        drawnow %Something to look at while ACG calculations are ongoing


    end


end
%% Now we make histograms for the Area Under the Curve scores
subplot(3,3,3)
hold on
h1 = histogram(FingerprintAUC(1,:));
h1.FaceColor = [0.25 0.25 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25 0.25 0];

h2 = histogram(FingerprintAUC(2,:));
h2.FaceColor = [0 0.5 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0 0.5 0.5];

h3 = histogram(FingerprintAUC(3,:));
h3.FaceColor = [0.5 0 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5 0 0.5];

line([median(FingerprintAUC(1,:)) median(FingerprintAUC(1,:))],get(gca,'ylim'),'Color',[0.25 0.25 0])
line([median(FingerprintAUC(2,:)) median(FingerprintAUC(2,:))],get(gca,'ylim'),'Color',[0 0.5 0.5])
line([median(FingerprintAUC(3,:)) median(FingerprintAUC(3,:))],get(gca,'ylim'),'Color',[0.5 0 0.5])

xlabel('AUC cross-correlation')
ylabel('Nr sessions')
makepretty

%% 
subplot(3,3,6)
hold on
h1 = histogram(ACGAUC(1,:));
h1.FaceColor = [0.25 0.25 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25 0.25 0];

h2 = histogram(ACGAUC(2,:));
h2.FaceColor = [0 0.5 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0 0.5 0.5];

h3 = histogram(ACGAUC(3,:));
h3.FaceColor = [0.5 0 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5 0 0.5];

line([median(ACGAUC(1,:)) median(ACGAUC(1,:))],get(gca,'ylim'),'Color',[0.25 0.25 0])
line([median(ACGAUC(2,:)) median(ACGAUC(2,:))],get(gca,'ylim'),'Color',[0 0.5 0.5])
line([median(ACGAUC(3,:)) median(ACGAUC(3,:))],get(gca,'ylim'),'Color',[0.5 0 0.5])

xlabel('AUC ACG correlation')
ylabel('Nr sessions')
makepretty

%%
subplot(3,3,9)
hold on
h1 = histogram(RFAUC(1,:));
h1.FaceColor = [0.25 0.25 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25 0.25 0];

h2 = histogram(RFAUC(2,:));
h2.FaceColor = [0 0.5 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0 0.5 0.5];

h3 = histogram(RFAUC(3,:));
h3.FaceColor = [0.5 0 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5 0 0.5];

line([median(RFAUC(1,:)) median(RFAUC(1,:))],get(gca,'ylim'),'Color',[0.25 0.25 0])
line([median(RFAUC(2,:)) median(RFAUC(2,:))],get(gca,'ylim'),'Color',[0 0.5 0.5])
line([median(RFAUC(3,:)) median(RFAUC(3,:))],get(gca,'ylim'),'Color',[0.5 0 0.5])

xlabel('AUC RF distance')
ylabel('Nr sessions')
makepretty
