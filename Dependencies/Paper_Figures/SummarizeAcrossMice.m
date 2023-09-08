% SummarizeAcrossMice
FSCoreFig = figure('name', 'Functional Scores');

% INitialize
bins = -1:0.1:1;
Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
FingerPrintRAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
FingerprintAUC = nan(3, length(MiceOpt));
ACGRAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
ACGAUC = nan(3, length(MiceOpt));
RFDistAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
RFAUC = nan(3, length(MiceOpt));
EPosAndNeg = nan(3, length(MiceOpt));
UseKSLabels = PrepareClusInfoparams.RunPyKSChronicStitched;
RecSesPerUIDAllMice = cell(1,length(MiceOpt));
for midx = 1:length(MiceOpt)
    tmpfile = dir(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end
    tmpFile = matfile(fullfile(tmpfile.folder, tmpfile.name));
    MatchTable = tmpFile.MatchTable; %Extract matchtable
    UMparam = tmpFile.UMparam; % Extract parameters
    UniqueIDConversion = tmpFile.UniqueIDConversion; % Extract UniqueID Conversion

    % Load AUCS
    if exist(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch', 'AUC.mat'))
        AUC = load(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch', 'AUC.mat'))';
        if midx == 1 || ~exist('AUCParams','var')
            AUCParams = AUC.AUCStruct.ParamNames;
            AUCVals = nan(length(AUCParams), length(MiceOpt));
        end
        AUCVals(:, midx) = AUC.AUCStruct.AUC;
    end

    % Extract groups
    if ~UseKSLabels
        WithinIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
        MatchIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
        NonMatchIdx = find((MatchTable.UID1 ~= MatchTable.UID2)); % Not the same unit
    else
        WithinIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
        MatchIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
        NonMatchIdx = find((MatchTable.ID1 ~= MatchTable.ID2)); % Not the same unit
    end
    % Extract cluster information
    if UMparam.GoodUnitsOnly
        GoodId = logical(UniqueIDConversion.GoodID);
    else
        GoodId = true(1, length(UniqueIDConversion.GoodID));
    end
    UniqueID = UniqueIDConversion.UniqueID(GoodId);
    OriID = UniqueIDConversion.OriginalClusID(GoodId);
    OriIDAll = UniqueIDConversion.OriginalClusID;
    recses = UniqueIDConversion.recsesAll(GoodId);
    recsesall = UniqueIDConversion.recsesAll;
    nRecs = length(unique(recsesall));
    AllKSDir = UMparam.KSDir; %original KS Dir
    AllRawDir = UMparam.AllRawPaths; %
    if isstruct(AllRawDir)
        AllRawDir = arrayfun(@(X) fullfile(X.folder,X.name),AllRawDir,'Uni',0);
    end
    nclus = length(UniqueID);
    ChannelPos = UMparam.channelpos;

    %% How many units vs number of units were tracked?
    TrackingPerformance = nan(4, 0); % MatchesExpected, Difference between recording day, % Tracked units %maximum possibility
    TrackingPerformanceKS = nan(4, 0); % MatchesExpected, Difference between recording day, % Tracked units %maximum possibility
    AllDeltaDays = nan(1,nRecs);
    ProbeSN = cell(1,nRecs); % To make sure the same probe was used, read the Serial number from the params file an store.
    for did1 = 1:nRecs

        % Read ProbeID
        if isempty(ProbeSN{did1})
            if any(strfind(AllRawDir{did1},'zinu')) && ~any(strfind(AllRawDir{did1},'\\zinu.cortexlab.net\'))
                 AllRawDir{did1} = strrep(AllRawDir{did1},AllRawDir{did1}(1:strfind(AllRawDir{did1},'zinu')+3),'\\zinu.cortexlab.net\subjects');
            end
            rawdir = dir(fullfile(AllRawDir{did1}));
            meta = ReadMeta2(fullfile(rawdir.folder,'*.ap.meta'));
            try
                ProbeSN{did1} = meta.imDatPrb_sn;
            catch ME
                disp(ME)
                keyboard
            end
        end

       
        % For difference in days
        tmpday1 = strsplit(fullfile(AllRawDir{did1}),'\');% relies on yyyy-mm-dd format
        tmpday1 = tmpday1(cellfun(@(X) any(strfind(X,'-')),tmpday1));
        try
        tmpday1 = datetime(tmpday1{1},'InputFormat','yyyy-MM-dd');
        catch
            tmpday1 = datetime(tmpday1{1},'InputFormat','dd-MM-yyyy');
        end

        if did1==1
            firstday = tmpday1;
        end
        AllDeltaDays(did1) = double(days(duration(tmpday1-firstday)));

        if length(ChannelPos)>1
            ChansD1 = ChannelPos{did1};
        else
            ChansD1 = ChannelPos{1};
        end

        for did2 = 1:nRecs
            if did2 <= did1
                continue
            end

            % Read ProbeID
            if isempty(ProbeSN{did2})
                if any(strfind(AllRawDir{did2},'zinu')) && ~any(strfind(AllRawDir{did2},'\\zinu.cortexlab.net\'))
              
                    AllRawDir{did2} = strrep(AllRawDir{did2},AllRawDir{did2}(1:strfind(AllRawDir{did2},'zinu')+3),'\\zinu.cortexlab.net\subjects');
                end
                rawdir = dir(fullfile(AllRawDir{did2}));
                meta = ReadMeta2(fullfile(rawdir.folder,'*.ap.meta'));
                try
                    ProbeSN{did2} = meta.imDatPrb_sn;
                catch ME
                    disp(ME)
                    keyboard
                end
            end
            tmpday2 = strsplit(fullfile(AllRawDir{did2}),'\');% relies on yyyy-mm-dd format
            tmpday2 = tmpday2(cellfun(@(X) any(strfind(X,'-')),tmpday2));
            try
                tmpday2 = datetime(tmpday2{1},'InputFormat','yyyy-MM-dd');
            catch
                tmpday2 = datetime(tmpday2{1},'InputFormat','dd-MM-yyyy');
            end

            % Delta Days
            DDay = days(duration(tmpday2-tmpday1));

            if length(ChannelPos)>1
                ChansD2 = ChannelPos{did2};
            else
                ChansD2 = ChannelPos{1};
            end

            if strcmp(RecordingType{midx},'Chronic') & ProbeSN{did1} == ProbeSN{did2}% How much did the IMRO overlay??
                MatchesExpected = sum(ChansD1(:) == ChansD2(:))./numel(ChansD1);
            else
                MatchesExpected = 0;
            end

            thesedaysidx = find(ismember(recses, [did1, did2]));
            % can possibly only track these many units:
            nMax = min([sum(recses(thesedaysidx) == did1), sum(recses(thesedaysidx) == did2)]); 
            if nMax < 25
                continue
            end
            if sum(recses(thesedaysidx) == did1)>sum(recses(thesedaysidx) == did2)
                nMatches = sum(ismember(UniqueID(recses==did2),UniqueID(recses==did1)));
            else
                nMatches = sum(ismember(UniqueID(recses==did1),UniqueID(recses==did2)));
            end
            if nMatches>nMax
                keyboard
            end
            TrackingPerformance = cat(2, TrackingPerformance, [MatchesExpected, double(DDay), nMatches, nMax]');
            if MatchesExpected == 1 && (nMatches/nMax)==0
                warning(['No matches unexpectedly for ' MiceOpt{midx}])
            end

            if sum(recses(thesedaysidx) == did1)>sum(recses(thesedaysidx) == did2)
                nMatches = sum(ismember(OriID(recses==did2),OriID(recses==did1)));
            else
                nMatches = sum(ismember(OriID(recses==did1),OriID(recses==did2)));
            end
            if UseKSLabels
                TrackingPerformanceKS = cat(2, TrackingPerformanceKS, [MatchesExpected, double(DDay), nMatches, nMax]');
            end
        end
    end

    UMTrackingPerformancePerMouse{midx} = TrackingPerformance;
    if UseKSLabels
        KSTrackingPerformancePerMouse{midx} = TrackingPerformanceKS;
    end



    %% UniqueID X tracking
    [UniqueIDOpt,idx1,idx2] = unique(UniqueID); %UID options
    DayOpt = unique(AllDeltaDays);
    dayses = AllDeltaDays(recses); % Get in days after first recording
    RecSesPerUID = arrayfun(@(X) ismember(DayOpt,dayses(idx2==X)),1:numel(UniqueIDOpt),'Uni',0); % Extract which recording sessions a unite appears in
    RecSesPerUID = cat(1,RecSesPerUID{:});

    if midx == 1
        MatrixFig = figure('name','TrackingAcrossDaysMatrix');
    else
        figure(MatrixFig)
    end
    subplot(ceil(sqrt(length(MiceOpt))),round(sqrt(length(MiceOpt))),midx)
    h = imagesc(RecSesPerUID(sum(RecSesPerUID,2)>1,:));
    colormap(flipud(gray))
    title(MiceOpt{midx})
    xlabel('\Deltadays')
    ylabel('Tracked Units')
    makepretty
    RecSesPerUIDAllMice{midx} = RecSesPerUID; % Save for later use

    %% Extra Positives/negatives
    MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
    %Extra Positive
    EPosAndNeg(1, midx) = sum(MatchProb(NonMatchIdx) > UMparam.ProbabilityThreshold) ./ length(NonMatchIdx);

    %Extra Negative
    EPosAndNeg(2, midx) = sum(MatchProb(WithinIdx) < UMparam.ProbabilityThreshold) ./ length(WithinIdx);

    if UseKSLabels
        EPosAndNeg(3, midx) = sum(MatchProb(MatchIdx) < UMparam.ProbabilityThreshold) ./ length(MatchIdx);
    end

    %% Fingerprint correlation
    if any(ismember(MatchTable.Properties.VariableNames, 'FingerprintCor'))
    figure(FSCoreFig)
    FingerprintCor = reshape(MatchTable.FingerprintCor, nclus, nclus);


    subplot(3, 3, 1)
    hold on

    hw = histcounts(FingerprintCor(WithinIdx), bins) ./ length(WithinIdx);
    hm = histcounts(FingerprintCor(MatchIdx), bins) ./ length(MatchIdx);
    hn = histcounts(FingerprintCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
    FingerPrintRAcrossMice(:, :, midx) = cat(1, hw, hm, hn)';
    plot(Vector, hw, 'color', [0, 0, 0.5])
    plot(Vector, hm, 'color', [0, 0.5, 0])
    plot(Vector, hn, 'color', [0.5, 0, 0])
    xlabel('Cross-correlation Fingerprint')
    ylabel('Proportion|Group')
    legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
    axis square
    makepretty

    subplot(3, 3, 2)
    hold on
    clear h
    if length(MatchIdx) > 10
        labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [FingerprintCor(MatchIdx)', FingerprintCor(NonMatchIdx)'];
        [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
        h(1) = plot(X, Y, 'color', [0.25, 0.25, 0]);
        hold all
        labels = [ones(1, numel(MatchIdx)), zeros(1, numel(WithinIdx))];
        scores = [FingerprintCor(MatchIdx)', FingerprintCor(WithinIdx)'];
        [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
        h(2) = plot(X, Y, 'color', [0, 0.5, 0.5]);
    else
        AUC1 = nan;
        AUC2 = nan;
    end

    labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
    scores = [FingerprintCor(WithinIdx)', FingerprintCor(NonMatchIdx)'];
    [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
    h(3) = plot(X, Y, 'color', [0.5, 0, 0.5]);
    axis square
    FingerprintAUC(:, midx) = [AUC1, AUC2, AUC3];

    plot([0, 1], [0, 1], 'k--')
    xlabel('False positive rate')
    ylabel('True positive rate')
    if length(MatchIdx) > 10
        legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
    else
        legend([h(3)], 'Within vs No Match', 'Location', 'best')
    end
    title('Cross-Correlation Fingerprint')
    makepretty
    drawnow %Something to look at while ACG calculations are ongoing
    end
    %% Autocorrelogram
    if any(ismember(MatchTable.Properties.VariableNames, 'ACGCorr'))
        ACGCor = reshape(MatchTable.ACGCorr, nclus, nclus);

        subplot(3, 3, 4)
        hold on

        hw = histcounts(ACGCor(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(ACGCor(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(ACGCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
        ACGRAcrossMice(:, :, midx) = cat(1, hw, hm, hn)';
        plot(Vector, hw, 'color', [0, 0, 0.5])
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0.5, 0, 0])
        xlabel('Autocorrelogram correlation')
        ylabel('Proportion|Group')
        axis square
        makepretty

        subplot(3, 3, 5)
        hold on
        labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [ACGCor(MatchIdx)', ACGCor(NonMatchIdx)'];
        [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
        h(1) = plot(X, Y, 'color', [0.25, 0.25, 0]);
        hold all
        labels = [ones(1, numel(MatchIdx)), zeros(1, numel(WithinIdx))];
        scores = [ACGCor(MatchIdx)', ACGCor(WithinIdx)'];
        [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
        h(2) = plot(X, Y, 'color', [0, 0.5, 0.5]);
        labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [ACGCor(WithinIdx)', ACGCor(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.5, 0, 0.5]);
        axis square
        ACGAUC(:, midx) = [AUC1, AUC2, AUC3];

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        title('Auto-correlogram correlations')
        makepretty
        drawnow %Something to look at while ACG calculations are ongoing
    end

    %% Receptive Field (?)
    if any(ismember(MatchTable.Properties.VariableNames, 'RFDist'))
        RFDist = reshape(MatchTable.RFDist, nclus, nclus);
        subplot(3, 3, 7)
        hold on

        bins = linspace(0, 50, length(Vector)+1);
        stepsz = unique(diff(bins));
        hw = histcounts(RFDist(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(RFDist(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(RFDist(NonMatchIdx), bins) ./ length(NonMatchIdx);
        RFDistAcrossMice(:, :, midx) = cat(1, hw, hm, hn)';
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hw, 'color', [0, 0, 0.5])
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hm, 'color', [0, 0.5, 0])
        plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hn, 'color', [0.5, 0, 0])
        xlabel('Receptive Field Distance')
        ylabel('Proportion|Group')
        axis square
        makepretty

        subplot(3, 3, 8)
        hold on
        labels = [zeros(1, numel(MatchIdx)), ones(1, numel(NonMatchIdx))];
        scores = [RFDist(MatchIdx)', RFDist(NonMatchIdx)'];
        [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
        h(1) = plot(X, Y, 'color', [0.25, 0.25, 0]);
        hold all
        labels = [ones(1, numel(MatchIdx)), zeros(1, numel(WithinIdx))];
        scores = [RFDist(MatchIdx)', RFDist(WithinIdx)'];
        [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
        h(2) = plot(X, Y, 'color', [0, 0.5, 0.5]);
        labels = [zeros(1, numel(WithinIdx)), ones(1, numel(NonMatchIdx))];
        scores = [RFDist(WithinIdx)', RFDist(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.5, 0, 0.5]);
        axis square
        RFAUC(:, midx) = [AUC1, AUC2, AUC3];

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        title('RF Distance')
        makepretty
        drawnow %Something to look at while ACG calculations are ongoing


    end


end

%% AUC

if exist('AUCVals')
    meanAUC = nanmean(AUCVals,2);
    [~,sortidx] = sort(meanAUC,'descend');
    figure; h=barwitherr(nanstd(AUCVals(sortidx,:),[],2),nanmean(AUCVals(sortidx,:),2));

    set(gca, 'XTick', 1:size(AUCVals, 1), 'XTickLabel', AUCParams(sortidx))
    ylabel('AUC')
    makepretty

    saveas(gcf,fullfile(SaveDir,'AUCParameters.fig'))
    saveas(gcf,fullfile(SaveDir,'AUCParameters.bmp'))
end
%% Now we make histograms for the Area Under the Curve scores
figure(FSCoreFig)
edges = [0:0.1:1];
subplot(3, 3, 3)
hold on
h1 = histogram(FingerprintAUC(1, :), edges);
h1.FaceColor = [0.25, 0.25, 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25, 0.25, 0];

h2 = histogram(FingerprintAUC(2, :), edges);
h2.FaceColor = [0, 0.5, 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0, 0.5, 0.5];


h3 = histogram(FingerprintAUC(3, :), edges);
h3.FaceColor = [0.5, 0, 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5, 0, 0.5];

line([nanmedian(FingerprintAUC(1, :)), nanmedian(FingerprintAUC(1, :))], get(gca, 'ylim'), 'Color', [0.25, 0.25, 0])
line([nanmedian(FingerprintAUC(2, :)), nanmedian(FingerprintAUC(2, :))], get(gca, 'ylim'), 'Color', [0, 0.5, 0.5])
line([nanmedian(FingerprintAUC(3, :)), nanmedian(FingerprintAUC(3, :))], get(gca, 'ylim'), 'Color', [0.5, 0, 0.5])

xlabel('AUC cross-correlation')
ylabel('Nr sessions')
makepretty

%%
subplot(3, 3, 6)
hold on
h1 = histogram(ACGAUC(1, :), edges);
h1.FaceColor = [0.25, 0.25, 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25, 0.25, 0];

h2 = histogram(ACGAUC(2, :), edges);
h2.FaceColor = [0, 0.5, 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0, 0.5, 0.5];

h3 = histogram(ACGAUC(3, :), edges);
h3.FaceColor = [0.5, 0, 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5, 0, 0.5];

line([nanmedian(ACGAUC(1, :)), nanmedian(ACGAUC(1, :))], get(gca, 'ylim'), 'Color', [0.25, 0.25, 0])
line([nanmedian(ACGAUC(2, :)), nanmedian(ACGAUC(2, :))], get(gca, 'ylim'), 'Color', [0, 0.5, 0.5])
line([nanmedian(ACGAUC(3, :)), nanmedian(ACGAUC(3, :))], get(gca, 'ylim'), 'Color', [0.5, 0, 0.5])

xlabel('AUC ACG correlation')
ylabel('Nr sessions')
makepretty

%%
subplot(3, 3, 9)
hold on
h1 = histogram(RFAUC(1, :), edges);
h1.FaceColor = [0.25, 0.25, 0];
% h1.FaceAlpha = 0.5;
h1.EdgeColor = [0.25, 0.25, 0];


h2 = histogram(RFAUC(2, :), edges);
h2.FaceColor = [0, 0.5, 0.5];
% h2.FaceAlpha = 0.5;
h2.EdgeColor = [0, 0.5, 0.5];

h3 = histogram(RFAUC(3, :), edges);
h3.FaceColor = [0.5, 0, 0.5];
% h2.FaceAlpha = 0.5;
h3.EdgeColor = [0.5, 0, 0.5];
line([nanmedian(RFAUC(1, :)), nanmedian(RFAUC(1, :))], get(gca, 'ylim'), 'Color', [0.25, 0.25, 0])
line([nanmedian(RFAUC(2, :)), nanmedian(RFAUC(2, :))], get(gca, 'ylim'), 'Color', [0, 0.5, 0.5])
line([nanmedian(RFAUC(3, :)), nanmedian(RFAUC(3, :))], get(gca, 'ylim'), 'Color', [0.5, 0, 0.5])

xlabel('AUC RF distance')
ylabel('Nr sessions')
makepretty
saveas(FSCoreFig, fullfile(SaveDir, 'FunctionScoreFigAcrossMice.fig'))
saveas(FSCoreFig, fullfile(SaveDir, 'FunctionScoreFigAcrossMice.bmp'))

%% Extra Positive/negative
EPosAndNeg(:, sum(isnan(EPosAndNeg), 1) == 3) = [];
EPosAndNeg(sum(isnan(EPosAndNeg), 2) == size(EPosAndNeg, 2), :) = [];

figure('name', 'ExtraPositiveAndNegative');
subplot(1,3,1)
bar(EPosAndNeg'.*100, 'EdgeColor', 'none', 'BarWidth', 1)
set(gca, 'XTick', 1:length(MiceOpt), 'XTickLabel', MiceOpt, 'YAxisLocation', 'right', 'XTickLabelRotation', 90)
ylabel('% (relative to spike sorting)')
legend('Extra Positives', 'Extra Negatives Within', 'Extra Negatives Across')
makepretty

subplot(1,3,2)
cols = lines(size(EPosAndNeg,2));
hold on
scatter(EPosAndNeg(1,:).*100,EPosAndNeg(2,:).*100,40,cols,'filled')

% xlims = get(gca,'xlim');
% ylims = get(gca,'ylim');
line([0 8],[0 8],'color',[0 0 0],'LineStyle','-')

title('Within')
xlabel('Units to be merged (%)')
ylabel('Units to be split (%)')
makepretty

if size(EPosAndNeg,1)>2
    subplot(1,3,3)
    scatter(EPosAndNeg(1,:).*100,EPosAndNeg(3,:).*100,40,cols,'filled')

    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    line([0 max([xlims ylims])],[0 max([xlims ylims])],'color',[0 0 0],'LineStyle','-')
    title('Across')
    xlabel('Units to be merged (%)')
    ylabel('Units to be split (%)')
    makepretty
end


nanmean(EPosAndNeg.*100,2)
nanstd(EPosAndNeg.*100,[],2)

%% Tracking performance
%             TrackingPerformance = cat(2,TrackingPerformance,[did2-did1,nMatches,nMax]');

% UMTrackingPerformancePerMouse{midx} = TrackingPerformance;
%     KSTrackingPerformancePerMouse{midx} = TrackingPerformanceKS;
tmpUM = cat(2, UMTrackingPerformancePerMouse{:});
tmpUM(:,tmpUM(4,:)<15) = [];
if UseKSLabels
    tmpKS = cat(2, KSTrackingPerformancePerMouse{:});
end
figure('name', 'Tracking Performance when tracking is expected')
% scatter(1:size(tmpUM,2), tmpUM(3, :)./tmpUM(3, :), 20, [0, 0, 0], 'filled')
hold on
if UseKSLabels
    scatter(1:sum(tmpUM(1,:)==1), tmpKS(3, tmpUM(1,:)==1)./tmpKS(4, tmpUM(1,:)==1).*100, 20, tmpUM(2,tmpUM(1,:)==1))
end
scatter(1:sum(tmpUM(1,:)==1), tmpUM(3, tmpUM(1,:)==1)./tmpUM(4,  tmpUM(1,:)==1).*100, 20, tmpUM(2, tmpUM(1,:)==1), 'filled')
line([0.5 sum(tmpUM(1,:)==1)+0.5],[100 100],'color',[0 0 0],'LineStyle','-')
% set(gca, 'XTick', 1:length(tmpUM), 'XTickLabel', 1:size(tmpUM,2), 'XTickLabelRotation', 90)
if UseKSLabels
    legend('Kilosort tracked', 'UnitMatch tracked (concatenated)','maximum possible')
else
    legend('UnitMatch tracked (not concatenated)','maximum possible')
end
xlim([0.5,sum(tmpUM(1,:)==1) + 0.5])
xlabel('SessionID')
ylabel('Units Traced (%)')
hc = colorbar;
hc.Label.String = '\Delta Recording days'
colormap(cool)
makepretty

saveas(gcf,fullfile(SaveDir,'TrackingPerformance.fig'))
saveas(gcf,fullfile(SaveDir,'TrackingPerformance.bmp'))

%% Split per type of recordings
NonEmptyIdx = ~cellfun(@isempty,UMTrackingPerformancePerMouse);
RecTypeOpt = {'Acute','Some IMRO differences','Chronic - same IMRO'}
UMTrackingPerformancePerMouse = UMTrackingPerformancePerMouse(NonEmptyIdx);
tmpmice = MiceOpt(NonEmptyIdx);
rectypeid = nan(1,size(tmpUM,2));
rectypeid(find(tmpUM(1,:) == 0)) = 1; % Acute
rectypeid(find(tmpUM(1,:) > 0 & tmpUM(1,:) <1 )) = 2; % Some IMRO differences
rectypeid(find(tmpUM(1,:) == 1)) = 3; %Chronic, same IMRO
OneDayDiffIdx = tmpUM(2,:)==1;

AllDat = tmpUM(3,OneDayDiffIdx)./tmpUM(4,OneDayDiffIdx).*100;
stepsz = 10;
Edges = floor(min(AllDat)):stepsz:round(max(AllDat)+stepsz);
histvec = floor(min(AllDat))+stepsz/2:stepsz:round(max(AllDat)+stepsz)-stepsz/2;

cols = [0 0 1; 0 1 0; 1 0 0];
figure('name','TrackedUnits (%)')
clear h
hold on
for tid = 1:length(RecTypeOpt)
    hc = histcounts(AllDat(rectypeid(OneDayDiffIdx)==tid),Edges)./sum(rectypeid(OneDayDiffIdx)==tid);
    h(tid) = plot(histvec,hc.*100,'color',cols(tid,:));
%     h(tid).FaceColor = cols(tid,:);
%     h(tid).EdgeColor = cols(tid,:);
%     h(tid).FaceAlpha = 0.5;
end
xlabel('Tracked units (%)')
ylabel('Pairs of days (%)')
legend(RecTypeOpt)
makepretty


%% Across days
figure('name','Tracking across days')
cols = lines(length(UMTrackingPerformancePerMouse));
for tid = 1:length(RecTypeOpt)
    subplot(1,length(RecTypeOpt),tid)
    hold on
    for midx = 1:length(UMTrackingPerformancePerMouse)
        if ~isempty(UMTrackingPerformancePerMouse{midx})
            if tid == 1
                Idx = UMTrackingPerformancePerMouse{midx}(1,:) == 0;
            elseif tid == 2
                Idx = UMTrackingPerformancePerMouse{midx}(1,:)~=0 & UMTrackingPerformancePerMouse{midx}(1,:)~=1;
            elseif tid == 3
                Idx = UMTrackingPerformancePerMouse{midx}(1,:) == 1;
            end
            scatter(UMTrackingPerformancePerMouse{midx}(2,Idx),UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100,20,cols(midx,:),'filled')
        end
    end
    ylabel('Tracked units (%)')
    xlabel('\Delta Days')
    ylim([0 100])
    title(RecTypeOpt{tid})
    
    makepretty
end
legend(tmpmice)

%% Tracking For chronic only
expFun = @(p,d) p(1)*exp(-p(2)*d);%+p(3); % For nneurons decay
opts = optimset('Display','off');

figure('name','ChronicTrackingAcrossDays')
cols = distinguishable_colors(length(UMTrackingPerformancePerMouse));
miceincluded = true(1,length(UMTrackingPerformancePerMouse));
for midx = 1:length(UMTrackingPerformancePerMouse)
    subplot(ceil(sqrt(length(UMTrackingPerformancePerMouse)+1)),round(sqrt(length(UMTrackingPerformancePerMouse)+1)),midx)
    Idx = UMTrackingPerformancePerMouse{midx}(1,:) == 1;
    if ~any(Idx)
        miceincluded(midx)=0;
        continue
    end
    hold on
    scatter(UMTrackingPerformancePerMouse{midx}(2,Idx),UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100,20,cols(midx,:),'filled')
    p = lsqcurvefit(expFun,[max(UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100) 0.05],UMTrackingPerformancePerMouse{midx}(2,Idx),UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100,[],[],opts);
    plot(unique(UMTrackingPerformancePerMouse{midx}(2,Idx)),expFun(p,unique(UMTrackingPerformancePerMouse{midx}(2,Idx))),'color',cols(midx,:))
    title(tmpmice{midx})
    xlabel('\Deltadays')
    ylabel('Tracked (%)')
    makepretty
end
subplot(ceil(sqrt(length(UMTrackingPerformancePerMouse)+1)),round(sqrt(length(UMTrackingPerformancePerMouse)+1)),length(UMTrackingPerformancePerMouse)+1)

for midx = 1:length(UMTrackingPerformancePerMouse)
    Idx = UMTrackingPerformancePerMouse{midx}(1,:) == 1;
    if ~any(Idx)
        miceincluded(midx)=0;
        continue
    end
    hold on
    p = lsqcurvefit(expFun,[max(UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100) 0.05],UMTrackingPerformancePerMouse{midx}(2,Idx),UMTrackingPerformancePerMouse{midx}(3,Idx)./UMTrackingPerformancePerMouse{midx}(4,Idx).*100,[],[],opts);
    plot(unique(UMTrackingPerformancePerMouse{midx}(2,Idx)),expFun(p,unique(UMTrackingPerformancePerMouse{midx}(2,Idx))),'color',cols(midx,:))

end
title('All chronic mice')
xlabel('\Deltadays')
ylabel('Tracked (%)')
makepretty
legend(tmpmice(miceincluded))



%% save figure
figure(MatrixFig)
saveas(gcf,fullfile(SaveDir,'TrackingPerformanceMatrix.fig'))
saveas(gcf,fullfile(SaveDir,'TrackingPerformanceMatrix.bmp'))

