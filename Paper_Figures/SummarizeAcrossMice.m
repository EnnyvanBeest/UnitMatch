% SummarizeAcrossMice
FSCoreFig = figure('name', 'Functional Scores');

% Initialize
TakeRank = 0 %if 0 , take cross-correlation scores (these may be less informative than rank)
if TakeRank
    stepsz = 1;
    bins = -20:stepsz:-1;
    Vector =  -20+stepsz/2:stepsz:-1-stepsz/2;
else
    stepsz = 0.1;
    bins = -1:stepsz:1;
    Vector = -1+stepsz/2:stepsz:1-stepsz/2;
end
FingerPrintRAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
FingerprintAUC = nan(3, length(MiceOpt));
ACGRAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
ACGAUC = nan(3, length(MiceOpt));
RFDistAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
RFAUC = nan(3, length(MiceOpt));
NatImgAcrossMice = nan(length(Vector), 3, length(MiceOpt)); % Vector / within,Match,Non-match / mice
NatImgAUC = nan(3, length(MiceOpt));

EPosAndNeg = nan(3, length(MiceOpt));
UseKSLabels = PipelineParams.RunPyKSChronicStitched;
RecSesPerUIDAllMice = cell(1,length(MiceOpt));

% Prepare these
AllFPCor = [];
AllCGCor = [];
AllRFCor = [];
AllNMCor = [];

AUCFPCurves = [];
AUCACGCurves = [];
AUCNImCurves = [];
AUCRFCurves = [];

minMatches = 10; % if less than 10, skip

AUCCols = [0 0.7 0; 1 0 0; 0 0 0.7]; %WIthin %Match %non-match
aucprecision = [0:0.05:1];
clear UMTrackingPerformancePerMouse
for midx = 1:length(MiceOpt)
    fprintf('Reference %s...\n', MiceOpt{midx})

    tmpfile = dir(fullfile(SaveDir, MiceOpt{midx}, '*','*', 'UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end
    for tmpid = 1:length(tmpfile)

        fprintf('Loading the data...\n')
        tic
        load(fullfile(tmpfile(tmpid).folder, tmpfile(tmpid).name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
        toc

        % Load AUCS
        if exist(fullfile(tmpfile(tmpid).folder, 'AUC.mat'))
            AUC = load(fullfile(tmpfile(tmpid).folder, 'AUC.mat'))';
            if (midx == 1 && tmpid==1)|| ~exist('AUCParams','var')
                AUCParams = AUC.AUCStruct.ParamNames;
                AUCVals = nan(length(AUCParams), 0);
            end
            AUCVals = cat(2,AUCVals,AUC.AUCStruct.AUC');
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
        if numel(MatchIdx)<minMatches
            continue
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
        if isstruct(AllRawDir{1})
            AllRawDir = cellfun(@(X) fullfile(X.folder,X.name),AllRawDir,'Uni',0);
        end
        nclus = length(UniqueID);
        ChannelPos = UMparam.AllChannelPos;

        %% How many units vs number of units were tracked?
        fprintf('Evaluating tracked units...\n')
        CrossCorrMatching = nan(0,4); % number of units matched with cross-correlation, % number of units matched of those % Number of units match with UM, % number of units mof those matched with Cross Corr
        TrackingPerformance = nan(5, 0); % MatchesExpected, Difference between recording day, % Tracked units %maximum possibility % Difference in number of units
        TrackingPerformanceKS = nan(5, 0); % MatchesExpected, Difference between recording day, % Tracked units %maximum possibility
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


                % Units found to match based on
                %cross-correlation activity
                rowidx = find(MatchTable.RecSes1 == did1 & MatchTable.RecSes2 == did2);
                % find all pairs with a significant fingerprint cross-correlation
                SigR = find(MatchTable.refPopCorr(rowidx) == 1);
                NXCorr = length(SigR);
                NXCorrUM = sum(MatchTable.MatchProb(rowidx(SigR))>0.5);

                % OTher way around
                SigUM = find(MatchTable.MatchProb(rowidx) > 0.5);
                NUM = length(SigUM);
                NUMXCorr = sum(MatchTable.refPopCorr(rowidx(SigUM))==1);

                CrossCorrMatching = cat(1, CrossCorrMatching, [NXCorr, NXCorrUM, NUM, NUMXCorr]);

                thesedaysidx = find(ismember(recses, [did1, did2]));
                % can possibly only track these many units:
                nMax = min([sum(recses(thesedaysidx) == did1), sum(recses(thesedaysidx) == did2)]);
                diffN = sum(recses(thesedaysidx) == did2)-sum(recses(thesedaysidx) == did1);
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
                TrackingPerformance = cat(2, TrackingPerformance, [MatchesExpected, double(DDay), nMatches, nMax, diffN]');
                if MatchesExpected == 1 && (nMatches/nMax)<0.15 && double(DDay)<2
                    warning(['No matches unexpectedly for ' MiceOpt{midx}])
                end
                if MatchesExpected == 0 && (nMatches/nMax)>0.15
                    warning(['Some matches unexpectedly for ' MiceOpt{midx}])
                end

                if sum(recses(thesedaysidx) == did1)>sum(recses(thesedaysidx) == did2)
                    nMatches = sum(ismember(OriID(recses==did2),OriID(recses==did1)));
                else
                    nMatches = sum(ismember(OriID(recses==did1),OriID(recses==did2)));
                end
                if UseKSLabels
                    TrackingPerformanceKS = cat(2, TrackingPerformanceKS, [MatchesExpected, double(DDay), nMatches, nMax, diffN]');
                end
            end
        end

        UMTrackingPerformancePerMouse{midx,tmpid} = TrackingPerformance;
        if UseKSLabels
            KSTrackingPerformancePerMouse{midx,tmpid} = TrackingPerformanceKS;
        end

        %% UniqueID X tracking
        [UniqueIDOpt,idx1,idx2] = unique(UniqueID); %UID options
        DayOpt = unique(AllDeltaDays);
        dayses = AllDeltaDays(recses); % Get in days after first recording
        RecSesPerUID = arrayfun(@(X) (ismember(DayOpt,dayses(idx2==X))),1:numel(UniqueIDOpt),'Uni',0); % Extract which recording sessions a unite appears in
        RecSesPerUID = double(cat(1,RecSesPerUID{:}));
        for nrec = 1:size(RecSesPerUID,2)
            RecSesPerUID(RecSesPerUID(:,nrec)>0,nrec) = RecSesPerUID(RecSesPerUID(:,nrec)>0,nrec)+(size(RecSesPerUID,2)-nrec)+1;
        end
        if tmpid==1
            MatrixFig = figure('name',[MiceOpt{midx} ' TrackingAcrossDaysMatrix']);
        else
            figure(MatrixFig)
        end
        %     [sortval,idx1,sortidx] = unique(RecSesPerUID,'rows');
        [sortval,sortidx] = sort(sum(RecSesPerUID,2),'descend');
        subplot(ceil(sqrt(length(tmpfile))),round(sqrt(length(tmpfile))),tmpid)
        h = imagesc(1:length(DayOpt),[],RecSesPerUID(sortidx,:));
        set(gca,'XTick',1:length(DayOpt),'XTickLabel',DayOpt)
        colormap(flipud(gray))
        title(tmpfile(tmpid).folder)
        xlabel('Days since first recording')
        ylabel('Tracked Units')
        makepretty
        RecSesPerUIDAllMice{midx,tmpid} = RecSesPerUID; % Save for later use

        %% Extra Positives/negatives
        MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
        %Extra Positive
        EPosAndNeg(1, midx, tmpid) = sum(MatchProb(NonMatchIdx) > UMparam.ProbabilityThreshold) ./ length(NonMatchIdx);

        %Extra Negative
        EPosAndNeg(2, midx, tmpid) = sum(MatchProb(WithinIdx) < UMparam.ProbabilityThreshold) ./ length(WithinIdx);

        if UseKSLabels
            EPosAndNeg(3, midx, tmpid) = sum(MatchProb(MatchIdx) < UMparam.ProbabilityThreshold) ./ length(MatchIdx);
        end

        %% Fingerprint correlation
        fprintf('Plotting fingerprints...\n')
        if any(ismember(MatchTable.Properties.VariableNames, 'refPopCorr')) && MatchesExpected == 1
            figure(FSCoreFig)
            if TakeRank
                FingerprintCor = -reshape(MatchTable.refPopRank, nclus, nclus);
            else
                FingerprintCor = reshape(MatchTable.refPopCorr, nclus, nclus);
            end

            subplot(4, 3, 1)
            hold on

            hw = histcounts(FingerprintCor(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(FingerprintCor(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(FingerprintCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
            FingerPrintRAcrossMice(:, :, midx,tmpid) = cat(1, hw, hm, hn)';
            plot(Vector, hw, 'color', AUCCols(1,:))
            plot(Vector, hm, 'color', AUCCols(2,:))
            plot(Vector, hn, 'color', AUCCols(3,:))
            xlabel('Cross-correlation Fingerprint')
            ylabel('Proportion|Group')
            %     legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
            axis square
            makepretty

            AllFPCor = cat(3,AllFPCor,[hw;hm;hn]);
            subplot(4, 3, 2)
            hold on
            clear h
            if length(MatchIdx) > 10
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [FingerprintCor(MatchIdx)', FingerprintCor(NonMatchIdx)'];
                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(1) = plot(X, Y, 'color', nanmean(AUCCols([2,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y1 = Y(idx1);
                hold all
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
                scores = [FingerprintCor(MatchIdx)', FingerprintCor(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(2) = plot(X, Y, 'color',nanmean(AUCCols([1,2],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y2 = Y(idx1);
                labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [FingerprintCor(WithinIdx)', FingerprintCor(NonMatchIdx)'];
                [X, Y, ~, AUC3] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(3) = plot(X, Y, 'color', nanmean(AUCCols([1,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y3 = Y(idx1);
                axis square
            else
                AUC1 = nan;
                AUC2 = nan;
                AUC3 = nan;
                Y1 = nan(numel(aucprecision),1);
                Y2 = nan(numel(aucprecision),1);
                Y3 = nan(numel(aucprecision),1);
            end
            FingerprintAUC(:, midx,tmpid) = [AUC1, AUC2, AUC3];
            AUCFPCurves = cat(3,AUCFPCurves,[Y1,Y2,Y3]);


            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')




            %     if length(MatchIdx) > 10
            %         legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
            %     else
            %         legend([h(3)], 'Within vs No Match', 'Location', 'best')
            %     end
            title('Cross-Correlation Fingerprint')
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing
        end
        %% Autocorrelogram
        if any(ismember(MatchTable.Properties.VariableNames, 'ACGCorr')) && MatchesExpected == 1
            if TakeRank
                ACGCor = -reshape(MatchTable.ACGRank, nclus, nclus);
            else
                ACGCor = reshape(MatchTable.ACGCorr, nclus, nclus);
            end

            subplot(4, 3, 4)
            hold on

            hw = histcounts(ACGCor(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(ACGCor(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(ACGCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
            ACGRAcrossMice(:, :, midx,tmpid) = cat(1, hw, hm, hn)';
            plot(Vector, hw, 'color', AUCCols(1,:))
            plot(Vector, hm, 'color', AUCCols(2,:))
            plot(Vector, hn, 'color', AUCCols(3,:))
            xlabel('Autocorrelogram')
            ylabel('Proportion|Group')
            axis square
            makepretty

            AllCGCor = cat(3,AllCGCor,[hw;hm;hn]);


            subplot(4, 3, 5)
            hold on
            clear h
            if ~isempty(MatchIdx)
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [ACGCor(MatchIdx)', ACGCor(NonMatchIdx)'];
                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(1) = plot(X, Y, 'color', nanmean(AUCCols([2,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y1 = Y(idx1);
                hold all
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
                scores = [ACGCor(MatchIdx)', ACGCor(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(2) = plot(X, Y, 'color', nanmean(AUCCols([1,2],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y2 = Y(idx1);
                labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [ACGCor(WithinIdx)', ACGCor(NonMatchIdx)'];
                [X, Y, ~, AUC3] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(3) = plot(X, Y, 'color', nanmean(AUCCols([1,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y3 = Y(idx1);
                axis square
            else
                AUC1 = nan;
                AUC2 = nan;
                AUC3 = nan;
                Y1 = nan(numel(aucprecision),1);
                Y2 = nan(numel(aucprecision),1);
                Y3 = nan(numel(aucprecision),1);
            end
            ACGAUC(:, midx,tmpid) = [AUC1, AUC2, AUC3];
            AUCACGCurves = cat(3,AUCACGCurves,[Y1,Y2,Y3]);

            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')
            title('Auto-correlogram')
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing



        end
        %% Natural Images correlation
        if any(ismember(MatchTable.Properties.VariableNames, 'natImCorr')) && MatchesExpected == 1
            figure(FSCoreFig)
            if TakeRank
                NatImCorr = -reshape(MatchTable.natImRank, nclus, nclus);
            else
                NatImCorr = reshape(MatchTable.natImCorr, nclus, nclus);
            end


            subplot(4, 3, 7)
            hold on

            hw = histcounts(NatImCorr(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(NatImCorr(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(NatImCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
            NatImgAcrossMice(:, :, midx,tmpid) = cat(1, hw, hm, hn)';
            plot(Vector, hw, 'color', AUCCols(1,:))
            plot(Vector, hm, 'color', AUCCols(2,:))
            plot(Vector, hn, 'color', AUCCols(3,:))
            xlabel('Natural Images Fingerprint')
            ylabel('Proportion|Group')
            %     legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
            axis square
            makepretty

            AllNMCor = cat(3,AllNMCor,[hw;hm;hn]);

            subplot(4, 3, 8)
            hold on
            clear h
            if ~isempty(MatchIdx) && ~all(isnan(NatImCorr(MatchIdx)))
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [NatImCorr(MatchIdx)', NatImCorr(NonMatchIdx)'];
                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(1) = plot(X, Y, 'color', nanmean(AUCCols([2,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y1 = Y(idx1);
                hold all
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
                scores = [NatImCorr(MatchIdx)', NatImCorr(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(2) = plot(X, Y, 'color',nanmean(AUCCols([1,2],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y2 = Y(idx1);
                labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [NatImCorr(WithinIdx)', NatImCorr(NonMatchIdx)'];
                [X, Y, ~, AUC3] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(3) = plot(X, Y, 'color', nanmean(AUCCols([1,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y3 = Y(idx1);
                axis square
            else
                AUC1 = nan;
                AUC2 = nan;
                AUC3 = nan;
                Y1 = nan(numel(aucprecision),1);
                Y2 = nan(numel(aucprecision),1);
                Y3 = nan(numel(aucprecision),1);
            end
            NatImgAUC(:, midx,tmpid) = [AUC1, AUC2, AUC3];
            AUCNImCurves = cat(3,AUCNImCurves,[Y1,Y2,Y3]);

            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')
            %     if length(MatchIdx) > 10
            %         legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
            %     else
            %         legend([h(3)], 'Within vs No Match', 'Location', 'best')
            %     end
            title('Natural Images Fingerprint')
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing


        end
        %% Receptive Field (?)
        if any(ismember(MatchTable.Properties.VariableNames, 'RFDist')) && MatchesExpected == 1

            RFDist = reshape(MatchTable.RFDist, nclus, nclus);
            subplot(4, 3, 10)
            hold on

            bins = linspace(0, 50, length(Vector)+1);
            stepsz = unique(diff(bins));
            hw = histcounts(RFDist(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(RFDist(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(RFDist(NonMatchIdx), bins) ./ length(NonMatchIdx);
            RFDistAcrossMice(:, :, midx,tmpid) = cat(1, hw, hm, hn)';
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hw, 'color', AUCCols(1,:))
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hm, 'color', AUCCols(2,:))
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hn, 'color', AUCCols(3,:))
            xlabel('Receptive Field Distance')
            ylabel('Proportion|Group')
            axis square
            makepretty
            AllRFCor = cat(3,AllRFCor,[hw;hm;hn]);


            subplot(4, 3, 11)
            hold on
            clear h
            if ~isempty(MatchIdx)
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(NonMatchIdx))];
                scores = [RFDist(MatchIdx)', RFDist(NonMatchIdx)'];
                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(1) = plot(X, Y, 'color', nanmean(AUCCols([2,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y1 = Y(idx1);
                hold all
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
                scores = [RFDist(MatchIdx)', RFDist(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(2) = plot(X, Y, 'color', nanmean(AUCCols([1,2],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y2 = Y(idx1);
                labels = [zeros(1, numel(WithinIdx)), ones(1, numel(NonMatchIdx))];
                scores = [RFDist(WithinIdx)', RFDist(NonMatchIdx)'];
                [X, Y, ~, AUC3] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(3) = plot(X, Y, 'color', nanmean(AUCCols([1,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y3 = Y(idx1);
                axis square
            else
                AUC1 = nan;
                AUC2 = nan;
                AUC3 = nan;
                Y1 = nan(numel(aucprecision),1);
                Y2 = nan(numel(aucprecision),1);
                Y3 = nan(numel(aucprecision),1);
            end
            RFAUC(:, midx,tmpid) = [AUC1, AUC2, AUC3];
            AUCRFCurves = cat(3,AUCRFCurves,[Y1,Y2,Y3]);

            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')
            title('RF Distance')
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing


        elseif any(ismember(MatchTable.Properties.VariableNames, 'FRDiff')) && MatchesExpected == 1
            if TakeRank
                RFDist = -reshape(MatchTable.FRRank, nclus, nclus);
            else
                RFDist = reshape(MatchTable.FRDiff, nclus, nclus);
                RFDist = -(RFDist-nanmin(RFDist(:)))./(nanmax(RFDist(:))-nanmin(RFDist(:))); %normalize between 0 and 1 to resemble correlation
            end
            subplot(4, 3, 10)
            hold on

            hw = histcounts(RFDist(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(RFDist(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(RFDist(NonMatchIdx), bins) ./ length(NonMatchIdx);
            RFDistAcrossMice(:, :, midx,tmpid) = cat(1, hw, hm, hn)';
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hw, 'color', AUCCols(1,:))
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hm, 'color', AUCCols(2,:))
            plot(bins(1)+stepsz/2:stepsz:bins(end)-stepsz/2, hn, 'color', AUCCols(3,:))
            xlabel('Firing Rate Distance')
            ylabel('Proportion|Group')
            axis square
            makepretty
            AllRFCor = cat(3,AllRFCor,[hw;hm;hn]);


            subplot(4, 3, 11)
            hold on
            clear h
            if ~isempty(MatchIdx)
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [RFDist(MatchIdx)', RFDist(NonMatchIdx)'];
                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(1) = plot(X, Y, 'color', nanmean(AUCCols([2,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y1 = Y(idx1);
                hold all
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(WithinIdx))];
                scores = [RFDist(MatchIdx)', RFDist(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(2) = plot(X, Y, 'color', nanmean(AUCCols([1,2],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y2 = Y(idx1);
                labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [RFDist(WithinIdx)', RFDist(NonMatchIdx)'];
                [X, Y, ~, AUC3] = perfcurve(labels, scores, 1,'XVals',aucprecision);
                h(3) = plot(X, Y, 'color', nanmean(AUCCols([1,3],:),1));
                % find closest precision value for every X
                [~,idx1] = min(abs(aucprecision-X),[],1);
                Y3 = Y(idx1);
                axis square
            else
                AUC1 = nan;
                AUC2 = nan;
                AUC3 = nan;
                Y1 = nan(numel(aucprecision),1);
                Y2 = nan(numel(aucprecision),1);
                Y3 = nan(numel(aucprecision),1);
            end
            RFAUC(:, midx,tmpid) = [AUC1, AUC2, AUC3];
            AUCRFCurves = cat(3,AUCRFCurves,[Y1,Y2,Y3]);

            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')
            title('Firing rate Difference')
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing


        end

        fprintf('Done!\n')
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
hc.Label.String = '\Delta Recording days';
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

%% Split per type of recordings
RecTypeOpt = {'Different IMRO tables','Chronic - same IMRO'}
rectypeid = nan(1,size(tmpUM,2));
rectypeid(find(tmpUM(1,:) == 0)) = 1; % Acute
% rectypeid(find(tmpUM(1,:) > 0 & tmpUM(1,:) <1 )) = 2; % Some IMRO differences
rectypeid(find(tmpUM(1,:) == 1)) = 2; %Chronic, same IMRO
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
ylim([0 100])

% subplot(ceil(sqrt(length(UMTrackingPerformancePerMouse)+1)),round(sqrt(length(UMTrackingPerformancePerMouse)+1)),length(UMTrackingPerformancePerMouse)+1)

figure('name','All Mice Tracking peformance')
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
ylim([0 100])



%% save figure
figure(MatrixFig)
saveas(gcf,fullfile(SaveDir,'TrackingPerformanceMatrix.fig'))
saveas(gcf,fullfile(SaveDir,'TrackingPerformanceMatrix.bmp'))

%% Average figure
figure('name','FunctionalScoresAverages')
subplot(4,3,1)
clear h
hold on
for hid = 1:3
    if numel(size(AllFPCor))>2
        h(hid) = shadedErrorBar(Vector, squeeze(nanmean(AllFPCor(hid,:,:),3)),squeeze(nanstd(AllFPCor(hid,:,:),[],3))./sqrt(size(AllFPCor,3)-1));

        h(hid).mainLine.Color = AUCCols(hid,:);
        h(hid).patch.FaceColor = AUCCols(hid,:);
        h(hid).edge(1).Color = AUCCols(hid,:);
        h(hid).edge(2).Color = AUCCols(hid,:);
    else
        h(hid) = plot(Vector, squeeze(nanmean(AllFPCor(hid,:,:),3)));
        h(hid).Color = AUCCols(hid,:);

    end

end
xlabel('Cross-correlation rank')
ylabel('Proportion|Group')
if numel(size(AllFPCor))>2

    legend([h(:).mainLine], 'Same unit within', 'Matches', 'Neighbors', 'Location', 'best')
else
    legend([h(:)], 'Same unit within', 'Matches', 'Neighbors', 'Location', 'best')

end
makepretty

subplot(4,3,4)
clear h
hold on
for hid = 1:3
    if numel(size(AllFPCor))>2

        h(hid) = shadedErrorBar(Vector, squeeze(nanmean(AllCGCor(hid,:,:),3)),squeeze(nanstd(AllCGCor(hid,:,:),[],3))./sqrt(size(AllCGCor,3)-1));
        h(hid).mainLine.Color = AUCCols(hid,:);
        h(hid).patch.FaceColor = AUCCols(hid,:);
        h(hid).edge(1).Color = AUCCols(hid,:);
        h(hid).edge(2).Color = AUCCols(hid,:);
    else
        h(hid) = plot(Vector, squeeze(nanmean(AllCGCor(hid,:,:),3)));
        h(hid).Color = AUCCols(hid,:);

    end
end
xlabel('ACG Correlation')
ylabel('Proportion|Group')
makepretty

if any(~isnan(AllNMCor(:)))
    subplot(4,3,7)
    clear h
    hold on
    for hid = 1:3
        if numel(size(AllFPCor))>2

            h(hid) = shadedErrorBar(Vector, squeeze(nanmean(AllNMCor(hid,:,:),3)),squeeze(nanstd(AllNMCor(hid,:,:),[],3))./sqrt(size(AllNMCor,3)-1));
            h(hid).mainLine.Color = AUCCols(hid,:);
            h(hid).patch.FaceColor = AUCCols(hid,:);
            h(hid).edge(1).Color = AUCCols(hid,:);
            h(hid).edge(2).Color = AUCCols(hid,:);
        else
            h(hid) = plot(Vector, squeeze(nanmean(AllNMCor(hid,:,:),3)));
            h(hid).Color = AUCCols(hid,:);

        end
    end
    xlabel('Natural Images Fingerprint')
    ylabel('Proportion|Group')
    makepretty
end

% AllRFCor
if any(~isnan(AllRFCor(:)))
    subplot(4,3,10)
    clear h
    hold on
    for hid = 1:3
        if numel(size(AllFPCor))>2

            h(hid) = shadedErrorBar(Vector, squeeze(nanmean(AllRFCor(hid,:,:),3)),squeeze(nanstd(AllRFCor(hid,:,:),[],3))./sqrt(size(AllRFCor,3)-1));
            h(hid).mainLine.Color = AUCCols(hid,:);
            h(hid).patch.FaceColor = AUCCols(hid,:);
            h(hid).edge(1).Color = AUCCols(hid,:);
            h(hid).edge(2).Color = AUCCols(hid,:);
        else
            h(hid) = plot(Vector, squeeze(nanmean(AllRFCor(hid,:,:),3)));
            h(hid).Color = AUCCols(hid,:);
        end
    end
    xlabel('RF Distance')
    ylabel('Proportion|Group')
    makepretty
end

%% AUC curves
% AUCFPCurves = [];
% AUCACGCurves = [];
% AUCNImCurves = [];
% AUCRFCurves = [];
% % aucprecision = [0:0.01:1];
AUCmixedCols = cat(1,nanmean(AUCCols([2,3],:),1),nanmean(AUCCols([2,1],:),1), nanmean(AUCCols([1,3],:),1));

subplot(4,3,2)
clear h
hold on
for hid = 1:3
    if numel(size(AUCFPCurves))>2

        h(hid) = shadedErrorBar(aucprecision, squeeze(nanmean(AUCFPCurves(:,hid,:),3)),squeeze(nanstd(AUCFPCurves(:,hid,:),[],3))./sqrt(size(AUCFPCurves,3)-1));
        h(hid).mainLine.Color = AUCmixedCols(hid,:);
        h(hid).patch.FaceColor = AUCmixedCols(hid,:);
        h(hid).edge(1).Color = AUCmixedCols(hid,:);
        h(hid).edge(2).Color = AUCmixedCols(hid,:);
    else
        h(hid) = plot(aucprecision, squeeze(nanmean(AUCFPCurves(:,hid,:),3)));
        h(hid).Color = AUCmixedCols(hid,:);
    end
end
plot([0, 1], [0, 1], 'k--')

xlabel('False positives')
ylabel('True positives')
makepretty
if numel(size(AUCFPCurves))>2
    legend([h(:).mainLine], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
else
    legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')

end

subplot(4,3,5)
clear h
hold on
for hid = 1:3
    if numel(size(AllFPCor))>2

        h(hid) = shadedErrorBar(aucprecision, squeeze(nanmean(AUCACGCurves(:,hid,:),3)),squeeze(nanstd(AUCACGCurves(:,hid,:),[],3))./sqrt(size(AUCACGCurves,3)-1));
        h(hid).mainLine.Color = AUCmixedCols(hid,:);
        h(hid).patch.FaceColor = AUCmixedCols(hid,:);
        h(hid).edge(1).Color = AUCmixedCols(hid,:);
        h(hid).edge(2).Color = AUCmixedCols(hid,:);
    else
        h(hid) = plot(aucprecision, squeeze(nanmean(AUCACGCurves(:,hid,:),3)));
        h(hid).Color = AUCmixedCols(hid,:);
    end

end
plot([0, 1], [0, 1], 'k--')

xlabel('False positives')
ylabel('True positives')
makepretty

if any(~isnan(NatImgAUC(:)))

    subplot(4,3,8)
    clear h
    hold on
    for hid = 1:3
        if numel(size(AllFPCor))>2

            h(hid) = shadedErrorBar(aucprecision, squeeze(nanmean(AUCNImCurves(:,hid,:),3)),squeeze(nanstd(AUCNImCurves(:,hid,:),[],3))./sqrt(size(AUCNImCurves,3)-1));
            h(hid).mainLine.Color = AUCmixedCols(hid,:);
            h(hid).patch.FaceColor = AUCmixedCols(hid,:);
            h(hid).edge(1).Color = AUCmixedCols(hid,:);
            h(hid).edge(2).Color = AUCmixedCols(hid,:);
        else
            h(hid) = plot(aucprecision, squeeze(nanmean(AUCNImCurves(:,hid,:),3)));
            h(hid).Color = AUCmixedCols(hid,:);
        end

    end
    plot([0, 1], [0, 1], 'k--')

    xlabel('False positives')
    ylabel('True positives')
    makepretty
end
% AllRFCor
if any(~isnan(AllRFCor(:)))
    subplot(4,3,11)
    clear h
    hold on
    for hid = 1:3
        if numel(size(AllFPCor))>2

            h(hid) = shadedErrorBar(aucprecision, squeeze(nanmean(AUCRFCurves(:,hid,:),3)),squeeze(nanstd(AUCRFCurves(:,hid,:),[],3))./sqrt(size(AUCRFCurves,3)-1));
            h(hid).mainLine.Color = AUCmixedCols(hid,:);
            h(hid).patch.FaceColor = AUCmixedCols(hid,:);
            h(hid).edge(1).Color = AUCmixedCols(hid,:);
            h(hid).edge(2).Color = AUCmixedCols(hid,:);
        else
            h(hid) = plot(aucprecision, squeeze(nanmean(AUCRFCurves(:,hid,:),3)));
            h(hid).Color = AUCmixedCols(hid,:);
        end
    end
    plot([0, 1], [0, 1], 'k--')

    xlabel('False positives')
    ylabel('True positives')
    makepretty
end
%% Now we make histograms for the Area Under the Curve scores
Incl1 = find(~cellfun(@isempty,UMTrackingPerformancePerMouse));
ChronicMice = find(cellfun(@(X) nanmean(X(1,:)),UMTrackingPerformancePerMouse(Incl1))>0);
ChronicMice = Incl1(ChronicMice);
pAUC = nan(4,3);

edges = [0:0.1:1];
subplot(4, 3, 3)
hold on
h1 = histogram(FingerprintAUC(1, ChronicMice), edges);
h1.FaceColor = nanmean(AUCCols([2,3],:),1);
% h1.FaceAlpha = 0.5;
h1.EdgeColor = nanmean(AUCCols([2,3],:),1);

h2 = histogram(FingerprintAUC(2, ChronicMice), edges);
h2.FaceColor = nanmean(AUCCols([1,2],:),1);
% h2.FaceAlpha = 0.5;
h2.EdgeColor = nanmean(AUCCols([1,2],:),1);


h3 = histogram(FingerprintAUC(3, ChronicMice), edges);
h3.FaceColor = nanmean(AUCCols([1,3],:),1);
% h2.FaceAlpha = 0.5;
h3.EdgeColor = nanmean(AUCCols([1,3],:),1);

line([nanmedian(FingerprintAUC(1, ChronicMice)), nanmedian(FingerprintAUC(1, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,3],:),1))
line([nanmedian(FingerprintAUC(2, ChronicMice)), nanmedian(FingerprintAUC(2, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,2],:),1))
line([nanmedian(FingerprintAUC(3, ChronicMice)), nanmedian(FingerprintAUC(3, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,3],:),1))

xlabel('AUC cross-correlation')
ylabel('Nr sessions')
makepretty

% Statistics
[~,pAUC(1,1)] = ttest(FingerprintAUC(1, ChronicMice),FingerprintAUC(2, ChronicMice));
[~,pAUC(1,2)] = ttest(FingerprintAUC(1, ChronicMice),FingerprintAUC(3, ChronicMice));
[~,pAUC(1,3)] = ttest(FingerprintAUC(2, ChronicMice),FingerprintAUC(3, ChronicMice));

disp(['Fingerprint correlation - Match/Non-match vs Within/Match: p=' num2str(round(pAUC(1,1).*3*100)./100)])
disp(['Fingerprint correlation - Match/Non-match vs Within/Non-match: p=' num2str(round(pAUC(1,2).*3*100)./100)])
disp(['Fingerprint correlation - Within/Match vs Within/Non-match: p=' num2str(round(pAUC(1,3).*3*100)./100)])

%%
subplot(4, 3, 6)
hold on
h1 = histogram(ACGAUC(1, :), edges);
h1.FaceColor = nanmean(AUCCols([2,3],:),1);
% h1.FaceAlpha = 0.5;
h1.EdgeColor = nanmean(AUCCols([2,3],:),1);

h2 = histogram(ACGAUC(2, :), edges);
h2.FaceColor = nanmean(AUCCols([1,2],:),1);
% h2.FaceAlpha = 0.5;
h2.EdgeColor = nanmean(AUCCols([1,2],:),1);

h3 = histogram(ACGAUC(3, :), edges);
h3.FaceColor = nanmean(AUCCols([1,3],:),1);
% h2.FaceAlpha = 0.5;
h3.EdgeColor = nanmean(AUCCols([1,3],:),1);

line([nanmedian(ACGAUC(1, :)), nanmedian(ACGAUC(1, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,3],:),1))
line([nanmedian(ACGAUC(2, :)), nanmedian(ACGAUC(2, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,1],:),1))
line([nanmedian(ACGAUC(3, :)), nanmedian(ACGAUC(3, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,3],:),1))

xlabel('AUC ACG correlation')
ylabel('Nr sessions')
makepretty

% Statistics
[~,pAUC(2,1)] = ttest(ACGAUC(1, ChronicMice),ACGAUC(2, ChronicMice));
[~,pAUC(2,2)] = ttest(ACGAUC(1, ChronicMice),ACGAUC(3, ChronicMice));
[~,pAUC(2,3)] = ttest(ACGAUC(2, ChronicMice),ACGAUC(3, ChronicMice));

disp(['Autocorrellogram correlation - Match/Non-match vs Within/Match: p=' num2str(round(pAUC(2,1).*3*100)./100)])
disp(['Autocorrellogram correlation - Match/Non-match vs Within/Non-match: p=' num2str(round(pAUC(2,2).*3*100)./100)])
disp(['Autocorrellogram correlation - Within/Match vs Within/Non-match: p=' num2str(round(pAUC(2,3).*3*100)./100)])

%% Natural images
if any(~isnan(NatImgAUC(:)))
    subplot(4, 3, 9)
    hold on
    h1 = histogram(NatImgAUC(1, ChronicMice), edges);
    h1.FaceColor = nanmean(AUCCols([2,3],:),1);
    % h1.FaceAlpha = 0.5;
    h1.EdgeColor = nanmean(AUCCols([2,3],:),1);

    h2 = histogram(NatImgAUC(2, ChronicMice), edges);
    h2.FaceColor = nanmean(AUCCols([1,2],:),1);
    % h2.FaceAlpha = 0.5;
    h2.EdgeColor = nanmean(AUCCols([1,2],:),1);


    h3 = histogram(NatImgAUC(3, ChronicMice), edges);
    h3.FaceColor = nanmean(AUCCols([1,3],:),1);
    % h2.FaceAlpha = 0.5;
    h3.EdgeColor = nanmean(AUCCols([1,3],:),1);

    line([nanmedian(NatImgAUC(1, ChronicMice)), nanmedian(NatImgAUC(1, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,3],:),1))
    line([nanmedian(NatImgAUC(2, ChronicMice)), nanmedian(NatImgAUC(2, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,2],:),1))
    line([nanmedian(NatImgAUC(3, ChronicMice)), nanmedian(NatImgAUC(3, ChronicMice))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,3],:),1))

    xlabel('AUC cross-correlation')
    ylabel('Nr sessions')
    makepretty

    % Statistics
    [~,pAUC(3,1)] = ttest(NatImgAUC(1, ChronicMice),NatImgAUC(2, ChronicMice));
    [~,pAUC(3,2)] = ttest(NatImgAUC(1, ChronicMice),NatImgAUC(3, ChronicMice));
    [~,pAUC(3,3)] = ttest(NatImgAUC(2, ChronicMice),NatImgAUC(3, ChronicMice));

    disp(['Fingerprint Natural Images - Match/Non-match vs Within/Match: p=' num2str(round(pAUC(3,1).*3*100)./100)])
    disp(['Fingerprint Natural Images - Match/Non-match vs Within/Non-match: p=' num2str(round(pAUC(3,2).*3*100)./100)])
    disp(['Fingerprint Natural Images - Within/Match vs Within/Non-match: p=' num2str(round(pAUC(3,3).*3*100)./100)])
end
%%
if any(~isnan(RFAUC(:)))
    subplot(4, 3, 12)
    hold on
    h1 = histogram(RFAUC(1, :), edges);
    h1.FaceColor = nanmean(AUCCols([2,3],:),1);
    % h1.FaceAlpha = 0.5;
    h1.EdgeColor = nanmean(AUCCols([2,3],:),1);


    h2 = histogram(RFAUC(2, :), edges);
    h2.FaceColor = nanmean(AUCCols([2,1],:),1);
    % h2.FaceAlpha = 0.5;
    h2.EdgeColor = nanmean(AUCCols([2,1],:),1);

    h3 = histogram(RFAUC(3, :), edges);
    h3.FaceColor = nanmean(AUCCols([1,3],:),1);
    % h2.FaceAlpha = 0.5;
    h3.EdgeColor = nanmean(AUCCols([1,3],:),1);
    line([nanmedian(RFAUC(1, :)), nanmedian(RFAUC(1, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,3],:),1))
    line([nanmedian(RFAUC(2, :)), nanmedian(RFAUC(2, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([2,1],:),1))
    line([nanmedian(RFAUC(3, :)), nanmedian(RFAUC(3, :))], get(gca, 'ylim'), 'Color', nanmean(AUCCols([1,3],:),1))

    xlabel('AUC RF distance')
    ylabel('Nr sessions')
    makepretty

    % Statistics
    [~,pAUC(4,1)] = ttest(RFAUC(1, ChronicMice),RFAUC(2, ChronicMice));
    [~,pAUC(4,2)] = ttest(RFAUC(1, ChronicMice),RFAUC(3, ChronicMice));
    [~,pAUC(4,3)] = ttest(RFAUC(2, ChronicMice),RFAUC(3, ChronicMice));

    disp(['Receptive field distance - Match/Non-match vs Within/Match: p=' num2str(round(pAUC(4,1).*3*100)./100)])
    disp(['Receptive field distance - Match/Non-match vs Within/Non-match: p=' num2str(round(pAUC(4,2).*3*100)./100)])
    disp(['Receptive field distance - Within/Match vs Within/Non-match: p=' num2str(round(pAUC(4,3).*3*100)./100)])

end
saveas(gcf, fullfile(SaveDir, 'FunctionScoreFigAcrossMice.fig'))
saveas(gcf, fullfile(SaveDir, 'FunctionScoreFigAcrossMice.bmp'))