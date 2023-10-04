function summaryFunctionalPlots(UMFiles, whichMetric, groupVector, UseKSLabels)
    %% Will plot summary plots: distribution, ROC and AUC. 
    % UMFiles: list cells contains path to UnitMatch.m files
    % whichMetric: will compute distributions/ROC/AUC on either 'Corr', 'Rank', or 'Sig'. 
    % groupVector: same size as UMFiles, to group e.g. recordings from the
    % same mouse together.
    % UseKSLabels: use KS labels for matching

    %% Define parameters
    
    % Initialize
    if ~exist('whichMetric','var') || isempty(whichMetric)
        whichMetric = 'Corr';
    end

    if ~exist('groupVector','var')
        groupVector = 1:length(UMFiles);
    end
    groups = unique(groupVector);
    groupColor = gray(length(groups)+1);

    if ~exist('UseKSLabels','var')
        UseKSLabels = 0;
    end

    switch whichMetric
        case 'Corr'
            fprintf("Taking the correlation values!\n")
            FPNames = {'FRDiff','ACGCorr','refPopCorr','natImCorr','natImRespCorr','natImScaledRespCorr'};
            stepsz = [0.1 0.1 0.1 0.1 0.1 0.1];
            minVal = [0 -1 -1 -1 -1 -1];
            maxVal = [15 1 1 1 1 1];
            flipROC = [0 1 1 1 1 1];
        case 'Rank'
            fprintf("Taking the rank!\n")
            FPNames = {'FRRank','ACGRank','refPopRank','natImRank','natImRespRank','natImScaledRespRank'};
            stepsz = [1 1 1 1 1 1];
            minVal = [1 1 1 1 1 1];
            maxVal = [21 21 21 21 21 21];
            flipROC = [0 0 0 0 0 0];
        case 'Sig'
            fprintf("Taking the rank!\n")
            FPNames = {'FRSig','ACGSig','refPopSig','natImSig','natImRespSig','natImScaledRespSig'};
            stepsz = [1 1 1 1 1 1];
            minVal = [1 1 1 1 1 1];
            maxVal = [21 21 21 21 21 21];
            flipROC = [0 0 0 0 0 0];
    end

    histBins = cell(1,numel(FPNames));
    histBinsCenter = cell(1,numel(FPNames));
    for fpIdx = 1:numel(FPNames)
        histBins{fpIdx} = minVal(fpIdx):stepsz(fpIdx):maxVal(fpIdx);
        histBinsCenter{fpIdx} = histBins{fpIdx}(1:end-1) + diff(histBins{fpIdx})/2;
    end
    ROCBins = 0:0.01:1;
    minMatches = 20;
    durLim = 10*60;

    
    %% Loop over mice to get all Distributions / ROCs / AUCs
    
    FPSum = struct();
    deltaDays = cell(1, length(UMFiles));
    numMatchedUnits = cell(1, length(UMFiles));
    InitialDrift = cell(1,length(UMFiles));
    FixedDrift =  cell(1,length(UMFiles));
    maxAvailableUnits = cell(1, length(UMFiles));
    for midx = 1:length(UMFiles)
        %% Load data
    
        fprintf('Reference %s...\n', UMFiles{midx})
    
        tmpfile = dir(UMFiles{midx});
        if isempty(tmpfile)
            continue
        end
    
        fprintf('Loading the data...\n')
        tic
        load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UniqueIDConversion', 'UMparam');
        toc
    
        sessIDs = unique(MatchTable.RecSes1);
    
        % Initialize
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            FPSum.(FPNameCurr).Distr{midx} = nan(numel(histBins{fpIdx})-1, 3, numel(sessIDs)-1, numel(sessIDs));
            FPSum.(FPNameCurr).AUC{midx} = nan(3, numel(sessIDs)-1, numel(sessIDs));
            FPSum.(FPNameCurr).ROC{midx} = nan(numel(ROCBins), 3, numel(sessIDs)-1, numel(sessIDs));
        end
    
        %%% HACK -- Can remove later
        if ~iscell(UMparam.AllRawPaths)
            for ii = 1:numel(UMparam.AllRawPaths)
                tmp{ii} = UMparam.AllRawPaths(ii);
            end
            UMparam.AllRawPaths = tmp;
        end
    
        %% Loop through pairs of sessions
    
        fprintf('Looping through days...\n')
        tic
        days = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.AllRawPaths, 'uni', 0), 'uni', 0);
        days = cell2mat(days) - days{1};
        deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        InitialDrift{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        FixedDrift{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        for sess1Idx = 1:numel(sessIDs)-1
    
            sess1 = sessIDs(sess1Idx);
            day1 = regexp(UMparam.AllRawPaths{sess1}.folder,'\d*-\d*-\d*','match'); day1 = datenum(day1{1});
            meta = ReadMeta2(UMparam.AllRawPaths{sess1}.folder);
            durSess1 = str2double(meta.fileTimeSecs);
            if durSess1 < durLim 
                continue
            end

            for sess2Idx = sess1Idx+1:numel(sessIDs)
    
                sess2 = sessIDs(sess2Idx);
                day2 = regexp(UMparam.AllRawPaths{sess2}.folder,'\d*-\d*-\d*','match'); day2 = datenum(day2{1});
                deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);
                meta = ReadMeta2(UMparam.AllRawPaths{sess2}.folder);
                durSess2 = str2double(meta.fileTimeSecs);
                if durSess2 < durLim
                    continue
                end
    
                %% Cut table to specific days
    
                MatchTable_2sess = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);
    
                %% Number of matches
    
                if ~UseKSLabels
                    %%% CHECK THAT THIS MAKES SENSE
                    %%% CHOOSE BASED ON UID
                    matchedUnitsIdx = (MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2); % using Unique ID
                    %%% OR RECOMPUTE
%                     [~,~,idx,~] = getPairsAcross2Sess(MatchTable_2sess, UMparam.ProbabilityThreshold);
%                     matchedUnitsIdx = zeros(size(MatchTable_2sess,1),1);
%                     matchedUnitsIdx(idx) = 1;
                else
                    matchedUnitsIdx = (MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2);
                end
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum(matchedUnitsIdx)/2; % Divided by two because looking both ways -- can be non-integer
                
                maxAvailableUnits{midx}(sess1Idx,sess2Idx) = min([length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess1))) length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess2)))]);%
                
                %% Extract drift if present
    
                if isfield(UMparam,'drift')
                    InitialDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,1),2); % Drift in recording 1 is 1 vs 2, etc.
                    FixedDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,2),2); % Drift in recording 1 is 1 vs 2, etc.
                end
    
                %% Looping through fingerprints
                for fpIdx = 1:numel(FPNames)
                    FPNameCurr = FPNames{fpIdx};
    
                    if ~any(ismember(MatchTable_2sess.Properties.VariableNames, FPNameCurr))
                        fprintf('%s has no %s in table...\n', UMFiles{midx}, FPNameCurr)
                        continue
                    end
    
                    %% select pairs within and across days
    
                    % Special criteria
                    if contains(FPNameCurr,'NatI') || contains(FPNameCurr,'NImg') %%% clean
                        % Select pairs of which the first neuron has good
                        % test-retest reliability
                        Unit1ID = double(MatchTable_2sess.ID1) + 10e9*MatchTable_2sess.RecSes1; % identifier for each pair
                        uUnit1ID = unique(Unit1ID); % list of units
                        reliability = MatchTable_2sess((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2),:).NatImCorr; % test-retest reliability of each unit
    
                        validPairs = ismember(Unit1ID, uUnit1ID(reliability > 0.2));
                    else
                        validPairs = ones(size(MatchTable_2sess,1),1);
                    end
    
                    % Extract groups: "within", "match", "non-match"
                    if ~UseKSLabels
                        % WithinIdx = find((MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                        % MatchIdx = find((MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                        % NonMatchIdx = find((MatchTable_2sess.UID1 ~= MatchTable_2sess.UID2) & validPairs); % Not the same unit

                        %%% CHECK THAT THIS MAKES SENSE
                        WithinIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                        MatchIdx = find(matchedUnitsIdx & validPairs); % Across session, same unit (cross-validation)
                        NonMatchIdx = find(~matchedUnitsIdx & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); % Across session, not the same unit
                    else
                        WithinIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                        MatchIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                        NonMatchIdx = find((MatchTable_2sess.ID1 ~= MatchTable_2sess.ID2) & validPairs); % Not the same unit
                    end
    
                    %% Check that passed the criteria
    
                    % Check recordings duration
                   

                    % Condition to go ahead
                    goAhead = numel(MatchIdx)/2 >= minMatches && ... % minimum number of matches
                        ~all(isnan(MatchTable_2sess.(FPNameCurr)(MatchIdx))); % not all nans
    
                    if goAhead
                        %% Compute fingerprint correlations and AUCs
    
                        FPCorr = MatchTable_2sess.(FPNameCurr);
    
                        % Compute distributions
                        hw = histcounts(FPCorr(WithinIdx), histBins{fpIdx}) ./ length(WithinIdx);
                        hm = histcounts(FPCorr(MatchIdx), histBins{fpIdx}) ./ length(MatchIdx);
                        hn = histcounts(FPCorr(NonMatchIdx), histBins{fpIdx}) ./ length(NonMatchIdx);
                        % Save
                        FPSum.(FPNameCurr).Distr{midx}(:,:,sess1Idx,sess2Idx) = [hw', hm', hn'];
    
                        % Compute ROCs and AUCs
                        [ROC1, AUC1] = getAUC(FPCorr,MatchIdx,NonMatchIdx,ROCBins,flipROC(fpIdx));
                        [ROC2, AUC2] = getAUC(FPCorr,WithinIdx,MatchIdx,ROCBins,flipROC(fpIdx));
                        [ROC3, AUC3] = getAUC(FPCorr,WithinIdx,NonMatchIdx,ROCBins,flipROC(fpIdx));
                        % Save
                        FPSum.(FPNameCurr).ROC{midx}(:,:,sess1Idx,sess2Idx) = [ROC1,ROC2,ROC3];
                        FPSum.(FPNameCurr).AUC{midx}(:,sess1Idx,sess2Idx) = [AUC1, AUC2, AUC3];
                    end
                end
            end
        end
        toc
    end
    
    %% Figure -- example mouse
    
    % Best example mouse
    % midx = find(cellfun(@(X) nansum(X(:)),numMatchedUnits) == max(cellfun(@(X)  nansum(X(:)),numMatchedUnits)),1,'first');
    
    % Number of matched units (matrix)
    for midx = 1:numel(UMFiles)
        figure('Position', [80 700 1700 170],'Name', fileparts(fileparts(UMFiles{midx})));
        s = subplot(1,numel(FPNames)+2,1);
        imagesc(numMatchedUnits{midx})
        xticks(1:size(deltaDays{midx},2))
        xticklabels(days)
        yticks(1:size(deltaDays{midx},1))
        yticklabels(days(1:end-1))
        axis equal tight
        colorbar
        c = colormap("gray"); colormap(s(1), flipud(c));
        subplot(1,numel(FPNames)+2,2)
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,groupColor(groupVector(midx),:),'filled')
        ylabel('Number of matches')
        xlabel('\Deltadays')
    
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            s = subplot(1,numel(FPNames)+2,fpIdx+2);
            imagesc(squeeze(FPSum.(FPNameCurr).AUC{midx}(1,:,:)))
            xticks(1:size(deltaDays{midx},2))
            xticklabels(days)
            yticks(1:size(deltaDays{midx},1))
            yticklabels(days(1:end-1))
            axis equal tight
            colormap(s, "RedBlue")
            clim([0 1])
            title(sprintf('Fingerprint %s', FPNameCurr))
        end
    end
    
    %% Figure -- all mice
    % Build matrices across mice -- average in groups according to groupVector
    distMatrix = cell(1,numel(FPNames));
    ROCMatrix = cell(1,numel(FPNames));
    AUCMatrix = cell(1,numel(FPNames));
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        distMatrix{fpIdx} = nan(numel(histBinsCenter{fpIdx}), 3, length(groups));
        ROCMatrix{fpIdx} = nan(numel(ROCBins), 3, length(groups));
        AUCMatrix{fpIdx} = nan(3, length(groups));
        for gg = 1:length(groups)
            % Distributions
            tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).Distr(groupVector == gg), 'uni', 0);
            distMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);
    
            % ROC
            tmp = cellfun(@(x) x(:,:,:), FPSum.(FPNameCurr).ROC(groupVector == gg), 'uni', 0);
            ROCMatrix{fpIdx}(:,:,gg) = nanmean(cat(3,tmp{:}),3);
    
            % AUC
            tmp = cellfun(@(x) x(:,:), FPSum.(FPNameCurr).AUC(groupVector == gg), 'uni', 0);
            AUCMatrix{fpIdx}(:,gg) = nanmean(cat(2,tmp{:}),2);
        end
    end
    
    % Plot
    distrCols = [0 0.7 0; 1 0 0; 0 0 0.7]; % Within / Match / Non-match
    ROCCols = [1 0 0; 0.5 0.5 0.5; 0 0 0]; % Match vs Non-match / Within vs Match / Within vs Non-match
    figure('Position', [400 270 800 700]);
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
    
        % Plot distribution
        subplot(4,numel(FPNames),0*numel(FPNames)+fpIdx); hold all
        for hid = 1:3
            h = shadedErrorBar(histBinsCenter{fpIdx}, nanmean(distMatrix{fpIdx}(:,hid,:),3), ...
                nanstd(distMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(distMatrix{fpIdx}(:,hid,:)),3)));
            h.mainLine.Color = distrCols(hid,:);
            if ~isempty(h.patch)
                h.patch.FaceColor = distrCols(hid,:);
                h.edge(1).Color = 'none';
                h.edge(2).Color = 'none';
            end
        end
        title(sprintf('%s', FPNameCurr))
        xlabel(whichMetric)
        if strcmp(whichMetric,'Rank'); xticks([1 10 maxVal(fpIdx)]); xticklabels({'1','10',sprintf('>%d',maxVal(fpIdx)-1)}); end
        ylabel('Proportion')
        offsetAxes
        makepretty
    
        % Plot ROC
        subplot(4,numel(FPNames),1*numel(FPNames)+fpIdx); hold all
        for hid = 3:-1:1
            h = shadedErrorBar(ROCBins, nanmean(ROCMatrix{fpIdx}(:,hid,:),3), ...
                nanstd(ROCMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(ROCMatrix{fpIdx}(:,hid,:)),3)));
            h.mainLine.Color = ROCCols(hid,:);
            if ~isempty(h.patch)
                h.patch.FaceColor = ROCCols(hid,:);
                h.edge(1).Color = 'none';
                h.edge(2).Color = 'none';
            end
        end
        plot([0 1], [0 1], 'k--')
        xlim([0 1])
        ylim([0 1])
        axis equal tight
        xlabel('False positives')
        ylabel('Hits')
        offsetAxes
        makepretty
    
        % Plot AUC
        subplot(4,numel(FPNames),2*numel(FPNames)+fpIdx); hold all
        for hid = 1:3
            scatter(hid + 0.1*randn(1,numel(groups)), AUCMatrix{fpIdx}(hid,:), 30, ROCCols(hid,:), 'filled')
        end
        ylim([0 1])
        hline(0.5, 'k--')
        xticks(1:3)
        xticklabels({'M v NM', 'W v M', 'W v NM'})
        xtickangle(45)
        ylabel('AUC')
        offsetAxes
        makepretty
    
        % Plot stability of AUC with delta days
        subplot(4,numel(FPNames),3*numel(FPNames)+fpIdx); hold all
        for midx = 1:length(UMFiles)
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:)),20,groupColor(groupVector(midx),:),'filled')
            end
        end
        xlabel('\Delta days')
        ylabel('AUC')
        ylim([0 1])
        hline(0.5)
        offsetAxes
        makepretty
    end
    
    %% Additional figures
    
    % Plot number of matches as a function of delta days
    figure;
    hold all
    for midx = 1:length(UMFiles)
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,groupColor(groupVector(midx),:),'filled')
    end
    ylabel('Number of matches')
    xlabel('\Deltadays')
    
    % AUC as a function of delta days
    figure;
    for aucIdx = 1:3
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),20,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('\Delta days')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
    
    % Dependence of AUC on number of matched units
    figure;
    for aucIdx = 1:3
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(numMatchedUnits{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),20,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('Number of matched neurons')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
    
    
    %% Dependence of number matched units on drift
    
%     figure('name','NrUnits versus drift')
%     for midx = 1:length(UMFiles)
%         subplot(1,2,1)
%         scatter(InitialDrift{midx},numMatchedUnits{midx}./MaxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
%         hold on
%         xlabel('Drift (Eucl Distance)')
%         ylabel('Number of matches')
%         xlim([0 UMparam.NeighbourDist])
%         title('Initial Drift')
%     
%         subplot(1,2,2)
%         scatter(FixedDrift{midx},numMatchedUnits{midx}./ MaxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
%         hold on
%         xlabel('Drift (Eucl Distance)')
%         ylabel('Proportion of matches')
%         xlim([0 UMparam.NeighbourDist])
%     
%     end
%     InitialDrift(cellfun(@isempty,InitialDrift)) = [];
%     InitialDrift = cellfun(@(X) X(:),InitialDrift,'uni',0);
%     InitialDrift(InitialDrift>UMparam.maxdist) = nan;
%     MaxAvailableUnits(cellfun(@isempty,MaxAvailableUnits)) = [];
%     MaxAvailableUnits = cellfun(@(X) X(:),MaxAvailableUnits,'uni',0);
%     numMatchedUnits(cellfun(@isempty,numMatchedUnits)) = [];
%     numMatchedUnits = cellfun(@(X) X(:),numMatchedUnits,'uni',0);
%     InitialDrift = cat(1,InitialDrift{:});
%     MaxAvailableUnits = cat(1,MaxAvailableUnits{:});
%     numMatchedUnits = cat(1,numMatchedUnits{:});
%     FixedDrift(cellfun(@isempty,FixedDrift)) = [];
%     FixedDrift = cellfun(@(X) X(:),FixedDrift,'uni',0);
%     FixedDrift = cat(1,FixedDrift{:});
%     
%     
%     FixedDrift(FixedDrift>UMparam.maxdist) = nan;
%     InitialDrift(InitialDrift>UMparam.maxdist) = nan;
%     
%     
%     PercNeurons = numMatchedUnits./MaxAvailableUnits;
%     [r,p] = corr(InitialDrift(~isnan(InitialDrift)),PercNeurons(~isnan(InitialDrift)));
%     subplot(1,2,1)
%     title(['Initial Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])
%     
%     subplot(1,2,2)
%     [r,p] = corr(FixedDrift(~isnan(FixedDrift)),PercNeurons(~isnan(FixedDrift)));
%     title(['Corrected Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])
end