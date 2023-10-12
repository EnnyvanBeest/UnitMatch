function [FPSum, days, deltaDays, numMatchedUnits, maxAvailableUnits] = summaryFunctionalPlots(UMFiles, whichMetric, groupVector, UseKSLabels)
    %% Will plot summary plots: distribution, ROC and AUC. 
    % UMFiles: list cells contains path to UnitMatch.m files
    % whichMetric: will compute distributions/ROC/AUC on either 'Corr', 'Rank', or 'Sig'. 
    % groupVector: same size as UMFiles, to group e.g. recordings from the
    % same mouse together.
    % UseKSLabels: use KS labels for matching

    %% Define parameters
    
    % Initialize
    if ~exist('whichMetric','var') || isempty(whichMetric)
        whichMetric = 'Rank';
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
            FPNames = {'FRDiff','ACGCorr','natImRespCorr','refPopCorr'};
            stepsz = [0.1 0.1 0.1 0.1];
            minVal = [0 -1 -1 -1];
            maxVal = [15 1 1 1];
            flipROC = [0 1 1 1];
        case 'Rank'
            fprintf("Taking the rank!\n")
            FPNames = {'FRRank','ACGRank','natImRespRank','refPopRank'};
            stepsz = [1 1 1 1];
            minVal = [0.5 0.5 0.5 0.5];
            maxVal = [20.5 20.5 20.5 20.5];
            flipROC = [0 0 0 0];
        case 'Sig'
            fprintf("Taking the sig!\n")
            FPNames = {'FRSig','ACGSig','natImRespSig','refPopSig'};
            stepsz = [0.1 0.1 0.1 0.1];
            minVal = [0 0 0 0];
            maxVal = [1 1 1 1];
            flipROC = [0 0 0 0];
    end

    histBins = cell(1,numel(FPNames));
    histBinsCenter = cell(1,numel(FPNames));
    for fpIdx = 1:numel(FPNames)
        histBins{fpIdx} = minVal(fpIdx):stepsz(fpIdx):maxVal(fpIdx);
        histBinsCenter{fpIdx} = histBins{fpIdx}(1:end-1) + diff(histBins{fpIdx})/2;
        if strcmp(whichMetric,'Rank')
            histBins{fpIdx}(end+1) = inf;
            histBinsCenter{fpIdx}(end+1) = histBinsCenter{fpIdx}(end)+1;
        end
    end
    ROCBins = 0:0.01:1;
    minMatches = 20;
    durLim = 10*60;

    
    %% Loop over mice to get all Distributions / ROCs / AUCs
    
    FPSum = struct();
    days = cell(1, length(UMFiles));
    deltaDays = cell(1, length(UMFiles));
    numMatchedUnits = cell(1, length(UMFiles));
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
        load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam');
        toc
    
        sessIDs = unique(MatchTable.RecSes1);
    
        % Initialize
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            FPSum.(FPNameCurr).Distr{midx} = nan(numel(histBins{fpIdx})-1, 3, numel(sessIDs)-1, numel(sessIDs));
            FPSum.(FPNameCurr).AUC{midx} = nan(2, numel(sessIDs)-1, numel(sessIDs));
            FPSum.(FPNameCurr).ROC{midx} = nan(numel(ROCBins), 2, numel(sessIDs)-1, numel(sessIDs));
        end
    
        %% Loop through pairs of sessions
    
        fprintf('Looping through days...\n')
        tic
        days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
        days{midx} = cell2mat(days{midx}) - days{midx}{1};
        deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
        for sess1Idx = 1:numel(sessIDs)-1
    
            sess1 = sessIDs(sess1Idx);
            day1 = days{midx}(sess1Idx);
            meta = ReadMeta2(UMparam.RawDataPaths{sess1}.folder);
            durSess1 = str2double(meta.fileTimeSecs);
            if durSess1 < durLim 
                continue
            end

            for sess2Idx = sess1Idx+1:numel(sessIDs)
    
                sess2 = sessIDs(sess2Idx);
                day2 = days{midx}(sess2Idx);
                deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);
                meta = ReadMeta2(UMparam.RawDataPaths{sess2}.folder);
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

                        withinMatchIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                        withinNonMatchIdx = find((MatchTable_2sess.ID1 ~= MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, different unit (cross-validation)
                        acrossMatchIdx = find(matchedUnitsIdx & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); % Across session, same unit (cross-validation)
                        acrossNonMatchIdx = find(~matchedUnitsIdx & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); % Across session, not the same unit
                    else
                        withinMatchIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                        withinNonMatchIdx = find((MatchTable_2sess.ID1 ~= MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 == MatchTable_2sess.RecSes2) & validPairs); %Within session, different unit (cross-validation)
                        acrossMatchIdx = find((MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                        acrossNonMatchIdx = find((MatchTable_2sess.ID1 ~= MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2) & validPairs); % Not the same unit
                    end
    
                    %% Check that passed the criteria
    
                    % Check recordings duration
                   

                    % Condition to go ahead
                    goAhead = numel(acrossMatchIdx)/2 >= minMatches && ... % minimum number of matches
                        ~all(isnan(MatchTable_2sess.(FPNameCurr)(acrossMatchIdx))); % not all nans
    
                    if goAhead
                        %% Compute fingerprint correlations and AUCs
    
                        FPCorr = MatchTable_2sess.(FPNameCurr);
    
                        % Compute distributions
                        hw = histcounts(FPCorr(withinMatchIdx), histBins{fpIdx}) ./ length(withinMatchIdx);
                        hm = histcounts(FPCorr(acrossMatchIdx), histBins{fpIdx}) ./ length(acrossMatchIdx);
                        hn = histcounts(FPCorr(withinNonMatchIdx), histBins{fpIdx}) ./ length(acrossNonMatchIdx);
                        % Save
                        FPSum.(FPNameCurr).Distr{midx}(:,:,sess1Idx,sess2Idx) = [hw', hm', hn'];
    
                        % Compute ROCs and AUCs
                        [ROC1, AUC1] = getAUC(FPCorr,acrossMatchIdx,acrossNonMatchIdx,ROCBins,flipROC(fpIdx));
                        [ROC2, AUC2] = getAUC(FPCorr,withinMatchIdx,withinNonMatchIdx,ROCBins,flipROC(fpIdx));
                        % Save
                        FPSum.(FPNameCurr).ROC{midx}(:,:,sess1Idx,sess2Idx) = [ROC1,ROC2];
                        FPSum.(FPNameCurr).AUC{midx}(:,sess1Idx,sess2Idx) = [AUC1, AUC2];
                    end
                end
            end
        end
        toc
    end
    
    %% Figure -- example mouse

    % Number of matched units (matrix)
    for midx = 1:numel(UMFiles)
        figure('Position', [80 700 1700 170],'Name', fileparts(fileparts(UMFiles{midx})));
        s = subplot(1,numel(FPNames)+2,1);
        imagesc(numMatchedUnits{midx})
        xticks(1:size(deltaDays{midx},2))
        xticklabels(days{midx})
        yticks(1:size(deltaDays{midx},1))
        yticklabels(days{midx}(1:end-1))
        axis equal tight
        colorbar
        c = colormap("gray"); colormap(s(1), flipud(c));
        subplot(1,numel(FPNames)+2,2); hold all
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),10,groupColor(groupVector(midx),:),'filled')
        hline(minMatches)
        ylabel('Number of matches')
        xlabel('Deltadays')
    
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            s = subplot(1,numel(FPNames)+2,fpIdx+2);
            imagesc(squeeze(FPSum.(FPNameCurr).AUC{midx}(1,:,:)))
            xticks(1:size(deltaDays{midx},2))
            xticklabels(days{midx})
            yticks(1:size(deltaDays{midx},1))
            yticklabels(days{midx}(1:end-1))
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
        ROCMatrix{fpIdx} = nan(numel(ROCBins), 2, length(groups));
        AUCMatrix{fpIdx} = nan(2, length(groups));
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
    ROCCols = [1 0 0; 0 0.7 0]; % across Match vs. non-match / within Match vs. non-match
    figure('Position', [400 100 1000 800]);
    clear bsave
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        
        % Plot distribution
        subplot(4,numel(FPNames),0*numel(FPNames)+fpIdx); hold all
        for hid = [3 1 2]
%             distr2plt = cumsum(distMatrix{fpIdx}(:,hid,:));
            distr2plt = distMatrix{fpIdx}(:,hid,:);
            plot(histBinsCenter{fpIdx}, nanmean(distr2plt,3), 'color', distrCols(hid,:))
%             h = shadedErrorBar(histBinsCenter{fpIdx}, nanmean(distr2plt,3), ...
%                 nanstd(distr2plt,[],3)./sqrt(sum(~isnan(distr2plt),3)));
%             h.mainLine.Color = distrCols(hid,:);
%             if ~isempty(h.patch)
%                 h.patch.FaceColor = distrCols(hid,:);
%                 h.edge(1).Color = 'none';
%                 h.edge(2).Color = 'none';
%             end
        end
        title(sprintf('%s', FPNameCurr))
        xlabel(whichMetric)
        if strcmp(whichMetric,'Rank'); xticks([1 10 histBinsCenter{fpIdx}(end)]); xticklabels({'1','10',sprintf('>%d',histBinsCenter{fpIdx}(end-1))}); end
        ylabel('Proportion')
        offsetAxes
        makepretty
    
        % Plot ROC
        subplot(4,numel(FPNames),1*numel(FPNames)+fpIdx); hold all
        for hid = 2:-1:1
            plot(ROCBins, nanmean(ROCMatrix{fpIdx}(:,hid,:),3), 'color', ROCCols(hid,:))
%             h = shadedErrorBar(ROCBins, nanmean(ROCMatrix{fpIdx}(:,hid,:),3), ...
%                 nanstd(ROCMatrix{fpIdx}(:,hid,:),[],3)./sqrt(sum(~isnan(ROCMatrix{fpIdx}(:,hid,:)),3)));
%             h.mainLine.Color = ROCCols(hid,:);
%             if ~isempty(h.patch)
%                 h.patch.FaceColor = ROCCols(hid,:);
%                 h.edge(1).Color = 'none';
%                 h.edge(2).Color = 'none';
%             end
            AUCtext = num2str(sprintf('%0.0f\x00B1%0.0f', nanmean(AUCMatrix{fpIdx}(hid,:)*100,2), ...
                nanstd(AUCMatrix{fpIdx}(hid,:)*100,[],2)./sqrt(sum(~isnan(AUCMatrix{fpIdx}(hid,:)),2))));
            text(0.6,hid*0.2,AUCtext,'Color',ROCCols(hid,:))
        end
        plot([0 1], [0 1], 'k--')
        xlim([0 1])
        ylim([0 1])
        xticks([0 1])
        yticks([0 1])
        axis equal tight
        xlabel('False positives')
        ylabel('Hits')
        offsetAxes
        makepretty
    
        % Plot stability of AUC with delta days
        subplot(4,numel(FPNames),2*numel(FPNames)+fpIdx); hold all
        for midx = 1:length(UMFiles)
            xDays = deltaDays{midx}(:);
            xDays(xDays == 0) = 10^(-0.1);
            yVal = mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:));
            nanIdx = isnan(yVal);
            xDays(nanIdx) = [];
            yVal(nanIdx) = [];
            if ~isempty(FPSum.(FPNameCurr).AUC{midx}) && numel(unique(xDays)) > 1
                scatter(log10(xDays),yVal,10,ROCCols(1,:),'filled')
                X = [ones(numel(xDays),1), xDays];
                b = (X\yVal);
                plot(log10(1:max(xDays)), b(1) + b(2)*(1:max(xDays)), 'color',ROCCols(1,:),'LineWidth',1);
                scatter(-0.1,nanmean(mat2vec(FPSum.(FPNameCurr).AUC{midx}(2,:,:))),20,ROCCols(2,:),'filled')
            end
        end
        xlabel('Delta days')
        ylabel('AUC')
        xticks([-0.1 log10([1 10 100])])
        xticklabels({'within','1','10','100'})
        ylim([0 1])
        hline(0.5)
        offsetAxes
        makepretty

        % Plot stability of AUC with delta days
        bsave{fpIdx} = nan(length(UMFiles),2);
        subplot(4,numel(FPNames),3*numel(FPNames)+fpIdx); hold all
        for midx = 1:length(UMFiles)
            xDays = deltaDays{midx}(:);
            xDays(xDays == 0) = 10^(-0.1);
            yVal = mat2vec(FPSum.(FPNameCurr).AUC{midx}(1,:,:));
            nanIdx = isnan(yVal);
            xDays(nanIdx) = [];
            yVal(nanIdx) = [];
            if ~isempty(FPSum.(FPNameCurr).AUC{midx}) && numel(unique(xDays)) > 1
                X = [ones(numel(xDays),1), xDays];
                b = (X\yVal);
                plot(log10(1:max(xDays)), b(1) + b(2)*(1:max(xDays)), 'color',ROCCols(1,:),'LineWidth',1);
                bsave{fpIdx}(midx,:) = b;
                scatter(-0.1,nanmean(mat2vec(FPSum.(FPNameCurr).AUC{midx}(2,:,:))),20,ROCCols(2,:),'filled')
            end
        end
        xlabel('Delta days')
        ylabel('AUC')
        xticks([-0.1 log10([1 10 100])])
        xticklabels({'within','1','10','100'})
        ylim([0 1])
        hline(0.5)
        offsetAxes
        makepretty
    end
    set(gcf,'Renderer','painters')
    
    %% Additional figures
    
    %
    slope = nan(numel(FPNames),length(groups));
    for fpIdx = 1:numel(FPNames)
        for gg = 1:length(groups)
            slope(fpIdx,gg) = nanmean(bsave{fpIdx}(groupVector == gg,2));
        end
    end

    %
    for fpIdx = 1:numel(FPNames)
        notNanIdx = ~isnan(slope(fpIdx,:));
        n = sum(notNanIdx);
        % n = sum(~isnan(AUCMatrix{fpIdx}(1,:)),2);
        mu = nanmean(AUCMatrix{fpIdx}(1,notNanIdx)*100,2);
        se = nanstd(AUCMatrix{fpIdx}(1,notNanIdx)*100,[],2)./sqrt(n);
        fprintf('%s across: %0.0f\x00B1%0.0f (n = %d mice)\n', FPNames{fpIdx}, mu, se, n)

        % n = sum(~isnan(AUCMatrix{fpIdx}(2,:)),2);
        mu = nanmean(AUCMatrix{fpIdx}(2,notNanIdx)*100,2);
        se = nanstd(AUCMatrix{fpIdx}(2,notNanIdx)*100,[],2)./sqrt(n);
        fprintf('%s within: %0.0f\x00B1%0.0f (n = %d mice)\n', FPNames{fpIdx}, mu, se, n)

        % n = sum(~isnan(slope(fpIdx,:)),2);
        mu = nanmedian(slope(fpIdx,notNanIdx),2);
        se = mad(slope(fpIdx,notNanIdx));
        fprintf('%s slope: %0.5f\x00B1%0.5f (n = %d mice)\n', FPNames{fpIdx}, mu, se, n)
    end


    % Plot number of matches as a function of delta days
    figure;
    hold all
    for midx = 1:length(UMFiles)
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),10,groupColor(groupVector(midx),:),'filled')
    end
    ylabel('Number of matches')
    xlabel('Deltadays')
    
    % AUC as a function of delta days
    figure;
    for aucIdx = 1:2
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),10,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('Delta days')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
    
    % Dependence of AUC on number of matched units
    figure;
    for aucIdx = 1:2
        for fpIdx = 1:numel(FPNames)
            FPNameCurr = FPNames{fpIdx};
            subplot(2,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
            hold all
            for midx = 1:length(UMFiles)
                if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                    scatter(numMatchedUnits{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),10,groupColor(groupVector(midx),:),'filled')
                end
            end
            xlabel('Number of matched neurons')
            ylabel('AUC')
            title(sprintf('Fingerprint %s', FPNameCurr))
            ylim([0 1])
            hline(0.5)
        end
    end
end