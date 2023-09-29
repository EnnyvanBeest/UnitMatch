%% Define parameters

% Initialize
TakeRank = 0; %if 0 , take cross-correlation scores (these may be less informative than rank)
if TakeRank
    fprintf("Taking the rank!\n")
    stepsz = 1;
    minVal = -20;
    maxVal = 1;
    FPNames = {'RankScore','ACGRankScore','FRRankScore','NImgRankScore'};
else
    fprintf("Taking the correlation values!\n")
    stepsz = 0.1;
    minVal = -1;
    maxVal = 1;
    FPNames = {'FingerprintCor','ACGCorr','FRDiff','NatImCorr'};
end
histBins = minVal:stepsz:maxVal;
UseKSLabels = PrepareClusInfoparams.RunPyKSChronicStitched;
AUCPrecision = 0:0.01:1;
minMatches = 10;

%% Loop over mice to get all Distributions / ROCs / AUCs

clear FPSum
deltaDays = cell(1, length(UMFolders));
numMatchedUnits = cell(1, length(UMFolders));
for midx = 1:length(UMFolders)
    %% Load data
    
    fprintf('Reference %s...\n', UMFolders{midx})

    tmpfile = dir(fullfile(SaveDir, UMFolders{midx}, 'UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end

    fprintf('Loading the data...\n')
    tic
    load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
    toc

    sessIDs = unique(MatchTable.RecSes1);

    % Initialize
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        FPSum.(FPNameCurr).Distr{midx} = nan(numel(histBins)-1, 3, numel(sessIDs)-1, numel(sessIDs));
        FPSum.(FPNameCurr).AUC{midx} = nan(3, numel(sessIDs)-1, numel(sessIDs));
        FPSum.(FPNameCurr).ROC{midx} = nan(numel(AUCPrecision), 3, numel(sessIDs)-1, numel(sessIDs));
    end
                
    %%% HACK
    if ~iscell(UMparam.AllRawPaths)
        for ii = 1:numel(UMparam.AllRawPaths)
            tmp{ii} = UMparam.AllRawPaths(ii);
        end
        UMparam.AllRawPaths = tmp;
    end

    %% Loop through pairs of sessions

    days = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.AllRawPaths, 'uni', 0), 'uni', 0);
    days = cell2mat(days) - days{1};
    deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    for sess1Idx = 1:numel(sessIDs)-1
    
        sess1 = sessIDs(sess1Idx);
        day1 = regexp(UMparam.AllRawPaths{sess1}.folder,'\d*-\d*-\d*','match'); day1 = datenum(day1{1});
        
        for sess2Idx = sess1Idx+1:numel(sessIDs)
        
            sess2 = sessIDs(sess2Idx);
            day2 = regexp(UMparam.AllRawPaths{sess2}.folder,'\d*-\d*-\d*','match'); day2 = datenum(day2{1});
            deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);

            %% Cut table to specific days

            MatchTable_pair = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);

            %% Number of matches

            if ~UseKSLabels
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum((MatchTable_pair.UID1 == MatchTable_pair.UID2) & ...
                    (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2))/2;
            else
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum((MatchTable_pair.ID1 == MatchTable_pair.ID2) & ...
                    (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2))/2;
            end

            %% Looping through fingerprints
            for fpIdx = 1:numel(FPNames)
                FPNameCurr = FPNames{fpIdx};

                if ~any(ismember(MatchTable_pair.Properties.VariableNames, FPNameCurr))
                    disp([UMFolders{midx} ' has no ' FPNameCurr ' in table...'])
                    continue
                end

                %% select pairs within and across days

                % Special criteria
                if contains(FPNameCurr,'NatI') || contains(FPNameCurr,'NImg') %%% clean
                    % Select pairs of which the first neuron has good
                    % test-retest reliability
                    Unit1ID = double(MatchTable_pair.ID1) + 10e9*MatchTable_pair.RecSes1; % identifier for each pair
                    uUnit1ID = unique(Unit1ID); % list of units
                    reliability = MatchTable_pair((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2),:).NatImCorr; % test-retest reliability of each unit
                    
                    validPairs = ismember(Unit1ID, uUnit1ID(reliability > 0.2));
                else 
                    validPairs = ones(size(MatchTable_pair,1),1);
                end 

                % Extract groups: "within", "match", "non-match"
                if ~UseKSLabels
                    WithinIdx = find((MatchTable_pair.UID1 == MatchTable_pair.UID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                    MatchIdx = find((MatchTable_pair.UID1 == MatchTable_pair.UID2) & (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                    NonMatchIdx = find((MatchTable_pair.UID1 ~= MatchTable_pair.UID2) & validPairs); % Not the same unit
                else
                    WithinIdx = find((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 == MatchTable_pair.RecSes2) & validPairs); %Within session, same unit (cross-validation)
                    MatchIdx = find((MatchTable_pair.ID1 == MatchTable_pair.ID2) & (MatchTable_pair.RecSes1 ~= MatchTable_pair.RecSes2) & validPairs); %Across session, same unit (cross-validation)
                    NonMatchIdx = find((MatchTable_pair.ID1 ~= MatchTable_pair.ID2) & validPairs); % Not the same unit
                end

                %% Check that passed the criteria

                % Condition to go ahead
                goAhead = numel(MatchIdx)/2 >= minMatches && ... % minimum number of matches
                    ~all(isnan(MatchTable_pair.(FPNameCurr)(MatchIdx))); % not all nans

                if ~goAhead
                    fprintf('Not enough matches (or all nans) for sessions %d & %d (%d). Skipping.\n', sess1, sess2, numel(MatchIdx))
                else
                    %% Compute fingerprint correlations and AUCs

                    FPCorr = MatchTable_pair.(FPNameCurr);
                
                    % Compute distributions
                    hw = histcounts(FPCorr(WithinIdx), histBins) ./ length(WithinIdx);
                    hm = histcounts(FPCorr(MatchIdx), histBins) ./ length(MatchIdx);
                    hn = histcounts(FPCorr(NonMatchIdx), histBins) ./ length(NonMatchIdx);
                    % Save
                    FPSum.(FPNameCurr).Distr{midx}(:,:,sess1Idx,sess2Idx) = [hw', hm', hn'];

                    % Compute ROCs and AUCs
                    [ROC1, AUC1] = getAUC(FPCorr,MatchIdx,NonMatchIdx,AUCPrecision);
                    [ROC2, AUC2] = getAUC(FPCorr,WithinIdx,MatchIdx,AUCPrecision);
                    [ROC3, AUC3] = getAUC(FPCorr,WithinIdx,NonMatchIdx,AUCPrecision);
                    % Save
                    FPSum.(FPNameCurr).AUC{midx}(:,sess1Idx,sess2Idx) = [AUC1, AUC2, AUC3];
                    FPSum.(FPNameCurr).ROC{midx}(:,:,sess1Idx,sess2Idx) = [ROC1,ROC2,ROC3];
                end
            end
        end
    end
end

%% Plots -- example mouse

mouseColor = winter(length(UMFolders));
% Best example mouse
midx = find(cellfun(@(X) nansum(X(:)),numMatchedUnits) == max(cellfun(@(X)  nansum(X(:)),numMatchedUnits)),1,'first')
% Number of matched units (matrix)
figure;
subplot(121)
imagesc(numMatchedUnits{midx})
xticks(1:size(deltaDays{midx},2))
xticklabels(days)
yticks(1:size(deltaDays{midx},1))
yticklabels(days(1:end-1))
axis equal tight
colorbar
c = colormap("gray"); colormap(flipud(c));
subplot(122)
scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,mouseColor(midx,:),'filled')
ylabel('Number of matches')
xlabel('\Deltadays')

figure
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};
    subplot(1,numel(FPNames),fpIdx)
    imagesc(squeeze(FPSum.(FPNameCurr).AUC{midx}(1,:,:)))
    xticks(1:size(deltaDays{midx},2))
    xticklabels(days)
    yticks(1:size(deltaDays{midx},1))
    yticklabels(days(1:end-1))
    axis equal tight
    colormap("RedBlue")
    clim([0 1])
    title(sprintf('Fingerprint %s', FPNameCurr))
end

%%



figure
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};

end
Vector = minVal+stepsz/2:stepsz:maxVal-stepsz/2;

%% Plots -- all mice

midx = length(UMFolders);
figure;
hold all
for midx = 1:length(UMFolders)
    scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}),20,mouseColor(midx,:),'filled')
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
        for midx = 1:length(UMFolders)
            if ~isempty(FPSum.(FPNameCurr).AUC{midx})
                scatter(deltaDays{midx}(:), mat2vec(FPSum.(FPNameCurr).AUC{midx}(aucIdx,:,:)),20,mouseColor(midx,:),'filled')
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