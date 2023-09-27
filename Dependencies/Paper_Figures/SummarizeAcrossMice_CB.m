%% Define parameters

% Initialize
TakeRank = 0; %if 0 , take cross-correlation scores (these may be less informative than rank)
if TakeRank
    fprintf("Taking the rank!\n")
    stepsz = 1;
    minVal = -20;
    maxVal = 1;
else
    fprintf("Taking the correlation values!\n")
    stepsz = 0.1;
    minVal = -20;
    maxVal = 1;
end
bins = minVal:stepsz:maxVal;
Vector = minVal+stepsz/2:stepsz:maxVal-stepsz/2;
UseKSLabels = PrepareClusInfoparams.RunPyKSChronicStitched;
AUCPrecision = 0:0.01:1;
maxMatches = 30;

FPNames = {'FingerprintCor','ACGCorr','FRDiff','NatImCorr'};
% FPNames = {'RankScore','ACGRankScore','FRRankScore','NImgRankScore'};

%% Loop over mice to get all Distributions / ROCs / AUCs

clear FPSum
deltaDays = cell(1, length(MiceOpt));
numMatchedUnits = cell(1, length(MiceOpt));
for midx = 1:length(MiceOpt)
    %% Load data
    
    fprintf('Reference %s...\n', MiceOpt{midx})

    tmpfile = dir(fullfile(SaveDir, MiceOpt{midx}, 'UnitMatch', 'UnitMatch.mat'));
    if isempty(tmpfile)
        continue
    end

    fprintf('Loading the data...\n')
    tic
    load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
    toc

    sessIDs = unique(MatchTable.RecSes1);

    %% Loop through pairs of sessions

    deltaDays{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    numMatchedUnits{midx} = nan(numel(sessIDs)-1,numel(sessIDs));
    for sess1Idx = 1:numel(sessIDs)-1
    
        sess1 = sessIDs(sess1Idx);
        day1 = regexp(UMparam.AllRawPaths(sess1).folder,'\d*-\d*-\d*','match'); day1 = datenum(day1{1});
        
        for sess2Idx = sess1Idx+1:numel(sessIDs)
        
            sess2 = sessIDs(sess2Idx);
            day2 = regexp(UMparam.AllRawPaths(sess2).folder,'\d*-\d*-\d*','match'); day2 = datenum(day2{1});
            deltaDays{midx}(sess1,sess2) = abs(day2 - day1);

            %% Cut table to specific days

            MatchTable_pair = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);

            %% Looping through fingerprints
            for fpIdx = 1:numel(FPNames)
                FPNameCurr = FPNames{fpIdx};

                if ~any(ismember(MatchTable_pair.Properties.VariableNames, FPNameCurr))
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

                numMatchedUnits{midx}(sess1,sess2) = numel(MatchIdx);

                %% Check that passed the criteria

                % Condition to go ahead
                goAhead = numel(MatchIdx) >= maxMatches && ... % minimum number of matches
                    ~all(isnan(MatchTable_pair.(FPNameCurr)(MatchIdx))); % not all nans

                if ~goAhead
                    fprintf('Not enough matches (or all nans) for sessions %d & %d (%d). Skipping.\n', sess1, sess2, numel(MatchIdx))
                    FPSum.(FPNameCurr).Distr{midx}(:,:,sess1,sess2) = nan(numel(bins)-1, 3);
                    FPSum.(FPNameCurr).AUC{midx}(:,sess1,sess2) = nan(1,3);
                    FPSum.(FPNameCurr).ROC{midx}(:,:,sess1,sess2) = nan(numel(AUCPrecision), 3);
                else
                    %% Compute fingerprint correlations and AUCs

                    FPCorr = MatchTable_pair.(FPNameCurr);
                
                    % Compute distributions
                    hw = histcounts(FPCorr(WithinIdx), bins) ./ length(WithinIdx);
                    hm = histcounts(FPCorr(MatchIdx), bins) ./ length(MatchIdx);
                    hn = histcounts(FPCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
                    % Save
                    FPSum.(FPNameCurr).Distr{midx}(:,:,sess1,sess2) = [hw', hm', hn'];

                    % Compute ROCs and AUCs
                    [ROC1, AUC1] = getAUC(FPCorr,MatchIdx,NonMatchIdx,AUCPrecision);
                    [ROC2, AUC2] = getAUC(FPCorr,WithinIdx,MatchIdx,AUCPrecision);
                    [ROC3, AUC3] = getAUC(FPCorr,WithinIdx,NonMatchIdx,AUCPrecision);
                    % Save
                    FPSum.(FPNameCurr).AUC{midx}(:,sess1,sess2) = [AUC1, AUC2, AUC3];
                    FPSum.(FPNameCurr).ROC{midx}(:,:,sess1,sess2) = [ROC1,ROC2,ROC3];
                end
            end
        end
    end
end

%% Plots

figure;
imagesc(numMatchedUnits{midx})
axis equal tight

figure;
mouseColor = parula(length(MiceOpt));
for aucIdx = 1:3
    for fpIdx = 1:numel(FPNames)
        FPNameCurr = FPNames{fpIdx};
        subplot(3,numel(FPNames),(aucIdx-1)*numel(FPNames)+fpIdx);
        hold all
        for midx = 1:length(MiceOpt)
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