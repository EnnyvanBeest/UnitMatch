function [MatchTable, UniqueIDConversion] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam)
%% Input:
% UniqueIDConversion (struct) with at east per cluster: OriginalClusID, GoodID & recsesAll
% A MatchTable (output table from UnitMatch)
% UMparam (struct) with at least: ProbabilityThreshold and GoodUnitsOnly
disp('Assigning correct Unique ID values now')


%% Extract cluster information
AllClusterIDs = UniqueIDConversion.OriginalClusID;
UniqueIDLiberal = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
OriUniqueID = UniqueIDLiberal; % original
UniqueIDConservative = UniqueIDLiberal;% Conservative equivalent
UniqueID = UniqueIDLiberal; %Intermediate;
if UMparam.GoodUnitsOnly
    Good_Idx = find(UniqueIDConversion.GoodID); %Only care about good units at this point
else
    Good_Idx = 1:length(UniqueIDConversion.GoodID);
    disp('Use all units including MUA and noise')
end
GoodRecSesID = UniqueIDConversion.recsesAll;
recOpt = unique(GoodRecSesID);
nRec = length(recOpt);
nClus = numel(Good_Idx);

%% Initial replacement
UniqueIDConversion.UniqueIDConservative = UniqueIDConservative;
UniqueIDConversion.UniqueIDLiberal = UniqueIDLiberal;
UniqueIDConversion.UniqueID = UniqueID;

%% Make variables ready for indexing only Good_Idx
UniqueIDLiberal = UniqueIDLiberal(Good_Idx);
OriUniqueID = OriUniqueID(Good_Idx);
UniqueIDConservative = UniqueIDConservative(Good_Idx);
UniqueID = UniqueID(Good_Idx);
GoodRecSesID = GoodRecSesID(Good_Idx);
% Extract match prob
MatchProb = reshape(MatchTable.MatchProb,nClus,nClus);

%% Find probability threshold to use
if UMparam.UseDatadrivenProbThrs
    stepsz = 0.1;
    binedges = [0:stepsz:1];
    plotvec = stepsz/2:stepsz:1-stepsz/2;
    hw = histcounts(MatchTable.MatchProb(MatchTable.ID1==MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2),binedges)./sum(MatchTable.ID1==MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);
    %     hn = histcounts(MatchTable.MatchProb(MatchTable.ID1~=MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2),binedges)./sum(MatchTable.ID1~=MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);
    %     ha = histcounts(MatchTable.MatchProb(MatchTable.RecSes1 ~= MatchTable.RecSes2),binedges)./sum(MatchTable.RecSes1 ~= MatchTable.RecSes2);

    %     find sudden increase and take that as probability threshold
    UMparam.UsedProbability = plotvec(diff(hw)>0.1);
    disp(['Used Probability = ' num2str(UMparam.UsedProbability)])
    %
    %     figure; plot(plotvec,hw,'g-')
    %     hold on; plot(plotvec,ha,'r-')
    %     plot(plotvec,hn,'b-')

else
    UMparam.UsedProbability = UMparam.ProbabilityThreshold;
end

%% Initial pairing based on matchscore
[r,c] = find(MatchProb>UMparam.UsedProbability);
Pairs = [r,c];
Pairs(diff(Pairs,[],2)==0,:)=[]; %Remove own matches
Pairs = sort(Pairs,2,'ascend'); % Only use each unique pair once
Pairs = unique(Pairs,'stable','rows');

%% Cut down on the pairs that do not exceed ProbabilityThreshold on both cross validations
% Average of two cross-validations and sort by that
MatchProbabilityOri = nan(0,1);
MatchProbabilityFlip = nan(0,1);
batchsz = 1000;
numBatch = ceil(size(Pairs,1)./batchsz);
for batchid = 1:numBatch % In batches to prevent memory issues
    idx = (batchid-1)*batchsz+1:batchid*batchsz;
    idx(idx>size(Pairs,1)) = [];

    MatchProbabilityOri = cat(1,MatchProbabilityOri,diag(MatchProb(Pairs(idx,1),Pairs(idx,2))));
    MatchProbabilityFlip = cat(1,MatchProbabilityFlip,diag(MatchProb(Pairs(idx,2),Pairs(idx,1))));
end
NonSurvivalIdx = MatchProbabilityOri<UMparam.UsedProbability|MatchProbabilityFlip<UMparam.UsedProbability;
Pairs(NonSurvivalIdx,:) = []; % don't bother with these
if ~isempty(Pairs)

    MatchProbabilityOri(NonSurvivalIdx) = [];
    MatchProbabilityFlip(NonSurvivalIdx) = [];
    %% Sort by average match probability
    MatchProbability = nanmean(cat(2,MatchProbabilityOri,MatchProbabilityFlip),2); %Average
    [~,sortidx] = sort(MatchProbability,'descend');
    Pairs = Pairs(sortidx,:); %Pairs as UID, but now sorted by match probability

    %% assigning of Unique ID
    fprintf(1,'Assigning pairs. Progress: %3d%%',0)
    if ~isempty(Pairs)
        nMatchesConservative = 0;
        nMatchesLiberal = 0;
        nMatches = 0;
        for id = 1:size(Pairs,1)
            fprintf(1,'\b\b\b\b%3.0f%%',id/size(Pairs,1)*100)

            %% Conservative - most stringent
            % Identify the rows in the table containing TheseOriUids with
            % either of the currently considered units in Pairs
            TheseOriUids = OriUniqueID(ismember(UniqueIDConservative,UniqueIDConservative(Pairs(id,:)))); % Find all units that are currently having the same UniqueIDConservative as either one of the current pairs, and their original unique ID as they are known in matchtable
            Idx = find(ismember(OriUniqueID,TheseOriUids)); % Find all units that have to be considered

            % For conservative, all pairs need to match with eachother
            AllRelevantProbs = nanmin(cat(3,MatchProb(Idx,Idx),MatchProb(Idx,Idx)'),[],3); % Relevant probs for Conservative
            AllRelevantProbs = AllRelevantProbs(triu(true(size(AllRelevantProbs)),1));
            if all(AllRelevantProbs>UMparam.UsedProbability) % All of these should survive the threshold
                % All UIDs with this UID should now become the new UID
                UniqueIDConservative(Idx) = min(UniqueIDConservative(Idx)); % Assign them all with the minimum UniqueIDLiberal
                nMatchesConservative = nMatchesConservative + 1;
            end

            %% Intermediate
            TheseOriUids = OriUniqueID(ismember(UniqueID,UniqueID(Pairs(id,:)))); % Find all units that are currently having the same UniqueIDConservative as either one of the current pairs, and their original unique ID as they are known in matchtable
            Idx = find(ismember(OriUniqueID,TheseOriUids)); % Find all units that have to be considered

            Ses2Consid = [GoodRecSesID(Pairs(id,:))-1 GoodRecSesID(Pairs(id,:)) GoodRecSesID(Pairs(id,:))+1];
            Ses2Consid = unique(Ses2Consid(:));
            Idx(~ismember(GoodRecSesID(Idx),Ses2Consid)) = []; % Don't consider these

            AllRelevantProbs = nanmin(cat(3,MatchProb(Idx,Idx),MatchProb(Idx,Idx)'),[],3); % Relevant probs for Conservative
            AllRelevantProbs = AllRelevantProbs(triu(true(size(AllRelevantProbs)),1));
            if all(AllRelevantProbs>UMparam.UsedProbability) % All of these should survive the threshold
                % All UIDs with this UID should now become the new UID
                UniqueID(Idx) = min(UniqueID(Idx)); % Assign them all with the minimum UniqueIDLiberal
                nMatches = nMatches+ 1;
            end

            %% Liberal
            TheseOriUids = OriUniqueID(ismember(UniqueIDLiberal,UniqueIDLiberal(Pairs(id,:)))); % Find all units that are currently having the same UniqueIDLiberal as either one of the current pairs, and their original unique ID as they are known in matchtable
            Idx = find(ismember(OriUniqueID,TheseOriUids)); % Find all units that have to be considered
            UniqueIDLiberal(Idx) = min(UniqueIDLiberal(Idx)); % Assign them all with the minimum UniqueIDLiberal
            nMatchesLiberal = nMatchesLiberal + 1;
        end
        disp('')
        disp([num2str(nMatchesLiberal) ' liberal matches/' num2str(length(UniqueIDLiberal))])
        disp([num2str(nMatchesConservative) ' conservative matches/' num2str(length(UniqueIDLiberal))])
        disp([num2str(nMatches) ' intermediate matches/' num2str(length(UniqueID))])
    end
end
%% Replace in table
[PairID3,PairID4] = meshgrid(UniqueIDConservative);
MatchTable.UID1Conservative = PairID3(:);
MatchTable.UID2Conservative = PairID4(:);
[PairID3,PairID4] = meshgrid(UniqueIDLiberal);
MatchTable.UID1Liberal = PairID3(:);
MatchTable.UID2Liberal = PairID4(:);
[PairID3,PairID4] = meshgrid(UniqueID);
MatchTable.UID1 = PairID3(:);
MatchTable.UID2 = PairID4(:);

%% Replace in UniqueIDConversion
UniqueIDConversion.UniqueIDConservative(Good_Idx) = UniqueIDConservative;
UniqueIDConversion.UniqueIDLiberal(Good_Idx) = UniqueIDLiberal;
UniqueIDConversion.UniqueID(Good_Idx) = UniqueID;

return