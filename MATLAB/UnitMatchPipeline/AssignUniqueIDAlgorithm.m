function [MatchTable, UniqueIDConversion] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam)
%% Input:
% UniqueIDConversion (struct) with at east per cluster: OriginalClusID, GoodID & recsesAll
% A MatchTable (output table from UnitMatch)
% UMparam (struct) with at least: ProbabilityThreshold and GoodUnitsOnly

%% Extract cluster information
AllClusterIDs = UniqueIDConversion.OriginalClusID;
UniqueIDLiberal = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
OriUniqueID = UniqueIDLiberal; % Liberal
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
% Re-initialize UID
[PairID3,PairID4]=meshgrid(OriUniqueID(Good_Idx));
MatchTable.UID1 = PairID3(:);
MatchTable.UID2 = PairID4(:);

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
    %
    %     figure; plot(plotvec,hw,'g-')
    %     hold on; plot(plotvec,ha,'r-')
    %     plot(plotvec,hn,'b-')

else
    UMparam.UsedProbability = UMparam.ProbabilityThreshold;
end

%% Initial pairing based on matchscore
Pairs = [MatchTable.UID1(MatchTable.MatchProb>UMparam.UsedProbability) MatchTable.UID2(MatchTable.MatchProb>UMparam.UsedProbability)]; % These are UIDs, not Indexes!
Pairs(diff(Pairs,[],2)==0,:)=[]; %Remove own matches
Pairs = sort(Pairs,2,'ascend'); % Only use each unique pair once
Pairs = unique(Pairs,'stable','rows');

%% Cut down on the pairs that do not exceed ProbabilityThreshold on both cross validations
% Average of two cross-validations and sort by that
[~,tblidx] = ismember(Pairs,[MatchTable.UID1 MatchTable.UID2],'rows'); %  For the first of two cross-validations, find which rows of the table contain the pairs
MatchProbabilityOri = MatchTable.MatchProb(tblidx);

[~,tblidx] = ismember(Pairs,[MatchTable.UID2 MatchTable.UID1],'rows'); % Do the same for the other way around
MatchProbabilityFlip = MatchTable.MatchProb(tblidx);

Pairs(MatchProbabilityOri<UMparam.UsedProbability|MatchProbabilityFlip<UMparam.UsedProbability,:) = []; % don't bother with these

%% Sort by average match probability
MatchProbability = nanmean(cat(2,MatchProbabilityOri,MatchProbabilityFlip),2); %Average
MatchProbability(MatchProbabilityOri<UMparam.UsedProbability|MatchProbabilityFlip<UMparam.UsedProbability) = []; % don't bother with these
[~,sortidx] = sort(MatchProbability,'descend');
Pairs = Pairs(sortidx,:); %Pairs as UID, but now sorted by match probability

%% assigning of Unique ID
disp('Assigning correct Unique ID values now')
if ~isempty(Pairs)
    nMatchesConservative = 0;
    nMatchesLiberal = 0;
    nMatches = 0;
    for id = 1:size(Pairs,1)

        %% Liberal
        TheseOriUids = OriUniqueID(ismember(UniqueIDLiberal,Pairs(id,:))); % Find all units that are currently having the same UniqueIDLiberal as either one of the current pairs, and their original unique ID as they are known in matchtable
        Idx = find(ismember(OriUniqueID,TheseOriUids)); % Find all units that have to be considered
        UniqueIDLiberal(Idx) = min(UniqueIDLiberal(Idx)); % Assign them all with the minimum UniqueIDLiberal
        nMatchesLiberal = nMatchesLiberal + 1;

        %% Conservative
        % Identify the rows in the table containing TheseOriUids with
        % either of the currently considered units in Pairs
        tblidx = find((ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,Pairs(id,2)) | ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,Pairs(id,1)) | ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,Pairs(id,1)) | ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,Pairs(id,2))) & ~(MatchTable.UID1==MatchTable.UID2)); % !

        if all(MatchTable.MatchProb(tblidx)>UMparam.UsedProbability) % All of these should survive the threshold
            % All UIDs with this UID should now become the new UID
            UniqueIDConservative(Idx) = min(UniqueIDConservative(Idx)); % Assign them all with the minimum UniqueIDLiberal
            nMatchesConservative = nMatchesConservative + 1;
        end

        %% Intermediate
        % Far away days can be ignored - this will be decided later
        TheseRecs = GoodRecSesID(Pairs(id,:));
        OtherRecs = GoodRecSesID(ismember(UniqueIDLiberal,UniqueIDLiberal(Pairs(id,:))));
        TheseOriUids(~ismember(OtherRecs,[TheseRecs-1; TheseRecs; TheseRecs+1;]))=[]; % Remove far days
        tblidx = find((ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,Pairs(id,2)) | ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,Pairs(id,1)) | ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,Pairs(id,1)) | ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,Pairs(id,2))) & ~(MatchTable.UID1==MatchTable.UID2)); % !
        if all(MatchTable.MatchProb(tblidx)>UMparam.UsedProbability) % All of these should survive the threshold
            % All UIDs with this UID should now become the new UID
            UniqueID(Idx) = min(UniqueID(Idx)); % Assign them all with the minimum UniqueIDLiberal
            nMatches = nMatches+ 1;
        end

    end
    disp([num2str(nMatchesLiberal) ' liberal matches/' num2str(length(UniqueIDLiberal))])
    disp([num2str(nMatchesConservative) ' conservative matches/' num2str(length(UniqueIDLiberal))])
    disp([num2str(nMatches) ' intermediate matches/' num2str(length(UniqueID))])


end

%% Replace in table
[PairID3,PairID4] = meshgrid(UniqueIDConservative(Good_Idx));
MatchTable.UID1Conservative = PairID3(:);
MatchTable.UID2Conservative = PairID4(:);
[PairID3,PairID4] = meshgrid(UniqueIDLiberal(Good_Idx));
MatchTable.UID1Liberal = PairID3(:);
MatchTable.UID2Liberal = PairID4(:);
[PairID3,PairID4] = meshgrid(UniqueID(Good_Idx));
MatchTable.UID1 = PairID3(:);
MatchTable.UID2 = PairID4(:);

%% Replace in UniqueIDConversion
UniqueIDConversion.UniqueIDConservative = UniqueIDConservative;
UniqueIDConversion.UniqueIDLiberal = UniqueIDLiberal;
UniqueIDConversion.UniqueID = UniqueID;



return