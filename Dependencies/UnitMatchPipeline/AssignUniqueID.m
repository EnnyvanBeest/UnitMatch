function [UniqueID, MatchTable] = AssignUniqueID(MatchTable,clusinfo,Path4UnitNPY,param)
AllClusterIDs = clusinfo.cluster_id;
% nses = length(AllDecompPaths);
% OriginalClusID = AllClusterIDs; % Original cluster ID assigned by KS
UniqueID = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
OriUniqueID = UniqueID;
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')
end
GoodRecSesID = clusinfo.RecSesID;
RecOpt = unique(GoodRecSesID);


%% Initial pairing based on matchscore
Pairs = [MatchTable.UID1(MatchTable.MatchProb>param.ProbabilityThreshold) MatchTable.UID2(MatchTable.MatchProb>param.ProbabilityThreshold)]; %
Pairs(diff(Pairs,[],2)==0,:)=[]; %Remove own matches
Pairs = sort(Pairs,2,'ascend'); % Only use each unique pair once
Pairs = unique(Pairs,'stable','rows');

%% ISI violations (for over splits matching)
if param.removeoversplits
    CurrentSPPath = [];
    ISIViolationsScore = nan(1,size(Pairs,1));
    fprintf(1,'Computing functional properties for determining ISI. Progress: %3d%%',0)
    for pairid = 1:size(Pairs,1)
        if GoodRecSesID(Pairs(pairid,1)) == GoodRecSesID(Pairs(pairid,2))
            tmppath = strsplit(Path4UnitNPY{find(Good_Idx==Pairs(pairid,1))},'RawWaveforms');
            % Load sp for correct day
            if isempty(CurrentSPPath) || ~strcmp(CurrentSPPath{1},tmppath{1})
                if length(param.KSDir)>1
                    tmp = matfile(fullfile(param.KSDir{GoodRecSesID(Pairs(pairid,1))},'PreparedData.mat'));
                else %Stitched
                    tmp = matfile(fullfile(param.KSDir{1},'PreparedData.mat'));
                end
                sp = tmp.sp;
                CurrentSPPath = tmppath;
            end
            idx1 = sp.spikeTemplates == AllClusterIDs((Pairs(pairid,1)));
            idx2 = sp.spikeTemplates == AllClusterIDs((Pairs(pairid,2)));
            DifScore = diff(sort([sp.st(idx1); sp.st(idx2)]));
            ISIViolationsScore(pairid) = sum(DifScore.*1000<1.5)./length(DifScore);
            fprintf(1,'\b\b\b\b%3.0f%%',pairid/size(Pairs,1)*100)
        end
    end
    fprintf('\n')
    disp(['Removing ' num2str(sum(ISIViolationsScore>0.05)) ' matched oversplits, as merging them will violate ISI >5% of the time'])
    Pairs(ISIViolationsScore>0.05,:)=[];
end

%% Serial assigning of Unique ID (Day by day)
disp('Assigning correct Unique ID values now')
for recid = 1:length(RecOpt)
    for recid2 = 1:length(RecOpt)
        if recid2<recid
            continue
        end
        SubPairs = Pairs(GoodRecSesID(Pairs(:,1)) == recid & GoodRecSesID(Pairs(:,2)) == recid2,:); %Select days

        % Average of two cross-validations and sort by that
        MatchProbability = arrayfun(@(X) MatchTable.MatchProb(ismember(MatchTable.UID1,SubPairs(X,1))&ismember(MatchTable.UID2,SubPairs(X,2))),1:size(SubPairs,1));
        MatchProbabilityFlip = arrayfun(@(X) MatchTable.MatchProb(ismember(MatchTable.UID1,SubPairs(X,2))&ismember(MatchTable.UID2,SubPairs(X,1))),1:size(SubPairs,1));
        MatchProbability = nanmean(cat(1,MatchProbability,MatchProbabilityFlip),1); 
        SubPairs(MatchProbability<0.5,:) = []; % don't bother with these
        MatchProbability(MatchProbability<0.5) = []; % don't bother with these


        [MatchProbability,sortidx] = sort(MatchProbability,'descend');
        SubPairs = SubPairs(sortidx,:); %Pairs, but now sorted by match probability

        for id = 1:size(SubPairs,1)
            % Matchprobability should be high enough
            tblidx1 = find(ismember(MatchTable.UID1,SubPairs(id,1))&ismember(MatchTable.UID2,SubPairs(id,2)));
            if ~(MatchTable.MatchProb(tblidx1) > param.ProbabilityThreshold) %Requirement 1, match probability should be high enough
                continue
            end
            % Find the cross-validated version of this pair, this should also have
            % high enough probability
            tblidx2 = find(ismember(MatchTable.UID1,SubPairs(id,2))&ismember(MatchTable.UID2,SubPairs(id,1)));
            if ~(MatchTable.MatchProb(tblidx2) > param.ProbabilityThreshold) %Requirement 1, match probability should be high enough
                continue
            end
            % Extra check: It should also match with all the other pairs that were
            % already assigned!
            % All units currently identified as this UniqueID
            TheseOriUids = OriUniqueID(ismember(UniqueID,UniqueID(SubPairs(id,1))));
            TheseOriUids(GoodRecSesID(TheseOriUids)<recid | GoodRecSesID(TheseOriUids)>recid2) = [];
            % All of these need to match with the new one, if added
            tblidx = find(((ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,SubPairs(id,2))) | (ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,SubPairs(id,2)))) & ~(MatchTable.UID1==MatchTable.UID2)); % !
            if ~all(MatchTable.MatchProb(tblidx)>param.ProbabilityThreshold)
                continue
            end
     
            UniqueID(SubPairs(id,2)) = UniqueID(SubPairs(id,1)); %Survived, assign
        end
    end
end
%% Replace in table
[PairID3,PairID4]=meshgrid(UniqueID(Good_Idx));
MatchTable.UID1 = PairID3(:);
MatchTable.UID2 = PairID4(:);



return