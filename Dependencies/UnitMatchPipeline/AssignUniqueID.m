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

%% Cut down on the pairs
% Average of two cross-validations and sort by that
[~,tblidx] = ismember(Pairs,[MatchTable.UID1 MatchTable.UID2],'rows');
MatchProbabilityOri = MatchTable.MatchProb(tblidx);

[~,tblidx] = ismember(Pairs,[MatchTable.UID2 MatchTable.UID1],'rows');
MatchProbabilityFlip = MatchTable.MatchProb(tblidx);
Pairs(MatchProbabilityOri<0.5|MatchProbabilityFlip<0.5,:) = []; % don't bother with these
MatchProbability = nanmean(cat(2,MatchProbabilityOri,MatchProbabilityFlip),2);

MatchProbability(MatchProbabilityOri<0.5|MatchProbabilityFlip<0.5) = []; % don't bother with these

[~,sortidx] = sort(MatchProbability,'descend');
Pairs = Pairs(sortidx,:); %Pairs, but now sorted by match probability

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
RecCombinations = nan(2,0);
for recid = 1:length(RecOpt)
    for recid2 = recid:length(RecOpt)
        RecCombinations = cat(2,RecCombinations,[recid recid2]');
    end
end
disp('Assigning correct Unique ID values now')
RMPair = [];
parfor recid = 1:size(RecCombinations,2)
    Idx = find(GoodRecSesID(Pairs(:,1)) == RecCombinations(1,recid) & GoodRecSesID(Pairs(:,2)) == RecCombinations(2,recid));
    SubPairs = Pairs(Idx,:); %Select pairs in this day combinations
    Utmp = UniqueID;
    if ~isempty(SubPairs)
        nMatches = 0;
        for id = 1:size(SubPairs,1)
            %  check: It should also match with all the other pairs that were
            % already assigned!
            % All units currently identified as this UniqueID
            TheseOriUids = OriUniqueID(ismember(Utmp,Utmp(SubPairs(id,:))));
            TheseOriUids(~ismember(GoodRecSesID(TheseOriUids),RecCombinations(:,recid))) = [];
            % All of these need to match with the new one, if added
            tblidx = find(((ismember(MatchTable.UID1,TheseOriUids)&ismember(MatchTable.UID2,SubPairs(id,2))) | (ismember(MatchTable.UID2,TheseOriUids)&ismember(MatchTable.UID1,SubPairs(id,2)))) & ~(MatchTable.UID1==MatchTable.UID2)); % !
            if ~all(MatchTable.MatchProb(tblidx)>UMparam.ProbabilityThreshold)
                RMPair = [RMPair Idx(id)];
            else
                Utmp(SubPairs(id,2)) = Utmp(SubPairs(id,1));
                nMatches = nMatches + 1;
            end
        end
        disp(['Recording ' num2str(RecCombinations(1,recid)) ' vs ' num2str( RecCombinations(2,recid)) ': ' num2str(nMatches)])
    end   
end
%% Final assignment
disp('Final assignment')
Pairs(RMPair,:) = []; % Surviving Pairs
Pairs = sortrows(Pairs);
for uid = 1:size(Pairs,1)
    UniqueID(Pairs(uid,:)) = min(UniqueID(Pairs(uid,:))); % Serial assignment;
end

%% Replace in table
[PairID3,PairID4]=meshgrid(UniqueID(Good_Idx));
MatchTable.UID1 = PairID3(:);
MatchTable.UID2 = PairID4(:);

return