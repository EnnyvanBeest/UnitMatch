function UniqueID = AssignUniqueID(MatchTable,clusinfo,sp,param)
AllClusterIDs = clusinfo.cluster_id;
% nses = length(AllDecompPaths);
% OriginalClusID = AllClusterIDs; % Original cluster ID assigned by KS
UniqueID = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
GoodRecSesID = clusinfo.RecSesID;

%% Initial pairing based on matchscore
Pairs = [MatchTable.ID1(MatchTable.MatchProb>param.ProbabilityThreshold)+1 MatchTable.ID2(MatchTable.MatchProb>param.ProbabilityThreshold)+1]; % 0-indexed
Pairs(diff(Pairs,[],2)==0,:)=[]; %Remove own matches
%% ISI violations (for over splits matching)
ISIViolationsScore = nan(1,size(Pairs,1));
fprintf(1,'Computing functional properties similarity. Progress: %3d%%',0)
for pairid= 1:size(Pairs,1)
    if GoodRecSesID(Pairs(pairid,1)) == GoodRecSesID(Pairs(pairid,2))
        idx1 = sp.spikeTemplates == AllClusterIDs((Pairs(pairid,1)))&sp.RecSes == GoodRecSesID(Pairs(pairid,1));
        idx2 = sp.spikeTemplates == AllClusterIDs((Pairs(pairid,2)))&sp.RecSes == GoodRecSesID(Pairs(pairid,2));
        DifScore = diff(sort([sp.st(idx1); sp.st(idx2)]));
        ISIViolationsScore(pairid) = sum(DifScore.*1000<1.5)./length(DifScore);
        fprintf(1,'\b\b\b\b%3.0f%%',pairid/size(Pairs,1)*100)
    end
end
fprintf('\n')
disp(['Removing ' num2str(sum(ISIViolationsScore>0.05)) ' matched oversplits, as merging them will violate ISI >5% of the time'])
Pairs(ISIViolationsScore>0.05,:)=[];


%%
MatchProbability = arrayfun(@(X) MatchTable.MatchProb(ismember(MatchTable.ID1,Pairs(X,1)-1)&ismember(MatchTable.ID2,Pairs(X,2)-1)),1:size(Pairs,1));
[~,sortidx] = sort(MatchProbability,'descend');
Pairs = Pairs(sortidx,:);


for id = 1:size(Pairs,1)
    AllUID = find(UniqueID==UniqueID(Pairs(id,1))); %find already existing units with this UID
    if all((MatchTable.MatchProb(ismember(MatchTable.ID1,AllClusterIDs(AllUID))&ismember(MatchTable.ID2,Pairs(id,2)-1))>param.ProbabilityThreshold) | (MatchTable.MatchProb(ismember(MatchTable.ID1,Pairs(id,2)-1)&ismember(MatchTable.ID2,AllClusterIDs(AllUID)))>param.ProbabilityThreshold)') %only if all UID have a high enough probability with this second pair, we will include it to have the same UID
        UniqueID(Pairs(id,2)) = UniqueID(Pairs(id,1));
    end
end

return