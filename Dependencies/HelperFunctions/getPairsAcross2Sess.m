function [pairs,sessions,tableIdx,proba] = getPairsAcross2Sess(MatchTable, matchCriterion)
    %% Will get pairs across 2 sessions. Doesn't process within days.
    % MatchTable: the classic UM output
    % matchCriterion: threshold for probability (usually .5)

    if nargin < 2 || isempty(matchCriterion)
        matchCriterion = 0.5;
    end

    %% Gets bidirectional proba for each pair
    

    %% Will extract the proper pairs
    matchedPairsIdx = MatchTable.MatchProb > matchCriterion & ...
        MatchTable.RecSes1 ~= MatchTable.RecSes2; % 
    pairs = [MatchTable(matchedPairsIdx,:).ID1, MatchTable(matchedPairsIdx,:).ID2];
    sessions = [MatchTable(matchedPairsIdx,:).RecSes1, MatchTable(matchedPairsIdx,:).RecSes2];
    proba = MatchTable(matchedPairsIdx,:).MatchProb;
    tableIdx = find(matchedPairsIdx);

    [pairs,sessions,tableIdx,proba] = removeDoublets(pairs,sessions,tableIdx,proba);
end

function [pairs,sessions,tableIdx,proba] = removeDoublets(pairs,sessions,tableIdx,proba)
    % Remove doublets -- take the matching with the highest probability
    %%% Would need to deal with that better but meh for now
    %%% Couldn't find a more elegant way! :(
    for p = 1:size(pairs,1)
        [~,idxSess] = sort(sessions(p,:),'ascend');
        pairs(p,:) = pairs(p,idxSess);
        sessions(p,:) = sessions(p,idxSess);
    end
    
    % Reorder the pairs by proba
    [proba,idxSortProba] = sort(proba,'descend');
    pairs = pairs(idxSortProba,:);
    sessions = sessions(idxSortProba,:);
    tableIdx = tableIdx(idxSortProba,:);
    
    %%% SHOULD CHECK THAT EACH PAIR APPEARS BOTH WAYS?
    %%% and take average proba?
    
    sessUni = unique(sessions(:));
    for ss1 = 1:numel(sessUni)-1
        for ss2 = ss1+1:numel(sessUni) % Not looking at within session matches for now
            % Subselect pairs for these two sessions
            idxSess = find(sessions(:,1) == sessUni(ss1) & sessions(:,2) == sessUni(ss2));
    
            rmIdx = [];
            for ii = 1:numel(idxSess)
                % Decimate the pairs where one of the two is already present
                if ismember(pairs(idxSess(ii),1), pairs(idxSess(1:ii-1),1)) || ismember(pairs(idxSess(ii),2), pairs(idxSess(1:ii-1),2))
                    rmIdx = [rmIdx; idxSess(ii)];
                end
            end
            pairs(rmIdx,:) = [];
            sessions(rmIdx,:) = [];
            proba(rmIdx) = [];
            tableIdx(rmIdx) = [];
        end
    end
end