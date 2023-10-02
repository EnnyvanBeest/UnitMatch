function RankScore = getRank(M,SessionSwitch)
    %% Find rank

    nRec = numel(SessionSwitch)-1;
    M_noNaN = M;
    M_noNaN(isnan(M_noNaN)) = nanmin(M_noNaN(:));

    RankScore = nan(size(M));
    for did1 = 1:nRec
        for did2 = 1:nRec
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1)-1;
            if ~all(isnan(mat2vec(M(clusIdxD1All,clusIdxD2All))))
                [~,idx] = sort(M_noNaN(clusIdxD1All,clusIdxD2All),2,'descend');
                for c = 1:numel(clusIdxD1All)
                    RankScore(clusIdxD1All(c),clusIdxD2All(idx(c,:))) = 1:numel(clusIdxD2All);
                end
            end
        end
    end
    RankScore(isnan(M)) = nan;

end