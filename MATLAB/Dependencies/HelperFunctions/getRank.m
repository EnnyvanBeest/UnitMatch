function [rank, sig] = getRank(M,SessionSwitch)
    %% Find rank and rough significance

    nRec = numel(SessionSwitch)-1;
    M_noNaN = M;
    M_noNaN(isnan(M_noNaN)) = nanmin(M_noNaN(:));

   
    rank = nan(size(M));
    sig = nan(size(M));
    for did1 = 1:nRec
        for did2 = 1:nRec
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1)-1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1)-1;
            if ~all(isnan(mat2vec(M(clusIdxD1All,clusIdxD2All))))
                % Get rank
                [~,idx] = sort(M_noNaN(clusIdxD1All,clusIdxD2All),2,'descend');
                for c = 1:numel(clusIdxD1All)
                    rank(clusIdxD1All(c),clusIdxD2All(idx(c,:))) = 1:numel(clusIdxD2All);
                end

                % Find sig
                M_cut = M(clusIdxD1All,clusIdxD2All);
                sig(clusIdxD1All,clusIdxD2All) = M_cut >= nanmedian(M_cut,1) + 2*nanstd(M_cut,[],1) & ...
                    M_cut >= nanmedian(M_cut,2) + 2*nanstd(M_cut,[],2);

            end
        end
    end
    rank(isnan(M)) = nan;
    sig(isnan(M)) = nan;



end