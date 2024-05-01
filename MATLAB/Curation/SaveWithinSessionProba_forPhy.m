function SaveWithinSessionProba_forPhy(SaveDir, recompute)
    %%% Will save the within session probability matrix for all the
    %%% recordings involved. This can be used by phy.

    if ~exist('recompute','var')
        recompute = 0;
    end

    load(fullfile(SaveDir,'UnitMatch.mat'), 'MatchTable', 'UniqueIDConversion', 'UMparam')
    recSesAll = unique(MatchTable.RecSes1);
    
    for rr = 1:numel(recSesAll)
        targetFile = fullfile(UMparam.KSDir{recSesAll(rr)},'probability_templates.npy');

        if exist(targetFile,'file') && ~recompute
            continue
        end

        idxTable = MatchTable.RecSes1 == recSesAll(rr) & MatchTable.RecSes2 == recSesAll(rr);
        idxClusters = ismember(UniqueIDConversion.OriginalClusID(UniqueIDConversion.recsesAll == recSesAll(rr)), unique(MatchTable(idxTable,:).ID1));
        probaMatrix = nan(numel(idxClusters), numel(idxClusters));
        probaMatrix(idxClusters,idxClusters) = reshape(MatchTable(idxTable,:).MatchProb, ...
            [sum(idxClusters) sum(idxClusters)]);

        % easier if symmetrical
        probaMatrix = .5*probaMatrix + .5*probaMatrix';

        % write it
        writeNPY(probaMatrix, targetFile);
    end
end