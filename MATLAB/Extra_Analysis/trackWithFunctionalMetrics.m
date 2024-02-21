function trackWithFunctionalMetrics(UMFiles)

    for midx = 1:length(UMFiles)
        %% Load data
    
        fprintf('Reference %s...\n', UMFiles{midx})
        
        tmpfile = dir(UMFiles{midx});
        if isempty(tmpfile)
            continue
        end
        
        fprintf('Loading the data...\n')
        tic
        load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam', 'UniqueIDConversion');
        toc

        %% Match units with functional measures 

        if any(strcmp('ISISig',MatchTable.Properties.VariableNames))
            pairsISI = MatchTable.ISISig > 2;
        else
            pairsISI = true(size(MatchTable,1),1);
        end
        if any(strcmp('natImRespSig',MatchTable.Properties.VariableNames))
            pairsNatImResp = MatchTable.natImRespSig > 2;
            if ~any(pairsNatImResp)
                pairsNatImResp = true(size(MatchTable,1),1);
            end
        else
            pairsNatImResp = true(size(MatchTable,1),1);
        end
        if any(strcmp('refPopSig',MatchTable.Properties.VariableNames))
            pairsRefPop = MatchTable.refPopSig > 2;
        else
            pairsRefPop = true(size(MatchTable,1),1);
        end

        pairsAll = pairsISI & pairsNatImResp & pairsRefPop;

        funcMatchIdx = pairsAll & (MatchTable.RecSes1 ~= MatchTable.RecSes2) & MatchTable.CentroidDist > 0.5;
        pairsAllTable = MatchTable(funcMatchIdx,:);

        %% Plot distributions

        failedMatch = pairsAllTable.UID1 ~= pairsAllTable.UID2;
        scoreBins = 0:0.05:1;
        scoreBinsCenter = scoreBins(1:end-1)+0.025;

        figure;
        for mm = 1:numel(UMparam.Scores2Include)
            subplot(1,numel(UMparam.Scores2Include),mm); hold all
            hFailed = histcounts(pairsAllTable(failedMatch,:).(UMparam.Scores2Include{mm}),scoreBins,'Normalization','probability');
            hFound = histcounts(pairsAllTable(~failedMatch,:).(UMparam.Scores2Include{mm}),scoreBins,'Normalization','probability');
            hNonMatch = histcounts(MatchTable(~funcMatchIdx & MatchTable.UID1 ~= MatchTable.UID2,:).(UMparam.Scores2Include{mm}),scoreBins,'Normalization','probability');
            plot(scoreBinsCenter,hFailed,'k')
            plot(scoreBinsCenter,hFound,'r')
            plot(scoreBinsCenter,hNonMatch,'Color',[.5 .5 .5])
            title(UMparam.Scores2Include{mm})
            xlabel('Sim score')
            ylabel('Proportion')
        end

    end

