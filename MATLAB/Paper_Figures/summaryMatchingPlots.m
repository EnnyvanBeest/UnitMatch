function [unitPresence, unitProbaMatch, days] = summaryMatchingPlots(UMFiles)

    days = cell(1, length(UMFiles));
    deltaDays = cell(1, length(UMFiles));
    deltaDaysUni = cell(1, length(UMFiles));
    unitPresence = cell(1, length(UMFiles));
    unitProbaMatch = cell(1, length(UMFiles));
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

        %% Clean up data

        % Remove the splits?

        %% For each cluster, find presence and proba of being matched in subsequent recordings

        UIDuni = unique([MatchTable.UID1]);
        days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
        days{midx} = cell2mat(days{midx}) - days{midx}{1};
        deltaDays{midx} = days{midx} - days{midx}';

        unitPresence{midx} = zeros(numel(days{midx}), numel(UIDuni));
        deltaDaysUni{midx} = unique(deltaDays{midx});
        deltaDaysUniBins = [deltaDaysUni{midx}-0.5; deltaDaysUni{midx}(end)+0.5];
        unitProbaMatch{midx} = zeros(numel(deltaDaysUni{midx}), numel(UIDuni));
        for uidx = 1:numel(UIDuni)
            sessUnitPresent = unique(MatchTable(find(MatchTable.UID1 == UIDuni(uidx)),:).RecSes1);

            % Get unit presence
            unitPresence{midx}(sessUnitPresent,uidx) = 1;

            % Get unit proba of being matched in prev and next days
            tmp = nan(numel(days{midx}),numel(days{midx}));
            tmp(sessUnitPresent,:) = unitPresence{midx}(sessUnitPresent,uidx)*unitPresence{midx}(:,uidx)';
            tmp(1 + (1+size(tmp,1))*[0:size(tmp,2)-1]) = nan; % remove diagonal
            unitProbaMatch{midx}(:,uidx) = histcounts(deltaDays{midx}(tmp == 1),deltaDaysUniBins)./ ...
                histcounts(deltaDays{midx}(ismember(tmp, [0 1])),deltaDaysUniBins);
        end

        %% Plots

        % Units lifetime 
        figure;
        imagesc(unitPresence{midx}')
        c = colormap('gray'); c = flipud(c); colormap(c)
        caxis([0 1])
        ylabel('Unit')
        xlabel('Days')
        xticks(1:numel(days{midx}));
        xticklabels(num2str(days{midx}'))

        % Probe of matching a unit 
        figure;
        plot(deltaDaysUni{midx},nanmean(unitProbaMatch{midx},2),'k')
        ylabel('P(match)')
        xlabel('Delta days')
        
        figure;
        [~,sortIdx] = sort(nanmean(unitProbaMatch{midx},1),'descend');
        imagesc(deltaDaysUni{midx},1:numel(sortIdx),unitProbaMatch{midx}(:,sortIdx)')
        colormap(c)
        hcb = colorbar; hcb.Title.String = 'P(match)';
        caxis([0 1])
        ylabel('Unit')
        xlabel('Delta days')
    end

    %% Summary plots