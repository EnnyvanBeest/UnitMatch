function summaryFunctionalPlots_Part2(UMFiles, groupVector, UseKSLabels)
    %% Will plot summary plots: distribution, ROC and AUC. 
    % UMFiles: list cells contains path to UnitMatch.m files
    % whichMetric: will compute distributions/ROC/AUC on either 'Corr', 'Rank', or 'Sig'. 
    % groupVector: same size as UMFiles, to group e.g. recordings from the
    % same mouse together.
    % UseKSLabels: use KS labels for matching

    %% Define parameters
    if ~exist('groupVector','var')
        groupVector = 1:length(UMFiles);
    end
    groups = unique(groupVector);
    groupColor = gray(length(groups)+1);

    if ~exist('UseKSLabels','var')
        UseKSLabels = 0;
    end

    ROCBins = 0:0.01:1;
    minMatches = 20;
    durLim = 10*60;
    minUnits = 25;
    WithinSesNoiseLim = 0.20; % if more than 20% of within session matches are made this shouldn't count

    
    %% Loop over mice to get all Distributions & q params
    deltaDays = cell(1, length(UMFiles));
    numMatchedUnits = cell(1, length(UMFiles));
    InitialDrift = cell(1,length(UMFiles));
    FixedDrift =  cell(1,length(UMFiles));
    maxAvailableUnits = cell(1, length(UMFiles));
    for midx = 1:length(UMFiles)
        %% Load data
    
        fprintf('Reference %s...\n', UMFiles{midx})
    
        tmpfile = dir(UMFiles{midx});
        if isempty(tmpfile)
            continue
        end
    
        fprintf('Loading the data...\n')
        tic
        load(fullfile(tmpfile.folder, tmpfile.name), 'MatchTable', 'UMparam','UniqueIDConversion');
        toc
    
        sessIDs = unique(MatchTable.RecSes1);
           
        %%% HACK -- Can remove later
        if ~isfield(UMparam,'RawDataPaths')
            UMparam.RawDataPaths = UMparam.AllRawPaths;
            rmfield(UMparam,'AllRawPaths')
            save(fullfile(tmpfile.folder, tmpfile.name), 'UMparam','-append')
        end
        if ~iscell(UMparam.RawDataPaths)
            for ii = 1:numel(UMparam.RawDataPaths)
                tmp{ii} = UMparam.RawDataPaths(ii);
            end
            UMparam.RawDataPaths = tmp;
        end
    
        %% Loop through pairs of sessions
        WithinSesNoise = 0;
        fprintf('Looping through days...\n')
        tic
        days{midx} = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
        days{midx} = cell2mat(days{midx}) - days{midx}{1};
        deltaDays{midx} = nan(numel(sessIDs),numel(sessIDs));
        numMatchedUnits{midx} = nan(numel(sessIDs),numel(sessIDs));
        maxAvailableUnits{midx} = nan(numel(sessIDs),numel(sessIDs));
        InitialDrift{midx} = nan(numel(sessIDs),numel(sessIDs));
        FixedDrift{midx} = nan(numel(sessIDs),numel(sessIDs));
        nUnitsSelected{midx} = arrayfun(@(X) sum(UniqueIDConversion.GoodID(UniqueIDConversion.recsesAll==X)),unique(UniqueIDConversion.recsesAll));
        nUnitsTotal{midx} = arrayfun(@(X) length(UniqueIDConversion.GoodID(UniqueIDConversion.recsesAll==X)),unique(UniqueIDConversion.recsesAll));

        for sess1Idx = 1:numel(sessIDs)
    
            sess1 = sessIDs(sess1Idx);
            day1 = days{midx}(sess1Idx);
            meta = ReadMeta2(UMparam.RawDataPaths{sess1}.folder);
            durSess1 = str2double(meta.fileTimeSecs);
            if durSess1 < durLim 
                continue
            end

            for sess2Idx = sess1Idx:numel(sessIDs)
    
                sess2 = sessIDs(sess2Idx);
                day2 = days{midx}(sess2Idx);
                deltaDays{midx}(sess1Idx,sess2Idx) = abs(day2 - day1);
                meta = ReadMeta2(UMparam.RawDataPaths{sess2}.folder);
                durSess2 = str2double(meta.fileTimeSecs);
                if durSess2 < durLim
                    continue
                end

                %% Cut table to specific days

                MatchTable_2sess = MatchTable(ismember(MatchTable.RecSes1, [sess1 sess2]) & ismember(MatchTable.RecSes2, [sess1 sess2]), :);
                if min([length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess1))) length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess2)))]) < minUnits
                    continue
                end
                %% Number of matches
    
                if ~UseKSLabels
                    %%% CHECK THAT THIS MAKES SENSE
                    %%% CHOOSE BASED ON UID
                    if sess1 ~= sess2
                        matchedUnitsIdx = (MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2); % using Unique ID
                    else
                        matchedUnitsIdx = (MatchTable_2sess.UID1 == MatchTable_2sess.UID2) & (MatchTable_2sess.ID1 ~= MatchTable_2sess.ID2); % using Unique ID
                    end
                    %%% OR RECOMPUTE
%                     [~,~,idx,~] = getPairsAcross2Sess(MatchTable_2sess, UMparam.ProbabilityThreshold);
%                     matchedUnitsIdx = zeros(size(MatchTable_2sess,1),1);
%                     matchedUnitsIdx(idx) = 1;
                else
                    matchedUnitsIdx = (MatchTable_2sess.ID1 == MatchTable_2sess.ID2) & (MatchTable_2sess.RecSes1 ~= MatchTable_2sess.RecSes2);
                end
                numMatchedUnits{midx}(sess1Idx,sess2Idx) = sum(matchedUnitsIdx)/2; % Divided by two because looking both ways -- can be non-integer
                
                maxAvailableUnits{midx}(sess1Idx,sess2Idx) = min([length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess1))) length(unique(MatchTable_2sess.ID1(MatchTable_2sess.RecSes1 == sess2)))]);%

                if sess1Idx == sess2Idx  && numMatchedUnits{midx}(sess1Idx,sess2Idx)/maxAvailableUnits{midx}(sess1Idx,sess2Idx) > WithinSesNoiseLim
                    warning([tmpfile.folder ', Sess' num2str(sess1Idx) 'x Sess' num2str(sess2Idx) ', ' num2str(round(numMatchedUnits{midx}(sess1Idx,sess2Idx)/maxAvailableUnits{midx}(sess1Idx,sess2Idx).*100)) '%'])
                    WithinSesNoise = 1;
                end
                
                %% Extract drift if present
    
                if isfield(UMparam,'drift') && sess2Idx>1
                    InitialDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,1),2); % Drift in recording 1 is 1 vs 2, etc.
                    FixedDrift{midx}(sess1Idx,sess2Idx) =  vecnorm(UMparam.drift(sess2Idx-1,:,2),2); % Drift in recording 1 is 1 vs 2, etc.
                end
            end
        end
        if WithinSesNoise % This file shouldn't count, too much within session noise
            numMatchedUnits{midx} = nan(numel(sessIDs),numel(sessIDs));
            maxAvailableUnits{midx} = nan(numel(sessIDs),numel(sessIDs));
            InitialDrift{midx} = nan(numel(sessIDs),numel(sessIDs));
            FixedDrift{midx} = nan(numel(sessIDs),numel(sessIDs));
        end
        %% qParams
        if ~exist(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))
            try
                QualityMetricsROCs(UMparam.SaveDir)
                close all
            catch ME
                keyboard
            end
        end
        load(fullfile(tmpfile.folder, 'qMetricAUCs.mat'))

        if ~exist('qParamNames')
            qParamNames = AUCqParams.qParamNames;
            AUCPerqParam = nan(length(qParamNames),0);
        end
        [takethese,puthere] = ismember(AUCqParams.qParamNames,qParamNames);
        tmpQP = nan(length(qParamNames),1);
        tmpQP(puthere(puthere~=0)) = AUCqParams.AUCMvNM(takethese);
        AUCPerqParam = cat(2,AUCPerqParam,tmpQP);
        toc
    end

     SelUnitsPerc = cat(1,nUnitsSelected{:})./cat(1,nUnitsTotal{:}).*100;
    %% AUC distributions for qMetrics
    figure('name','AUC Distr')
    stepsz = 0.05;
    binvec = [0:stepsz:1];
    plotvec = stepsz/2:stepsz:1-stepsz/2;
    for qid = 1:length(qParamNames)
        subplot(ceil(sqrt(length(qParamNames))),round(sqrt(length(qParamNames))),qid)
        if nanmedian(AUCPerqParam(qid,:))<0.5
            AUCPerqParam(qid,:) = 1-AUCPerqParam(qid,:);
        end
        nums = histcounts(AUCPerqParam(qid,:),binvec);
        plot(plotvec,nums,'k-');
        hold on
        line([nanmedian(AUCPerqParam(qid,:)) nanmedian(AUCPerqParam(qid,:))],get(gca,'ylim'),'color',[1 0 0])
        [h,p(qid)] = ttest(AUCPerqParam(qid,:),0.5);

        title([qParamNames{qid} ', p=' num2str(round(p(qid)*100)./100)])
        xlim([0 1])
        makepretty
        offsetAxes
    end

    
    %% Additional figure 
    % Plot number of matches as a function of delta days
    figure;
    hold all
    for midx = 1:length(UMFiles)
        scatter(deltaDays{midx}(:), mat2vec(numMatchedUnits{midx}(:))./mat2vec(maxAvailableUnits{midx}(:)),10,groupColor(groupVector(midx),:),'filled')
    end
    ylabel('Number of matches')
    xlabel('\Deltadays')
    

    %% Dependence of number matched units on drift
    figure('name','NrUnits versus drift')
    for midx = 1:length(UMFiles)
        subplot(1,2,1)
        scatter(InitialDrift{midx},numMatchedUnits{midx}./maxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
        hold on
        xlabel('Drift (Eucl Distance)')
        ylabel('Proportion of matches')
        xlim([0 UMparam.NeighbourDist])
        title('Initial Drift')
    
        subplot(1,2,2)
        scatter(FixedDrift{midx},numMatchedUnits{midx}./ maxAvailableUnits{midx},20,groupColor(groupVector(midx),:),'filled')
        hold on
        xlabel('Drift (Eucl Distance)')
        ylabel('Proportion of matches')
        xlim([0 UMparam.NeighbourDist])
    end
    InitialDrift(cellfun(@isempty,InitialDrift)) = [];
    InitialDrift = cellfun(@(X) X(:),InitialDrift,'uni',0);
    maxAvailableUnits(cellfun(@isempty,maxAvailableUnits)) = [];
    maxAvailableUnitstmp = cellfun(@(X) X(:),maxAvailableUnits,'uni',0);
    numMatchedUnits(cellfun(@isempty,numMatchedUnits)) = [];
    numMatchedUnitstmp = cellfun(@(X) X(:),numMatchedUnits,'uni',0);
    deltaDays(cellfun(@isempty,deltaDays)) = [];
    deltaDaystmp = cellfun(@(X) X(:),deltaDays,'uni',0);

    InitialDrift = cat(1,InitialDrift{:});

    maxAvailableUnitstmp = cat(1,maxAvailableUnitstmp{:});
    numMatchedUnitstmp = cat(1,numMatchedUnitstmp{:});
    deltaDaystmp = cat(1,deltaDaystmp{:});
    FixedDrift(cellfun(@isempty,FixedDrift)) = [];
    FixedDrift = cellfun(@(X) X(:),FixedDrift,'uni',0);
    FixedDrift = cat(1,FixedDrift{:});
    
    
    FixedDrift(FixedDrift>UMparam.maxdist) = nan;
    InitialDrift(InitialDrift>UMparam.maxdist) = nan;
    
    PercNeurons =     numMatchedUnitstmp./maxAvailableUnitstmp;
    [r,p] = corr(InitialDrift(~isnan(InitialDrift)),PercNeurons(~isnan(InitialDrift)));
    subplot(1,2,1)
    title(['Initial Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])
    
    subplot(1,2,2)
    [r,p] = corr(FixedDrift(~isnan(FixedDrift)),PercNeurons(~isnan(FixedDrift)));
    title(['Corrected Drift, r=' num2str(round(r*100)/100) ', p=' num2str(round(p*100)/100)])

    %% Histogram
    TakeThese = find(deltaDaystmp == 2)% Only take days close to each other (1 day apart)
    PercNeurons = numMatchedUnitstmp(TakeThese)./maxAvailableUnitstmp(TakeThese);

    stepsz = 5;
    binvec = 0:stepsz:100;
    plotvec = stepsz/2:stepsz:100-stepsz/2;
    hc = histcounts(PercNeurons.*100,binvec)./sum(~isnan(PercNeurons)).*100;
    figure('name','Matching performance')
    plot(plotvec,hc)
    xlabel('Matched pairs / maximally available (%)')
    ylabel('pairs of days (%)')
    makepretty
    offsetAxes
    save(fullfile('C:\Users\EnnyB\OneDrive - University College London\UnitMatch_Manuscript','chronicMice.mat'),'PercNeurons','numMatchedUnits','maxAvailableUnits')
    disp([num2str(nanmedian(PercNeurons.*100)) '+/- ' num2str(mad(PercNeurons.*100)) ' in ' num2str(sum(~isnan(PercNeurons))) ' recordings'])
     
    %% 
    tmpAcute = load(fullfile('C:\Users\EnnyB\OneDrive - University College London\UnitMatch_Manuscript','AcuteMice.mat')) % Load other one
    disp([num2str(nanmedian(tmpAcute.PercNeurons.*100)) '+/- ' num2str(mad(tmpAcute.PercNeurons.*100)) ' in ' num2str(sum(~isnan(tmpAcute.PercNeurons))) ' recordings'])

    figure('name','Acute vs Chronic')
    ha = histcounts(tmpAcute.PercNeurons.*100,binvec)./sum(~isnan(tmpAcute.PercNeurons)).*100;
    plot(plotvec,ha,'color',[0.5 0.5 0.5])
    hold on
        plot(plotvec,hc,'color',[0 0 0])
         xlabel('Matched pairs / maximally available (%)')
    ylabel('pairs of days (%)')
    makepretty
    offsetAxes
    legend('Acute','Chronic')
    w = ranksum(tmpAcute.PercNeurons,PercNeurons)


end