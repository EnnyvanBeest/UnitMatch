function ComputeFunctionalScores(SaveDir, saveFig, recompute)

if nargin < 2 || isempty(saveFig)
    saveFig = 0;
end

if nargin < 3 || isempty(recompute)
    recompute = 0;
end

load(fullfile(SaveDir, 'UnitMatch.mat'), 'MatchTable', 'UMparam', 'UniqueIDConversion');
UMparam.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows

if all(ismember({'refPopCorr','ACGCorr','FRDiff','natImRespCorr'},MatchTable.Properties.VariableNames)) && ~recompute
    disp('Already computed functional scores')
    return
end

% Extract cluster information
if UMparam.GoodUnitsOnly
    GoodId = logical(UniqueIDConversion.GoodID);
else
    GoodId = true(1, length(UniqueIDConversion.GoodID));
end
UniqueID = UniqueIDConversion.UniqueID(GoodId);
OriID = UniqueIDConversion.OriginalClusID(GoodId);
OriIDAll = UniqueIDConversion.OriginalClusID;
recses = UniqueIDConversion.recsesAll(GoodId);
recsesall = UniqueIDConversion.recsesAll;

AllKSDir = UMparam.KSDir; %original KS Dir
nclus = length(UniqueID);

%%% TO REMOVE
if ~isfield(UMparam,'AllRawPaths') % For now no to disturb Célian
    UMparam.AllRawPaths = UMparam.RawDataPaths;
end


if length(UMparam.AllRawPaths{1}) > 1 %Reshape for Stitched
    UMparam.AllRawPaths = arrayfun(@(X) UMparam.AllRawPaths{1}(X),1:length(UMparam.AllRawPaths{1}),'uni',0);
end

% Load SP
disp('Loading spike information...')
sp = getSpikesFromPrepData(AllKSDir);

nRec = length(unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID))));
RecOpt = unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID)));
EuclDist = reshape(MatchTable.EucledianDistance, nclus, nclus);
SessionSwitch = arrayfun(@(X) find(recses == X, 1, 'first'), unique(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID))), 'Uni', 0);
SessionSwitch(cellfun(@isempty, SessionSwitch)) = [];
SessionSwitch = [cell2mat(SessionSwitch); nclus + 1];

% Plotting order (sort units based on distance)
[~, SortingOrder] = arrayfun(@(X) sort(EuclDist(1, SessionSwitch(X):SessionSwitch(X+1)-1)), 1:nRec, 'Uni', 0);
SortingOrder = arrayfun(@(X) squeeze(SortingOrder{X}+SessionSwitch(X)-1), 1:nRec, 'Uni', 0);
if size(SortingOrder{1}, 1) == 1
    SortingOrder = cat(2, SortingOrder{:});
else
    SortingOrder = cat(1, SortingOrder{:});
end

%% Get cross-correlation fingerprints correlations

if ~any(ismember(MatchTable.Properties.VariableNames, 'refPopCorr')) || recompute% If it already exists in table, skip this entire thing

    %% Compute cross-correlation matrices for individual recordings
    disp('Computing cross-correlation fingerprint')
    SessionCorrelations = cell(1, nRec);
    for rid = 1:nRec
        % Define edges for this dataset
        edges = floor(min(sp.st(sp.RecSes == RecOpt(rid)))) - UMparam.binsz / 2:UMparam.binsz:ceil(max(sp.st(sp.RecSes == RecOpt(rid)))) + UMparam.binsz / 2;
        Good_Idx = find(GoodId & recsesall' == RecOpt(rid)); % Only care about good units at this point

        % bin data to create PSTH
        sr = nan(numel(Good_Idx), numel(edges)-1);
        for uid = 1:numel(Good_Idx)
            sr(uid, :) = histcounts(sp.st(sp.spikeTemplates == OriIDAll(Good_Idx(uid)) & sp.RecSes == recsesall(Good_Idx(uid))), edges);
        end

        % Define folds (two halves)
        idx_fold1 = 1:floor(size(sr, 2)./2);
        idx_fold2 = floor(size(sr, 2)./2) + 1:floor(size(sr, 2)./2) * 2;

        % Find cross-correlation in first and second half of session
        SessionCorrelations{rid}.fold1 = corr(sr(:, idx_fold1)', sr(:, idx_fold1)')';
        SessionCorrelations{rid}.fold2 = corr(sr(:, idx_fold2)', sr(:, idx_fold2)')';

        % Nan the diagonal
        SessionCorrelations{rid}.fold1(logical(eye(size(SessionCorrelations{rid}.fold1)))) = nan;
        SessionCorrelations{rid}.fold2(logical(eye(size(SessionCorrelations{rid}.fold2)))) = nan;
    end

    %% Get correlation matrics for fingerprint correlations
    sessionCorrelationsAll = cell(1, nRec);
    for did = 1:nRec
        if length(SessionCorrelations) == nRec
            sessionCorrelationsAll{did} = SessionCorrelations{did};
        elseif length(SessionCorrelations) == 1 %Normal situation
            if iscell(SessionCorrelations)
                sessionCorrelationsAll{did} = SessionCorrelations{1};
            else
                sessionCorrelationsAll{did} = SessionCorrelations;
            end
        else
            disp('This is a weird situation...')
            keyboard
        end
    end

    %% Compute functional score (cross-correlation fingerprint)
    if nRec < 5
        drawdrosscorr = 1;
    else
        drawdrosscorr = 0;
    end
    MatchProbability = reshape(MatchTable.MatchProb, nclus, nclus);
    [r, c] = find(MatchProbability > UMparam.ProbabilityThreshold); %Find matches
    Pairs = cat(2, r, c);
    Pairs = sortrows(Pairs);
    Pairs = unique(Pairs, 'rows');
    [refPopCorr,~,AllSessionCorrelations] = CrossCorrelationFingerPrint(sessionCorrelationsAll, Pairs, OriID, recses, drawdrosscorr);

    % get rank
    [refPopRank, refPopSig] = getRank(refPopCorr,SessionSwitch);
    refPopSig(refPopCorr==0) = nan; % Correlation of 0 means nothing
    refPopRank(refPopCorr==0) = nan; % Correlation of 0 means nothing
    
    % Output needs transposing to be properly stored in table
    refPopCorr = refPopCorr';
    refPopRank = refPopRank'; 
    refPopSig = refPopSig';

    % Save in table
    MatchTable.refPopCorr = refPopCorr(:);
    MatchTable.refPopRank = refPopRank(:); % What goes in the table should give ndays for every output when you do sum(refPopRank==1,1), if it's not, transpose!
    MatchTable.refPopSig = refPopSig(:);

    %% Compare to functional scores
    if saveFig % Otherwise this takes way tooo long
        figure;

        subplot(1, 3, 1)
        imagesc(refPopRank(SortingOrder, SortingOrder) == 1 & refPopSig(SortingOrder, SortingOrder) == 1)
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        colormap(flipud(gray))
        title('rank == 1*')
        makepretty

        subplot(1, 3, 2)
        imagesc(MatchProbability(SortingOrder, SortingOrder) > UMparam.ProbabilityThreshold)
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        colormap(flipud(gray))
        title(['Match Probability>', num2str(UMparam.ProbabilityThreshold)])
        makepretty

        subplot(1, 3, 3)
        imagesc(MatchProbability(SortingOrder, SortingOrder) >= UMparam.ProbabilityThreshold | (MatchProbability(SortingOrder, SortingOrder) > 0.05 & refPopRank(SortingOrder, SortingOrder) == 1 & refPopSig(SortingOrder, SortingOrder) == 1));
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        colormap(flipud(gray))
        title('Matching probability + rank')
        makepretty
        if saveFig
            saveas(gcf, fullfile(SaveDir, 'RankVSProbability.fig'))
            saveas(gcf, fullfile(SaveDir, 'RankVSProbability.bmp'))
        end

        tmpf = triu(refPopCorr);
        tmpm = triu(MatchProbability);
        tmpr = triu(refPopRank);
        tmpr = tmpr(tmpf ~= 0);
        tmpm = tmpm(tmpf ~= 0);
        tmpf = tmpf(tmpf ~= 0);
        figure;
        scatter(tmpm, tmpf, 14, tmpr, 'filled')
        colormap(cat(1, [0, 0, 0], winter))
        xlabel('Match Probability')
        ylabel('Cross-correlation fingerprint')
        makepretty
        saveas(gcf, fullfile(SaveDir, 'RankVSProbabilityScatter.fig'))
        saveas(gcf, fullfile(SaveDir, 'RankVSProbabilityScatter.bmp'))


        % Check these: should be z1, z2, z3, z12, z13, z23, z123
        VenFig = figure('name','Venn, r=Rank, b=sig, g=Match');
        subplot(2,2,1)
        Idx = MatchTable.refPopRank(:) == 1 | MatchTable.refPopSig(:) == 1 | MatchTable.MatchProb(:) > 0.5;
        h = venn([sum(MatchTable.refPopRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.refPopSig(Idx)==0) sum(MatchTable.refPopRank(Idx)>1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.refPopSig(Idx)==1) ...
            sum(MatchTable.refPopRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.refPopSig(Idx)==0) sum(MatchTable.refPopRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.refPopSig(Idx)==1) ...
            sum(MatchTable.refPopRank(Idx)==1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.refPopSig(Idx)==0) sum(MatchTable.refPopRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.refPopSig(Idx)==1) ...
            sum(MatchTable.refPopRank(Idx)==1 & MatchTable.refPopSig(Idx)==1 & MatchTable.MatchProb(Idx)>0.5)] );
        axis square
        axis off
        makepretty
        title('RefPopCor')

    end
end

%% Get ACG fingerprints correlations

if ~any(ismember(MatchTable.Properties.VariableNames, 'ACGCorr')) || recompute % If it already exists in table, skip this entire thing
    %% Compute ACG and correlate them between units
    % This is very time consuming
    disp('Computing ACG, this will take some time...')
    tvec = -UMparam.ACGduration / 2:UMparam.ACGbinSize:UMparam.ACGduration / 2;
    ACGMat = nan(length(tvec), 2, nclus);
    FR = nan(2, nclus);
    for clusid = 1:nclus %parfot QQ
        for cv = 1:2
            idx1 = find(sp.spikeTemplates == OriID(clusid) & sp.RecSes == recses(clusid));
            if ~isempty(idx1)
                if cv == 1
                    idx1 = idx1(1:floor(length(idx1)/2));
                else
                    idx1 = idx1(ceil(length(idx1)/2):end);
                end

                % Compute Firing rate
                nspkspersec = histcounts(sp.st(idx1), [min(sp.st(idx1)):1:max(sp.st(idx1))]);
                FR(cv, clusid) = nanmean(nspkspersec);


                % compute ACG
                [ccg, t] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
                    ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', UMparam.ACGbinSize, 'duration', UMparam.ACGduration, 'norm', 'rate'); %function
                ACGMat(:, cv, clusid) = ccg(:, 1, 1);
            end
        end
    end

    %% Correlation between ACG
    ACGCorr = corr(squeeze(ACGMat(:, 1, :)), squeeze(ACGMat(:, 2, :)));
    ACGCorr = tanh(.5*atanh(ACGCorr) + .5*atanh(ACGCorr)); %%% added after biorxiv
    ACGCorr = ACGCorr'; % getRank expects different input

    [ACGRank, ACGSig] = getRank(ACGCorr,SessionSwitch);

    % Transpose
    ACGCorr = ACGCorr';
    ACGRank = ACGRank';
    ACGSig = ACGSig';

    % Save in table
    MatchTable.ACGCorr = ACGCorr(:);
    MatchTable.ACGRank = ACGRank(:);
    MatchTable.ACGSig = ACGSig(:);
end

if saveFig
    % Check these: should be z1, z2, z3, z12, z13, z23, z123
    figure(VenFig)
    subplot(2,2,2)
    Idx = MatchTable.ACGRank(:) == 1 | MatchTable.ACGSig(:) == 1 | MatchTable.MatchProb(:) > 0.5;
    h = venn([sum(MatchTable.ACGRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.ACGSig(Idx)==0) sum(MatchTable.ACGRank(Idx)>1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.ACGSig(Idx)==1) ...
        sum(MatchTable.ACGRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.ACGSig(Idx)==0) sum(MatchTable.ACGRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.ACGSig(Idx)==1) ...
        sum(MatchTable.ACGRank(Idx)==1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.ACGSig(Idx)==0) sum(MatchTable.ACGRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.ACGSig(Idx)==1) ...
        sum(MatchTable.ACGRank(Idx)>1 & MatchTable.ACGSig(Idx)==1 & MatchTable.MatchProb(Idx)>0.5)] );
    axis square
    axis off
    makepretty
    title('ACGCor')
end

%% Get FR difference

if ~any(ismember(MatchTable.Properties.VariableNames, 'FRDiff')) || recompute
    FR = repmat(permute(FR, [2, 1]), [1, 1, nclus]);
    FRDiff = abs(squeeze(FR(:, 2, :)-permute(FR(:, 1, :), [3, 2, 1])));
    [FRRank, FRSig] = getRank(-FRDiff,SessionSwitch);

    % Transpose
    FRDiff = FRDiff';
    FRRank = FRRank';
    FRSig = FRSig';

    % Save in table
    MatchTable.FRDiff = FRDiff(:);
    MatchTable.FRRank = FRRank(:);
    MatchTable.FRSig = FRSig(:);
end
if saveFig
    % Check these: should be z1, z2, z3, z12, z13, z23, z123
    figure(VenFig)
    subplot(2,2,3)
    Idx = MatchTable.FRRank(:) == 1 | MatchTable.FRSig(:) == 1 | MatchTable.MatchProb(:) > 0.5;
    h = venn([sum(MatchTable.FRRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.FRSig(Idx)==0) sum(MatchTable.FRRank(Idx)>1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.FRSig(Idx)==1) ...
        sum(MatchTable.FRRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.FRSig(Idx)==0) sum(MatchTable.FRRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.FRSig(Idx)==1) ...
        sum(MatchTable.FRRank(Idx)==1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.FRSig(Idx)==0) sum(MatchTable.FRRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.FRSig(Idx)==1) ...
        sum(MatchTable.FRRank(Idx)==1 & MatchTable.FRSig(Idx)==1 & MatchTable.MatchProb(Idx)>0.5)] );
    axis square
    axis off
    makepretty
    title('FRDiff')
   
end


%% Get natural images fingerprints correlations

if ~any(ismember(MatchTable.Properties.VariableNames, 'natImRespCorr')) || recompute % If it already exists in table, skip this entire thing
    % Param for processing
    proc.window = [-0.3 0.5 ... % around onset
        0.0 0.5]; % around offset
    proc.binSize = 0.005; % in ms
    proc.smoothSize = 5; % PSTH smoothing filter
    gw = gausswin(proc.smoothSize,3);
    proc.smWin = gw./sum(gw);

    binsOn = (proc.window(1)+proc.binSize/2):proc.binSize:proc.window(2);
    binsOff =  proc.window(2)+((proc.binSize/2):proc.binSize:(proc.window(4)-proc.window(3)));
    bins = [binsOn binsOff];

    nRec = length(RecOpt); % not always the same!! numel(UMparam.AllRawPaths);
    nClu = nan(1,nRec); for ss = 1:nRec; nClu(ss) = numel(unique(MatchTable.ID1(MatchTable.RecSes1 == RecOpt(ss)))); end
    spikeData_cv = cell(2,nRec);
    for ss = 1:nRec
        % Get the original binFile (also for stitched?)
        Flag = 1;
        if ~isempty(UMparam.AllRawPaths{RecOpt(ss)}) % When no raw data is available
            if iscell(UMparam.AllRawPaths{RecOpt(ss)})
                binFileRef = fullfile(UMparam.AllRawPaths{RecOpt(ss)});
            else
                binFileRef = fullfile(UMparam.AllRawPaths{RecOpt(ss)}.folder,UMparam.AllRawPaths{RecOpt(ss)}.name);
            end
        else
            Flag = 0;
        end

        % Find the associated experiments
        if ispc && Flag
            exp2keep = getNatImExpRef(binFileRef);
        else
            exp2keep = [];
            % Julie Fabre: on linux the servers are mounted differently, and the above
            % function needs to be substantially changed to work. Since I
            % don't need to run it, just added an "ispc" flag.
        end

        if ~isempty(exp2keep)
            % Get the spikes
            st = sp.st(sp.RecSes == RecOpt(ss));
            clu = sp.spikeTemplates(sp.RecSes == RecOpt(ss));

            spikesAll.times = st;
            spikesAll.clusters = clu;
            spikesAll.clusterIDs = unique(clu); % follow the same units across days

            % Get the natIm responses
            spikeData = getNatImResp(spikesAll,exp2keep,binFileRef,proc);
            clusterIDs = spikesAll.clusterIDs;

            % Split in two halves and subselect units
            cluIdx = ismember(clusterIDs,unique(MatchTable.ID1(MatchTable.RecSes1 == RecOpt(ss))));
            spikeData_cv{1,ss} = spikeData(:,:,cluIdx,1:2:end); % odd
            spikeData_cv{2,ss} = spikeData(:,:,cluIdx,2:2:end); % even
        end
    end


    % ---
    % Second fingerprint with mean response to NI + temporal responses
    corrMeanResp_big = nan(sum(nClu),sum(nClu));
    corrTimecourse_big = nan(sum(nClu),sum(nClu));
    for ss1 = 1:nRec
        if ~isempty(spikeData_cv{1,ss1})
            meanResp1 = zscore(squeeze(nanmean(spikeData_cv{1,ss1}(:,bins<proc.window(2)+0.2 & bins>0,:,:),[2 4])));
            timecourse1 = zscore(squeeze(nanmean(spikeData_cv{1,ss1}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),[1 4])));
            for ss2 = 1:nRec
                if ~isempty(spikeData_cv{2,ss2})
                    meanResp2 = zscore(squeeze(nanmean(spikeData_cv{2,ss2}(:,bins<proc.window(2)+0.2 & bins>0,:,:),[2 4])));
                    timecourse2 = zscore(squeeze(nanmean(spikeData_cv{2,ss2}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),[1 4])));
                    corrMeanResp_big(sum(nClu(1:ss1-1))+1:sum(nClu(1:ss1)), sum(nClu(1:ss2-1))+1:sum(nClu(1:ss2))) = corr(meanResp1,meanResp2);
                    corrTimecourse_big(sum(nClu(1:ss1-1))+1:sum(nClu(1:ss1)), sum(nClu(1:ss2-1))+1:sum(nClu(1:ss2))) = corr(timecourse1, timecourse2);
                end
            end
        end
    end
    corrResp_big = tanh(.5*atanh(corrMeanResp_big) + .5*atanh(corrTimecourse_big));

    corrResp_big = tanh(.5*atanh(corrResp_big) + .5*atanh(corrResp_big')); % not sure that's needed?

    % Get rank
    [natImRespRank, natImRespSig] = getRank(corrResp_big, SessionSwitch);

    % Output needs transposing to be properly stored in table
    corrResp_big = corrResp_big';
    natImRespRank = natImRespRank';
    natImRespSig = natImRespSig';

    % Save in table
    MatchTable.natImRespCorr = corrResp_big(:);
    MatchTable.natImRespRank = natImRespRank(:); % What goes in the table should give ndays for every output when you do sum(refPopRank==1,1), if it's not, transpose!
    MatchTable.natImRespSig = natImRespSig(:);

    if saveFig & ~isempty(exp2keep)
        % Check these: should be z1, z2, z3, z12, z13, z23, z123
        figure(VenFig)
        subplot(2,2,3)
        Idx = MatchTable.natImRespRank(:) == 1 | MatchTable.natImRespSig(:) == 1 | MatchTable.MatchProb(:) > 0.5;
        h = venn([sum(MatchTable.natImRespRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.natImRespSig(Idx)==0) sum(MatchTable.natImRespRank(Idx)>1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.natImRespSig(Idx)==1) ...
            sum(MatchTable.natImRespRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.natImRespSig(Idx)==0) sum(MatchTable.natImRespRank(Idx)==1 & MatchTable.MatchProb(Idx)<=0.5 & MatchTable.natImRespSig(Idx)==1) ...
            sum(MatchTable.natImRespRank(Idx)==1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.natImRespSig(Idx)==0) sum(MatchTable.natImRespRank(Idx)>1 & MatchTable.MatchProb(Idx)>0.5 & MatchTable.natImRespSig(Idx)==1) ...
            sum(MatchTable.natImRespRank(Idx)==1 & MatchTable.natImRespSig(Idx)==1 & MatchTable.MatchProb(Idx)>0.5)] );
        axis square
        axis off
        makepretty
        title('NatImg')



    end


%     % ---
%     % Second fingerprint with CCA
%     
%     % Perform CCA across recordings
%     [corrMat, ~] = computeNatImCorr(spikeData_cv(:)); % default 75
%     % [corrMat, ~] = computeNatImCorr(spikeData_cv(:), 1:ceil(sqrt(size(spikeData_cv{1},1)*size(spikeData_cv{1},2))/2));
% 
%     % Reshape the matrix to a single one with correct clusters
%     corrWCCA_big = nan(sum(nClu),sum(nClu));
%     corrWCCA_1x2 = corrMat(1:2:end, 2:2:end);
%     for ss1 = 1:nRec
%         for ss2 = 1:nRec
%             if ~all(isnan(corrWCCA_1x2{ss1,ss2}(:)))
%                 corrWCCA_big(sum(nClu(1:ss1-1))+1:sum(nClu(1:ss1)), sum(nClu(1:ss2-1))+1:sum(nClu(1:ss2))) = corrWCCA_1x2{ss1,ss2};
%             end
%         end
%     end
%     corrWCCA_big = tanh(.5*atanh(corrWCCA_big) + .5*atanh(corrWCCA_big')); % not sure that's needed?
% 
%     % Get rank
%     [natImRank, natImSig] = getRank(corrWCCA_big, SessionSwitch);
%     
%     % Save in table
%     MatchTable.natImCorr = corrWCCA_big(:);
%     MatchTable.natImRank = natImRank(:);
%     MatchTable.natImSig = natImSig(:);
%
%     % ---
%     % Third fingerprint with scaled response
%     corrScaledResp_big = nan(sum(nClu),sum(nClu));
%     for ss1 = 1:nRec
%         if ~isempty(spikeData_cv{1,ss1})
%             mat1 = nanmean(spikeData_cv{1,ss1}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),4);
%             timecourse1 = zscore(nanmean(spikeData_cv{1,ss1}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),[1 4]));
%             scaledResp1 = squeeze(sum(mat1.*timecourse1,2));
%             for ss2 = 1:nRec
%                 if ~isempty(spikeData_cv{2,ss2})
%                     mat2 = nanmean(spikeData_cv{2,ss2}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),4);
%                     timecourse2 = zscore(nanmean(spikeData_cv{2,ss2}(:,bins<proc.window(2)+0.2 & bins>-0.2,:,:),[1 4]));
%                     scaledResp2 = squeeze(sum(mat2.*timecourse2,2));
%                     corrScaledResp_big(sum(nClu(1:ss1-1))+1:sum(nClu(1:ss1)), sum(nClu(1:ss2-1))+1:sum(nClu(1:ss2))) = corr(scaledResp1,scaledResp2);
%                 end
%             end
%         end
%     end
%     corrScaledResp_big = tanh(.5*atanh(corrScaledResp_big) + .5*atanh(corrTimecourse_big));
%     corrScaledResp_big = tanh(.5*atanh(corrScaledResp_big) + .5*atanh(corrScaledResp_big')); % not sure that's needed?
% 
%     % Get rank
%     [natImScaledRespRank, natImScaledRespSig] = getRank(corrScaledResp_big, SessionSwitch);
% 
%     % Save in table
%     MatchTable.natImScaledRespCorr = corrScaledResp_big(:);
%     MatchTable.natImScaledRespRank = natImScaledRespRank(:);
%     MatchTable.natImScaledRespSig = natImScaledRespSig(:);
end
%% Write to table

save(fullfile(SaveDir, 'UnitMatch.mat'), 'MatchTable', 'UMparam', 'UniqueIDConversion', '-append');
if exist('AllSessionCorrelations','var')
    save(fullfile(SaveDir, 'UnitMatch.mat'), 'AllSessionCorrelations', '-append');
end

%% Extract groups
if UMparam.RunPyKSChronicStitched
    ntimes = 2;
else
    ntimes = 1;
end
for id = 1:ntimes
    if id == 1 %Normal situation
        WithinIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
        MatchIdx = find((MatchTable.UID1 == MatchTable.UID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
        NonMatchIdx = find((MatchTable.UID1 ~= MatchTable.UID2)); % Not the same unit
        addname = [];
    elseif id == 2 % USE KS labels instead
        WithinIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 == MatchTable.RecSes2)); %Within session, same unit (cross-validation)
        MatchIdx = find((MatchTable.ID1 == MatchTable.ID2) & (MatchTable.RecSes1 ~= MatchTable.RecSes2)); %Across session, same unit (cross-validation)
        NonMatchIdx = find((MatchTable.ID1 ~= MatchTable.ID2)); % Not the same unit
        addname = 'KSlabels';
    end
    MatchIdx(isnan(MatchTable.refPopCorr(MatchIdx))) = [];
    if isempty(MatchIdx)
        disp('No Matches found... return')
        return
    end

    %% Cross-correlation ROC?
    refPopCorr = reshape(MatchTable.refPopCorr, nclus, nclus);

    if saveFig

        saveas(VenFig, fullfile(SaveDir, ['VennFunctionalScores.fig']))
        saveas(VenFig, fullfile(SaveDir, ['VennFunctionalScores.bmp']))

        figure('name', ['Functional score separatability ', addname])

        subplot(4, 3, 1)
        imagesc(refPopCorr)
        hold on
        colormap(flipud(gray))
        makepretty
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        title('Cross-correlation Fingerprint')
        axis square
        freezeColors

        subplot(4, 3, 2)
        bins = min(refPopCorr(:)):0.1:max(refPopCorr(:));
        Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
        hw = histcounts(refPopCorr(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(refPopCorr(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(refPopCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
        plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
        hold on
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0, 0, 0])
        xlabel('Cross-correlation Fingerprint')
        ylabel('Proportion|Group')
        legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
        axis square
        makepretty
        subplot(4, 3, 3)
        if any(MatchIdx)
            labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
            scores = [refPopCorr(MatchIdx)', refPopCorr(NonMatchIdx)'];
            [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
            h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
            hold all
            labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
            scores = [refPopCorr(MatchIdx)', refPopCorr(WithinIdx)'];
            [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
            h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
        end

        labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [refPopCorr(WithinIdx)', refPopCorr(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);
        axis square

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
        title(sprintf('Cross-Correlation Fingerprint AUC: %.3f, %.3f, %.3f', AUC1, AUC2, AUC3))
        makepretty
        drawnow %Something to look at while ACG calculations are ongoing

    
    end


    %% Plot ACG
    ACGCor = reshape(MatchTable.ACGCorr, nclus, nclus);
    if saveFig
        subplot(4, 3, 4)
        imagesc(ACGCor)
        hold on
        colormap(flipud(gray))
        makepretty
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        title('Autocorrelogram Correlation')
        axis square

        freezeColors

        subplot(4, 3, 5)
        % Subtr = repmat(diag(ACGCor),1,size(ACGCor,1));
        % ACGCor = ACGCor - Subtr; % Subtract diagonalcorrelations
        bins = min(ACGCor(:)):0.1:max(ACGCor(:));
        Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
        hw = histcounts(ACGCor(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(ACGCor(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(ACGCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
        plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
        hold on
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0, 0, 0])
        xlabel('Autocorrelogram Correlation')
        ylabel('Proportion|Group')
        %     legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
        axis square

        makepretty

        subplot(4, 3, 6)
        if any(MatchIdx)
            labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
            scores = [ACGCor(MatchIdx)', ACGCor(NonMatchIdx)'];
            [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
            h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
            hold all
            labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
            scores = [ACGCor(MatchIdx)', ACGCor(WithinIdx)'];
            [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
            h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
        end
        labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [ACGCor(WithinIdx)', ACGCor(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        %     legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
        title(sprintf('Autocorrelogram AUC: %.3f, %.3f, %.3f', AUC1, AUC2, AUC3))
        makepretty
        axis square

        freezeColors
    end
    %% Plot FR
    FRDiff = reshape(MatchTable.FRDiff, nclus, nclus);
    Subtr = repmat(diag(FRDiff), 1, size(FRDiff, 1));
    FRDiff = FRDiff - Subtr; % Subtract diagonalcorrelations
    if saveFig
        subplot(4, 3, 7)
        imagesc(FRDiff)
        hold on
        colormap(gray)
        makepretty
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        title('Firing rate differences')
        axis square



        subplot(4, 3, 8)
        bins = min(FRDiff(:)):0.1:5;
        Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
        hw = histcounts(FRDiff(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(FRDiff(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(FRDiff(NonMatchIdx), bins) ./ length(NonMatchIdx);
        plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
        hold on
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0, 0, 0])
        xlabel('Firing rate differences')
        ylabel('Proportion|Group')
        %     legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
        makepretty
        axis square

        subplot(4, 3, 9)
        if any(MatchIdx(:))
            labels = [zeros(1, numel(MatchIdx)), ones(1, numel(NonMatchIdx))];
            scores = [FRDiff(MatchIdx)', FRDiff(NonMatchIdx)'];
            [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
            h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
        end
        hold all
        labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
        scores = [FRDiff(MatchIdx)', FRDiff(WithinIdx)'];
        [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
        h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
        labels = [zeros(1, numel(WithinIdx)), ones(1, numel(NonMatchIdx))];
        scores = [FRDiff(WithinIdx)', FRDiff(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);
        axis square

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        %     legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
        title(sprintf('Firing rate differences AUC: %.3f, %.3f, %.3f', AUC1, AUC2, AUC3))
        makepretty
    end

    %% Plot Natural images

    natImRespCorr = reshape(MatchTable.natImRespCorr, nclus, nclus);
    if saveFig
        subplot(4, 3, 10)
        imagesc(natImRespCorr)
        hold on
        colormap(flipud(gray))
        makepretty
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
        title('natIm Fingerprint')
        axis square
        freezeColors


        if ~all(isnan(natImRespCorr(MatchIdx)))
            subplot(4, 3, 11)
            bins = min(natImRespCorr(:)):0.1:max(natImRespCorr(:));
            Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
            hw = histcounts(natImRespCorr(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(natImRespCorr(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(natImRespCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
            plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
            hold on
            plot(Vector, hm, 'color', [0, 0.5, 0])
            plot(Vector, hn, 'color', [0, 0, 0])
            xlabel('natIm Fingerprint')
            ylabel('Proportion|Group')
            %         legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
            axis square
            makepretty

            subplot(4, 3, 11)
            bins = min(natImRespCorr(:)):0.1:max(natImRespCorr(:));
            Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
            hw = histcounts(natImRespCorr(WithinIdx), bins) ./ length(WithinIdx);
            hm = histcounts(natImRespCorr(MatchIdx), bins) ./ length(MatchIdx);
            hn = histcounts(natImRespCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
            plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
            hold on
            plot(Vector, hm, 'color', [0, 0.5, 0])
            plot(Vector, hn, 'color', [0, 0, 0])
            xlabel('natIm Fingerprint')
            ylabel('Proportion|Group')
            %         legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
            axis square
            makepretty
            subplot(4, 3, 12)
            if any(MatchIdx) && ~all(isnan(natImRespCorr(MatchIdx)))
                labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
                scores = [natImRespCorr(MatchIdx)', natImRespCorr(NonMatchIdx)'];

                [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
                h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
                hold all
                labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];

                scores = [natImRespCorr(MatchIdx)', natImRespCorr(WithinIdx)'];
                [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
                h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
            end

            labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
            scores = [natImRespCorr(WithinIdx)', natImRespCorr(NonMatchIdx)'];
            [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
            h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);
            axis square

            plot([0, 1], [0, 1], 'k--')
            xlabel('False positive rate')
            ylabel('True positive rate')
            %         legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
            title(sprintf('natIm Fingerprint AUC: %.3f, %.3f, %.3f', AUC1, AUC2, AUC3))
            makepretty
            drawnow %Something to look at while ACG calculations are ongoing
        end
    end

    %% save
    set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
    if saveFig
        try
            saveas(gcf, fullfile(SaveDir, [addname, 'FunctionalScoreSeparability.fig']))
        catch 
        end
        try
            saveas(gcf, fullfile(SaveDir, [addname, 'FunctionalScoreSeparability.png']))
        catch 
        end
    end
end