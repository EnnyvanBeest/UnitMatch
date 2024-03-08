function ComputeFunctionalScores(SaveDir, loadMATsToSave)

if nargin < 2 || isempty(loadMATsToSave)
    TmpFile = matfile(fullfile(SaveDir, 'UnitMatch.mat')); % Access saved file
else
    TmpFile = load(fullfile(SaveDir, 'UnitMatch.mat'));
end

load(fullfile(SaveDir, 'UnitMatch.mat'), 'MatchTable', 'UMparam', 'UniqueIDConversion');
UMparam.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows

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

if length(UMparam.AllRawPaths{1}) > 1 %Reshape for Stitched
    UMparam.AllRawPaths = arrayfun(@(X) UMparam.AllRawPaths{1}(X), 1:length(UMparam.AllRawPaths{1}), 'uni', 0);
end

% Load SP
disp('Loading spike information...')
nKSFiles = length(AllKSDir);
sp = cell(1, nKSFiles);
for did = 1:nKSFiles
    tmp = matfile(fullfile(AllKSDir{did}, 'PreparedData.mat'));
    sptmp = tmp.sp;
    clear tmp

    % Only keep parameters used
    sp{did}.st = sptmp.st;
    sp{did}.spikeTemplates = sptmp.spikeTemplates;
    sp{did}.spikeAmps = sptmp.spikeAmps;
    % Replace recsesid with subsesid
    if UMparam.RunPyKSChronicStitched
        sp{did}.RecSes = sptmp.RecSes;
    else
        sp{did}.RecSes = repmat(did, size(sp{did}.st));
    end
end
clear sptmp


% Add all spikedata in one spikes struct - can be used for further analysis
sp = [sp{:}];
spnew = struct;
fields = fieldnames(sp(1));
for fieldid = 1:length(fields)
    try
        eval(['spnew.', fields{fieldid}, '= cat(1,sp(:).', fields{fieldid}, ');'])
    catch ME
        if strcmp(ME.message, 'Out of memory.')
            eval(['spnew.', fields{fieldid}, ' = sp(1).', fields{fieldid}, ';'])
            for tmpid = 2:length(sp)
                eval(['spnew.', fields{fieldid}, ' = cat(1,spnew.', fields{fieldid}, ', sp(tmpid).', fields{fieldid}, ');'])
            end
        else
            eval(['spnew.', fields{fieldid}, '= cat(2,sp(:).', fields{fieldid}, ');'])
        end
    end
end
sp = spnew;
clear spnew

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

if ~any(ismember(MatchTable.Properties.VariableNames, 'FingerprintCor')) % If it already exists in table, skip this entire thing


    % for paper a figure:
    %     figure
    %     MatchProbability = reshape(MatchTable.MatchProb, nclus, nclus);
    %     [r, c] = find(MatchProbability > UMparam.ProbabilityThreshold); %Find matches
    %     Pairs = cat(2, r, c);
    %     Pairs = sortrows(Pairs);
    %     Pairs = unique(Pairs, 'rows');
    %
    %     pairidx = recses(Pairs(:,1)) == 1 & recsesGood(Pairs(:,2))==2;
    %     PairsTmp = Pairs(pairidx,:);
    %     % Only use every 'unit' once --> take the highest scoring matches
    %     [~,id1,~]=unique(PairsTmp(:,1),'stable');
    %     PairsTmp = PairsTmp(id1,:);
    %     [~,id1,~]=unique(PairsTmp(:,2),'stable');
    %     PairsTmp = PairsTmp(id1,:);

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

        %         Unit2TakeIdx = PairsTmp(1:60,rid); % Only take each unit once

        %         subplot(1,nRec,rid)
        %         imagesc(sr(Unit2TakeIdx-SessionSwitch(rid)+1,1:10000),[0 1]); hold on
        %         colormap(flipud(gray))
        %         xlabel('time')
        %         ylabel('Reference units')
        %         makepretty

        % Define folds (two halves)
        idx_fold1 = 1:floor(size(sr, 2)./2);
        idx_fold2 = floor(size(sr, 2)./2) + 1:floor(size(sr, 2)./2) * 2;

        % Find cross-correlation in first and second half of session
        try
            SessionCorrelations{rid}.fold1 = corr(sr(:, idx_fold1)', sr(:, idx_fold1)')';
            SessionCorrelations{rid}.fold2 = corr(sr(:, idx_fold2)', sr(:, idx_fold2)')';
        catch
            keyboard
        end

        % Nan the diagonal
        SessionCorrelations{rid}.fold1(logical(eye(size(SessionCorrelations{rid}.fold1)))) = nan;
        SessionCorrelations{rid}.fold2(logical(eye(size(SessionCorrelations{rid}.fold2)))) = nan;
    end

    %% Get correlation matrics for fingerprint correlations
    sessionCorrelationsAll = cell(1, nRec);
    for did = 1:nRec
        % Load sp for correct day
        if length(AllKSDir) > 1
            tmp = matfile(fullfile(AllKSDir{did}, 'PreparedData.mat'));
        else %Stitched
            tmp = matfile(fullfile(AllKSDir{1}, 'PreparedData.mat'));
        end
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
    [FingerprintR, RankScoreAll, SigMask, AllSessionCorrelations] = CrossCorrelationFingerPrint(sessionCorrelationsAll, Pairs, OriID, recses, drawdrosscorr);

    % Save in table
    MatchTable.FingerprintCor = FingerprintR(:);
    MatchTable.RankScore = RankScoreAll(:);
    MatchTable.SigFingerprintR = SigMask(:);

    if nargin < 2 || isempty(loadMATsToSave)
        TmpFile.Properties.Writable = true;
    end
    TmpFile.MatchTable = MatchTable; % Overwrite
    TmpFile.AllSessionCorrelations = AllSessionCorrelations;

    if nargin < 2 || isempty(loadMATsToSave)

        TmpFile.Properties.Writable = false;
    else

        movefile(fullfile(SaveDir, 'UnitMatch.mat'), fullfile(SaveDir, 'UnitMatch_prev.mat'))
        save(fullfile(SaveDir, 'UnitMatch.mat'), '-struct', 'TmpFile');
        delete(fullfile(SaveDir, 'UnitMatch_prev.mat'))
    end

    %% Compare to functional scores

    figure;

    subplot(1, 3, 1)
    imagesc(RankScoreAll(SortingOrder, SortingOrder) == 1 & SigMask(SortingOrder, SortingOrder) == 1)
    hold on
    arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    colormap(flipud(gray))
    title('Rankscore == 1*')
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
    imagesc(MatchProbability(SortingOrder, SortingOrder) >= UMparam.ProbabilityThreshold | (MatchProbability(SortingOrder, SortingOrder) > 0.05 & RankScoreAll(SortingOrder, SortingOrder) == 1 & SigMask(SortingOrder, SortingOrder) == 1));
    % imagesc(MatchProbability>=0.99 | (MatchProbability>=0.05 & RankScoreAll==1 & SigMask==1))
    hold on
    arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    colormap(flipud(gray))
    title('Matching probability + rank')
    makepretty
    saveas(gcf, fullfile(SaveDir, 'RankScoreVSProbability.fig'))
    saveas(gcf, fullfile(SaveDir, 'RankScoreVSProbability.bmp'))

    tmpf = triu(FingerprintR);
    tmpm = triu(MatchProbability);
    tmpr = triu(RankScoreAll);
    tmpr = tmpr(tmpf ~= 0);
    tmpm = tmpm(tmpf ~= 0);
    tmpf = tmpf(tmpf ~= 0);
    figure;
    scatter(tmpm, tmpf, 14, tmpr, 'filled')
    colormap(cat(1, [0, 0, 0], winter))
    xlabel('Match Probability')
    ylabel('Cross-correlation fingerprint')
    makepretty
    saveas(gcf, fullfile(SaveDir, 'RankScoreVSProbabilityScatter.fig'))
    saveas(gcf, fullfile(SaveDir, 'RankScoreVSProbabilityScatter.bmp'))
end

%% Get ACG fingerprints correlations

if ~any(ismember(MatchTable.Properties.VariableNames, 'ACGCorr')) % If it already exists in table, skip this entire thing

    %% Compute ACG and correlate them between units
    % This is very time consuming
    disp('Computing ACG, this will take some time...')
    tvec = -UMparam.ACGduration / 2:UMparam.ACGbinSize:UMparam.ACGduration / 2;
    ACGMat = nan(length(tvec), 2, nclus);
    FR = nan(2, nclus);
    for clusid = 1:nclus %parfot QQ
        for cv = 1:2
            idx1 = find(sp.spikeTemplates == OriID(clusid) & sp.RecSes == recses(clusid));
            if ~isempty(idx1) && length(idx1) > UMparam.sampleamount
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
    MatchTable.ACGCorr = ACGCorr(:);

    ACGCorr(isnan(ACGCorr)) = nanmin(ACGCorr(:)); % should not participate to the rank
    % Find rank
    for did1 = 1:nRec
        for did2 = 1:nRec
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1) - 1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1) - 1;
            [~, idx] = sort(ACGCorr(clusIdxD1All, clusIdxD2All), 2, 'descend');
            for c = 1:numel(clusIdxD1All)
                ACGRankScore(clusIdxD1All(c), clusIdxD2All(idx(c, :))) = 1:numel(clusIdxD2All);
            end
        end
    end
    MatchTable.ACGRankScore = ACGRankScore(:);
end

%% Get FR difference

if ~any(ismember(MatchTable.Properties.VariableNames, 'FRDiff'))
    FR = repmat(permute(FR, [2, 1]), [1, 1, nclus]);
    FRDiff = abs(squeeze(FR(:, 2, :)-permute(FR(:, 1, :), [3, 2, 1])));
    MatchTable.FRDiff = FRDiff(:);

    FRDiff(isnan(FRDiff)) = nanmax(FRDiff(:)); % should not participate to the rank
    % Find rank
    for did1 = 1:nRec
        for did2 = 1:nRec
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1) - 1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1) - 1;
            [~, idx] = sort(FRDiff(clusIdxD1All, clusIdxD2All), 2, 'ascend');
            for c = 1:numel(clusIdxD1All)
                FRRankScore(clusIdxD1All(c), clusIdxD2All(idx(c, :))) = 1:numel(clusIdxD2All);
            end
        end
    end
    MatchTable.FRRankScore = FRRankScore(:);

end

%% Get natural images fingerprints correlations

if ~any(ismember(MatchTable.Properties.VariableNames, 'NatImCorr')) || all(isnan(MatchTable.NatImCorr)) % If it already exists in table, skip this entire thing
    % Param for processing
    proc.window = [-0.3, 0.5, ... % around onset
        0.0, 0.5]; % around offset
    proc.binSize = 0.002; % in ms
    proc.smoothSize = 5; % PSTH smoothing filter
    gw = gausswin(proc.smoothSize, 3);
    proc.smWin = gw ./ sum(gw);

    nRec = length(RecOpt); % not always the same!! numel(UMparam.AllRawPaths);
    nClu = nan(1, nRec);
    for ss = 1:nRec;
        nClu(ss) = numel(unique(MatchTable.ID1(MatchTable.RecSes1 == RecOpt(ss))));
    end
    spikeData_cv = cell(1, 2*nRec);
    for ss = 1:nRec
        % Get the original binFile (also for stitched?)
        if iscell(UMparam.AllRawPaths(RecOpt(ss)))
            fileThis = UMparam.AllRawPaths{RecOpt(ss)};
            binFileRef = fullfile([fileThis.folder, fileThis.name]);
        else
            binFileRef = fullfile(UMparam.AllRawPaths{RecOpt(ss)}.folder, UMparam.AllRawPaths{RecOpt(ss)}.name);
        end

        % Find the associated experiments
        if ispc
            exp2keep = getNatImExpRef(binFileRef);
        else
            exp2keep = [];
            % Julie Fabre: on linux the servers are mounted differently, and the above
            % function needs to be substantially changed to work. Since I
            % don't need to run it, just added an "ispc" flag.
        end

        if ~isempty(exp2keep)
            % Get the spikes
            try
                st = sp.st(sp.RecSes == RecOpt(ss));
                clu = sp.spikeTemplates(sp.RecSes == RecOpt(ss));
            catch ME
                keyboard
            end

            spikesAll.times = st;
            spikesAll.clusters = clu;
            spikesAll.clusterIDs = unique(clu); % follow the same units across days

            % Get the natim responses
            try
                spikeData = getNatImResp(spikesAll, exp2keep, binFileRef, proc);
                clusterIDs = spikesAll.clusterIDs;

                % Split in two halves and subselect units
                cluIdx = ismember(clusterIDs, unique(MatchTable.ID1(MatchTable.RecSes1 == RecOpt(ss))));
                currIdx = (ss - 1) * 2;
                spikeData_cv{currIdx+1} = spikeData(:, :, cluIdx, 1:2:end); % odd
                spikeData_cv{currIdx+2} = spikeData(:, :, cluIdx, 2:2:end); % even
            catch ME
                disp(ME)
            end
        end
    end

    % Perform CCA across recordings

    [corrMat, ~] = computeNatImCorr(spikeData_cv);

    % Reshape the matrix to a single one with correct clusters
    %%% SIMILAR TO WHAT ENNY IS DOING (e.g., 1x2 and 2x1, but not 1x1 and 2x2);
    %%% NOT OPTIMAL SINCE NOT USING ALL THE POSSIBLE PAIRS
    corrWCCA_big = nan(sum(nClu), sum(nClu));
    corrWCCA_1x2 = corrMat(1:2:end, 2:2:end);
    for ss1 = 1:nRec
        for ss2 = 1:nRec
            if ~all(isnan(corrWCCA_1x2{ss1, ss2}(:)))
                try
                    corrWCCA_big(sum(nClu(1:ss1-1))+1:sum(nClu(1:ss1)), sum(nClu(1:ss2-1))+1:sum(nClu(1:ss2))) = corrWCCA_1x2{ss1, ss2};
                catch ME
                    keyboard
                end
            end
        end
    end
    corrWCCA_big = .5 * corrWCCA_big + .5 * corrWCCA_big'; % not sure that's needed?

    MatchTable.NatImCorr = corrWCCA_big(:);


    % RANK
    corrWCCA_big_noNaN = corrWCCA_big;
    corrWCCA_big_noNaN(isnan(corrWCCA_big)) = nanmin(corrWCCA_big(:)); % should not participate to the rank
    % Find rank
    NImgRankScore = nan(size(corrWCCA_big));
    for did1 = 1:nRec
        for did2 = 1:nRec
            % Get the indices
            clusIdxD1All = SessionSwitch(did1):SessionSwitch(did1+1) - 1;
            clusIdxD2All = SessionSwitch(did2):SessionSwitch(did2+1) - 1;
            if ~all(isnan(mat2vec(corrWCCA_big(clusIdxD1All, clusIdxD2All))))
                [~, idx] = sort(corrWCCA_big_noNaN(clusIdxD1All, clusIdxD2All), 2, 'descend');
                for c = 1:numel(clusIdxD1All)
                    NImgRankScore(clusIdxD1All(c), clusIdxD2All(idx(c, :))) = 1:numel(clusIdxD2All);
                end
            end
        end
    end
    NImgRankScore(isnan(corrWCCA_big)) = nan;
    MatchTable.NImgRankScore = NImgRankScore(:);

end

%% Write to table

save(fullfile(SaveDir, 'UnitMatch.mat'), 'MatchTable', 'UMparam', 'UniqueIDConversion', '-append');
if exist('AllSessionCorrelations', 'var')
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
    MatchIdx(isnan(MatchTable.FingerprintCor(MatchIdx))) = [];
    if isempty(MatchIdx)
        disp('No Matches found... return')
        return
    end

    %% Cross-correlation ROC?
    figure('name', ['Functional score separatability ', addname])

    subplot(4, 3, 1)
    imagesc(reshape(MatchTable.FingerprintCor, nclus, nclus))
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

    FingerprintCor = reshape(MatchTable.FingerprintCor, nclus, nclus);
    % Subtr = repmat(diag(FingerprintCor),1,size(FingerprintCor,1));
    % FingerprintCor = FingerprintCor - Subtr; % Subtract diagonalcorrelations

    subplot(4, 3, 2)
    bins = min(FingerprintCor(:)):0.1:max(FingerprintCor(:));
    Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
    hw = histcounts(FingerprintCor(WithinIdx), bins) ./ length(WithinIdx);
    hm = histcounts(FingerprintCor(MatchIdx), bins) ./ length(MatchIdx);
    hn = histcounts(FingerprintCor(NonMatchIdx), bins) ./ length(NonMatchIdx);
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
        scores = [FingerprintCor(MatchIdx)', FingerprintCor(NonMatchIdx)'];
        [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
        h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
        hold all
        labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];
        scores = [FingerprintCor(MatchIdx)', FingerprintCor(WithinIdx)'];
        [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
        h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
    end

    labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
    scores = [FingerprintCor(WithinIdx)', FingerprintCor(NonMatchIdx)'];
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

    %% Plot ACG
    subplot(4, 3, 4)
    imagesc(reshape(MatchTable.ACGCorr, nclus, nclus))
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
    ACGCor = reshape(MatchTable.ACGCorr, nclus, nclus);
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

    %% Plot FR
    subplot(4, 3, 7)
    imagesc(reshape(MatchTable.FRDiff, nclus, nclus))
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

    FRDiff = reshape(MatchTable.FRDiff, nclus, nclus);
    Subtr = repmat(diag(FRDiff), 1, size(FRDiff, 1));
    FRDiff = FRDiff - Subtr; % Subtract diagonalcorrelations

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

    %% Plot Natural images

    NatImCorr = reshape(MatchTable.NatImCorr, nclus, nclus);

    subplot(4, 3, 10)
    imagesc(NatImCorr)
    hold on
    colormap(flipud(gray))
    makepretty
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X), SessionSwitch(X)], get(gca, 'ylim'), 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    arrayfun(@(X) line(get(gca, 'xlim'), [SessionSwitch(X), SessionSwitch(X)], 'color', [1, 0, 0]), 2:length(SessionSwitch), 'Uni', 0)
    title('NatIm Fingerprint')
    axis square
    freezeColors


    if ~all(isnan(NatImCorr(MatchIdx)))
        subplot(4, 3, 11)
        bins = min(NatImCorr(:)):0.1:max(NatImCorr(:));
        Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
        hw = histcounts(NatImCorr(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(NatImCorr(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(NatImCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
        plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
        hold on
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0, 0, 0])
        xlabel('NatIm Fingerprint')
        ylabel('Proportion|Group')
        %         legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
        axis square
        makepretty

        subplot(4, 3, 11)
        bins = min(NatImCorr(:)):0.1:max(NatImCorr(:));
        Vector = [bins(1) + 0.1 / 2:0.1:bins(end) - 0.1 / 2];
        hw = histcounts(NatImCorr(WithinIdx), bins) ./ length(WithinIdx);
        hm = histcounts(NatImCorr(MatchIdx), bins) ./ length(MatchIdx);
        hn = histcounts(NatImCorr(NonMatchIdx), bins) ./ length(NonMatchIdx);
        plot(Vector, hw, 'color', [0.5, 0.5, 0.5])
        hold on
        plot(Vector, hm, 'color', [0, 0.5, 0])
        plot(Vector, hn, 'color', [0, 0, 0])
        xlabel('NatIm Fingerprint')
        ylabel('Proportion|Group')
        %         legend('i=j; within recording', 'matches', 'non-matches', 'Location', 'best')
        axis square
        makepretty
        subplot(4, 3, 12)
        if any(MatchIdx) & ~all(isnan(NatImCorr(MatchIdx)))
            labels = [ones(1, numel(MatchIdx)), zeros(1, numel(NonMatchIdx))];
            scores = [NatImCorr(MatchIdx)', NatImCorr(NonMatchIdx)'];

            [X, Y, ~, AUC1] = perfcurve(labels, scores, 1);
            h(1) = plot(X, Y, 'color', [0, 0.25, 0]);
            hold all
            labels = [zeros(1, numel(MatchIdx)), ones(1, numel(WithinIdx))];

            scores = [NatImCorr(MatchIdx)', NatImCorr(WithinIdx)'];
            [X, Y, ~, AUC2] = perfcurve(labels, scores, 1);
            h(2) = plot(X, Y, 'color', [0, 0.5, 0]);
        end

        labels = [ones(1, numel(WithinIdx)), zeros(1, numel(NonMatchIdx))];
        scores = [NatImCorr(WithinIdx)', NatImCorr(NonMatchIdx)'];
        [X, Y, ~, AUC3] = perfcurve(labels, scores, 1);
        h(3) = plot(X, Y, 'color', [0.25, 0.25, 0.25]);
        axis square

        plot([0, 1], [0, 1], 'k--')
        xlabel('False positive rate')
        ylabel('True positive rate')
        %         legend([h(:)], 'Match vs No Match', 'Match vs Within', 'Within vs No Match', 'Location', 'best')
        title(sprintf('NatIm Fingerprint AUC: %.3f, %.3f, %.3f', AUC1, AUC2, AUC3))
        makepretty
        drawnow %Something to look at while ACG calculations are ongoing
    end

    %% save
    set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1])
    try
        saveas(gcf, fullfile(SaveDir, [addname, 'FunctionalScoreSeparability.fig']))
    catch ME
    end
    try
        saveas(gcf, fullfile(SaveDir, [addname, 'FunctionalScoreSeparability.png']))
    catch ME
    end
end
return