function EvaluatingUnitMatch(SaveDir)

%% Evaluating UnitMatch
% DirToEvaluate = 'H:\MatchingUnits\Output\AL032\UnitMatch';
stepsize = 0.01;
% Load UnitMatch Output
Output2Evaluate = dir(fullfile(SaveDir, 'UnitMatch.mat'));
Output2Evaluate = matfile(fullfile(Output2Evaluate.folder, Output2Evaluate.name));
UMparam = Output2Evaluate.UMparam;

ShowScores = 0;
for id = 1 %:2 % Loop: first use model that was used, then see if standard model would work better
    if id == 1
        Model2Evaluate = dir(fullfile(SaveDir, 'UnitMatchModel.mat'));
        disp('Evaluating the model as it was ran by user')
        extraname = 'Used Model';
    else
        disp('Now assuming the standard model provided by UnitMatch')
        % Compare to standard model
        P = mfilename('fullpath');
        P = fileparts(P);
        Model2Evaluate = dir(fullfile(P, 'UnitMatchModel.mat'));
        extraname = 'Standard Model';
    end
    Model2Evaluate = matfile(fullfile(Model2Evaluate.folder, Model2Evaluate.name));

    % Extract MatchTable
    MatchTable = Output2Evaluate.MatchTable;
    UniqueIDConversion = Output2Evaluate.UniqueIDConversion;
    nclus = length(UniqueIDConversion.Path4UnitNPY);
    ndays = length(unique(UniqueIDConversion.recsesAll));
    % Extract Model
    BestMdl = Model2Evaluate.BestMdl;
    if ~isfield(BestMdl, 'Priors')
        priorMatch = 1 - ((nclus + nclus .* sqrt(ndays-1) * UMparam.ExpectMatches) ./ (nclus * nclus)); %Punish multiple days (unlikely to find as many matches after a few days)
        BestMdl.Priors = [priorMatch, 1 - priorMatch];
    end
    if id == 1
        VariableNames = BestMdl.VariableNames;
    end
    VariableIndex = find(ismember(BestMdl.VariableNames, VariableNames));
    nRows = ceil(sqrt(length(VariableNames)+1));
    nCols = round(sqrt(length(VariableNames)+1));

    %% Assuming UnitMatch worked, how many units were tracked?
    if id == 1
        if UMparam.GoodUnitsOnly
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
        end
        TrackingPerformance = nan(3, 0); % Difference between recording number, % Tracked units %maximum possibility
        MatchProb = reshape(MatchTable.MatchProb, nclus, nclus);
        for did1 = 1:ndays
            for did2 = 1:ndays
                if did2 <= did1
                    continue
                end
                thesedaysidx = find(ismember(recses, [did1, did2]));
                % can possibly only track these many units:
                nMax = min([sum(recses(thesedaysidx) == did1), sum(recses(thesedaysidx) == did2)]);
                % only matches across days
                recsesthesedays = recses(thesedaysidx);
                nMatches = sum(arrayfun(@(X) length(unique(recsesthesedays(UniqueID(thesedaysidx) == X))) > 1, unique(UniqueID(thesedaysidx))));
                TrackingPerformance = cat(2, TrackingPerformance, [did2 - did1, nMatches, nMax]');
            end
        end

        figure('name', ['TrackingPerformance ', UMparam.SaveDir])
        subplot(1, 2, 1)
        scatter(TrackingPerformance(3, :), TrackingPerformance(2, :), 20, [0, 0, 0], 'filled')
        hold on
        xlims = get(gca, 'xlim');
        line([0, max(xlims)], [0, max(xlims)], 'color', [0, 0, 0])
        xlabel('Number Units Available')
        ylabel('Number Units Tracked')

        subplot(1, 2, 2)
        scatter(TrackingPerformance(1, :), TrackingPerformance(2, :)./TrackingPerformance(3, :), 20, [0, 0, 0], 'filled')
        xlabel('\Delta Recordings')
        ylabel('Proportion tracked cells')
    end

    %% 'Ground truth' (or as best as we can): Take the set where ID1 == ID2 (False Negatives)
    if UMparam.RunPyKSChronicStitched
        disp('Across and within recording cross-validation: PyKS was ran stitched')
        GTMatchidx = find(MatchTable.ID1 == MatchTable.ID2); % & MatchTable.RecSes1 ~= MatchTable.RecSes2);
        GTNonMatch = find(MatchTable.ID1 ~= MatchTable.ID2);
    else
        disp('Within recording cross-validation: PyKS was ran on individual sessions')
        GTMatchidx = find(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);
        GTNonMatch = find(MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);
    end

    % What is the match probability of
    Edges = [0:stepsize:1];
    Vect = stepsize / 2:stepsize:1 - stepsize / 2;
    hM = histcounts(MatchTable.MatchProb(GTMatchidx), Edges) ./ length(GTMatchidx);
    hNM = histcounts(MatchTable.MatchProb(GTNonMatch), Edges) ./ length(GTNonMatch);

    if ShowScores
        figure('name', ['Scores ', extraname, ' False Negatives']);
        subplot(nRows, nCols, 1)
        plot(Vect, hNM, 'b-')
        hold on
        plot(Vect, hM, 'r-')
        legend('KS not same unit', 'KS Same unit')
        title('Match probability of KS identified clusters')
        ylabel('nSamples')
        makepretty

        % Now show distributions for the other parameters
        for vid = 1:length(VariableNames)
            eval(['hM = histcounts(MatchTable.', VariableNames{vid}, '(GTMatchidx),Edges)./length(GTMatchidx);'])
            eval(['hNM = histcounts(MatchTable.', VariableNames{vid}, '(GTNonMatch),Edges)./length(GTNonMatch);'])
            subplot(nRows, nCols, 1+vid)
            plot(Vect, hNM, 'b-')
            hold on
            plot(Vect, hM, 'r-')
            title(['Score Distribution ', VariableNames{vid}])
            ylabel('nSamples')
            makepretty
        end
    end
    % Leave-One-Out Analysis; would we have performed better when we left out one of the variables?
    figure('name', ['Leave Out Analysis ', extraname])
    if id == 1
        FN = nan(length(VariableNames)+1, 2); %False negatives
        FP = FN; % False positives
        FalseNegativeChanges = nan(2, length(VariableNames)+1);
        FalsePositiveChanges = nan(2, length(VariableNames)+1);
    end

    for vid = 1:length(VariableNames) + 1
        VarNameSet = VariableNames;
        if vid > 1
            VarNameSet(vid-1) = []; % Leave this one out
        end
        % Re-run model
        Tbl = MatchTable(:, ismember(MatchTable.Properties.VariableNames, VarNameSet));
        [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl, BestMdl.Parameterkernels(:, ismember(BestMdl.VariableNames, VarNameSet), :), [0, 1], BestMdl.Priors);

        % Were these neurons assigned as final matches?
        FN(vid, id) = 1 - (sum(Fakelabel(GTMatchidx)) ./ length(GTMatchidx));
        if vid == 1
            FalseNegativeChanges(id, vid) = FN(vid, id) - FN(1, 1); % Negative is better
        else
            FalseNegativeChanges(id, vid) = FN(vid, id) - FN(1, id);
        end
        FP(vid, id) = sum(Fakelabel(GTNonMatch)) ./ length(GTNonMatch);
        if vid == 1
            FalsePositiveChanges(id, vid) = FP(vid, id) - FP(1, 1);
        else
            FalsePositiveChanges(id, vid) = FP(vid, id) - FP(1, id);
        end
        %         if vid==1
        %             disp(['Full model UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         else
        %             disp(['Leave out ' VariableNames{vid-1} ' UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         end

        % What is the match probability of
        hM = histcounts(LeaveOutProb(GTMatchidx, 2), Edges) ./ length(GTMatchidx);
        hNM = histcounts(LeaveOutProb(GTNonMatch, 2), Edges) ./ length(GTNonMatch);

        subplot(nRows, nCols, vid)
        plot(Vect, hNM, 'b-')
        hold on
        plot(Vect, hM, 'r-')
        if vid == 1
            title(['Full model: FN=', num2str(round(FN(vid, id)*1000)./10), '%, FP=', num2str(round(FP(vid, id)*1000)./10), '%'])
            legend('Not same (KS)', 'Same (KS)')
        else
            title(['wo ', VariableNames{vid-1}, ': FN=', num2str(round(FN(vid, id)*1000)./10), '%, PF=', num2str(round(FP(vid, id)*1000)./10), '%'])
        end
        ylabel('nSamples')
        xlabel('p|Match')
        makepretty

    end
    linkaxes

end

%% Results:
disp(['False negative rate: ', num2str(round(FN(1, 1)*1000)/10), '%'])
disp(['False positive rate: ', num2str(round(FP(1, 1)*1000)/10), '%'])

%% Some Advice:
if FalseNegativeChanges(2, 1) < FalseNegativeChanges(1, 1) & FalsePositiveChanges(2, 1) < FalsePositiveChanges(1, 1)
    disp('We advise to use the standard model on this data:')
    disp('PrepareClusInfoparams.ApplyExistingBayesModel = 1');
    id2take = 2;
else
    disp('We advise to use the model you just used on this data:')
    id2take = 1;
end
if any((FalseNegativeChanges(id2take, 2:end) < 0 & FalsePositiveChanges(id2take, 2:end) < 0))
    disp(['To optimize detection, consider taking out:'])
    VariableNames((FalseNegativeChanges(id2take, 2:end) > 0 & FalsePositiveChanges(id2take, 2:end) < 0))

    disp(['So please try running the model with: PrepareClusInfoparams.Scores2Include = '])
    VariableNames(~(FalseNegativeChanges(id2take, 2:end) > 0 & FalsePositiveChanges(id2take, 2:end) < 0))
else
    disp('This parameter set is probably the best')
end

%% AUC
if exist(fullfile(SaveDir, 'AUC.mat'))
    load(fullfile(SaveDir, 'AUC.mat'));
    AUC = AUCStruct.AUC;
    ParamNames = AUCStruct.ParamNames(~isnan(AUC));
    AUC = AUC(~isnan(AUC));

    figure('name', 'AUC versus FN Change')
    subplot(1, 2, 1)
    scatter(AUC(cellfun(@(X) find(ismember(ParamNames, X)), VariableNames)), FN(2:end, 1), 25, [0, 0, 0], 'filled');
    xlabel('AUC')
    ylabel('Leave-out FN (%)')
    title('False Negatives')
    makepretty

    subplot(1, 2, 2)
    scatter(AUC(cellfun(@(X) find(ismember(ParamNames, X)), VariableNames)), FP(2:end, 1), 25, [0, 0, 0], 'filled');
    xlabel('AUC')
    ylabel('Leave-out FP (%)')
    title('False Positives')
    makepretty

    idx = AUC(~ismember(ParamNames, VariableNames)) > min(AUC(ismember(ParamNames, VariableNames)));
    [~, sortid] = sort(AUC, 'descend');
    disp('AUC order:')
    ParamNames(sortid)
    AUC(sortid)
end
return
