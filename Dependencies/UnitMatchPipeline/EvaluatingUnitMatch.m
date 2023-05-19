function EvaluatingUnitMatch(DirToEvaluate)
%% Evaluating UnitMatch
% DirToEvaluate = 'H:\MatchingUnits\Output\AL032\UnitMatch';
stepsize = 0.01;
% Load UnitMatch Output
Output2Evaluate = dir(fullfile(DirToEvaluate,'UnitMatch.mat'));
Output2Evaluate = matfile(fullfile(Output2Evaluate.folder,Output2Evaluate.name));
ShowScores = 0;
for id = 1:2 % Loop: first use model that was used, then see if standard model would work better
    if id==1
        Model2Evaluate = dir(fullfile(DirToEvaluate,'UnitMatchModel.mat'));
        disp('Evaluating the model as it was ran by user')
        extraname = 'Used Model';
    else
        disp('Now assuming the standard model provided by UnitMatch')
        % Compare to standard model
        P = mfilename('fullpath');
        P = fileparts(P);
        Model2Evaluate = dir(fullfile(P,'UnitMatchModel.mat'));
        extraname = 'Standard Model';
    end
    Model2Evaluate = matfile(fullfile(Model2Evaluate.folder,Model2Evaluate.name));

    % Extract MatchTable
    MatchTable = Output2Evaluate.MatchTable;
    UniqueIDConversion = Output2Evaluate.UniqueIDConversion;
    nclus = length(UniqueIDConversion.Path4UnitNPY);
    ndays = length(unique(UniqueIDConversion.recsesAll));
    % Extract Model
    BestMdl = Model2Evaluate.BestMdl;
    if ~isfield(BestMdl,'Priors')
        priorMatch = 1-(nclus*ndays)./(nclus*nclus); %Now use a slightly more lenient prior
        BestMdl.Priors = [priorMatch 1-priorMatch];
    end
    if id == 1
        VariableNames = BestMdl.VariableNames;
    end
    VariableIndex = find(ismember(BestMdl.VariableNames,VariableNames));
    nRows = ceil(sqrt(length(VariableNames)+1));
    nCols = round(sqrt(length(VariableNames)+1));

    UMparam = Output2Evaluate.UMparam;
    %% 'Ground truth' (or as best as we can): Take the set where ID1 == ID2 (False Negatives)
    if UMparam.RunPyKSChronicStitched
        disp('Across recording cross-validation: PyKS was ran stitched')
        GTidx = find(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 ~= MatchTable.RecSes2);
    else
        disp('Within recording cross-validation: PyKS was ran on individual sessions')
        GTidx = find(MatchTable.ID1 == MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);
    end

    % What is the match probability of
    MatchIdx = GTidx(find(MatchTable.UID1(GTidx)==MatchTable.UID2(GTidx)));
    NonMatchIDx = GTidx(find(MatchTable.UID1(GTidx)~=MatchTable.UID2(GTidx)));
    Edges = [0:stepsize:1];
    Vect = stepsize/2:stepsize:1-stepsize/2;
    hM = histcounts(MatchTable.MatchProb(MatchIdx),Edges)./length(MatchIdx);
    hNM = histcounts(MatchTable.MatchProb(NonMatchIDx),Edges)./length(NonMatchIDx);

    if ShowScores
        figure('name',['Scores ' extraname ' False Negatives']);
        subplot(nRows,nCols,1)
        plot(Vect,hNM,'b-')
        hold on
        plot(Vect,hM,'r-')
        legend('False Negatives','Correct')
        title('Match probability of KS identified clusters')
        ylabel('nSamples')
        makepretty

        % Now show distributions for the other parameters
        for vid = 1:length(VariableNames)
            eval(['hM = histcounts(MatchTable.'  VariableNames{vid} '(MatchIdx),Edges)./length(MatchIdx);'])
            eval(['hNM = histcounts(MatchTable.'  VariableNames{vid} '(NonMatchIDx),Edges)./length(NonMatchIDx);'])
            subplot(nRows,nCols,1+vid)
            plot(Vect,hNM,'b-')
            hold on
            plot(Vect,hM,'r-')
            title(['Score Distribution ' VariableNames{vid}])
            ylabel('nSamples')
            makepretty
        end
    end
    % Leave-One-Out Analysis; would we have performed better when we left out one of the variables?
    figure('name',['Leave Out Analysis ' extraname ' False Negatives'])
    if id==1
        FoundAsMatchN = nan(length(VariableNames)+1,2);
        FalseNegativeChanges = nan(2,length(VariableNames)+1);
    end

    for vid = 1:length(VariableNames)+1
        VarNameSet = VariableNames;
        if vid>1
            VarNameSet(vid-1) = []; % Leave this one out
        end
        % Re-run model
        Tbl = MatchTable(:,ismember(MatchTable.Properties.VariableNames,VarNameSet));
        [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(BestMdl.VariableNames,VarNameSet),:),[0,1],BestMdl.Priors);

        % Were these neurons assigned as final matches?
        FoundAsMatchN(vid,id) = sum(Fakelabel(GTidx))./length(GTidx);
        if vid==1
            FalseNegativeChanges(id,vid) = FoundAsMatchN(vid,id) - FoundAsMatchN(1,1);
        else
            FalseNegativeChanges(id,vid) = FoundAsMatchN(vid,id) - FoundAsMatchN(1,id);
        end
        %         if vid==1
        %             disp(['Full model UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         else
        %             disp(['Leave out ' VariableNames{vid-1} ' UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         end

        % What is the match probability of
        MatchIdx = GTidx(find(Fakelabel(GTidx)));
        NonMatchIDx = GTidx(find(~Fakelabel(GTidx)));
        hM = histcounts(LeaveOutProb(MatchIdx,2),Edges)./length(MatchIdx);
        hNM = histcounts(LeaveOutProb(NonMatchIDx,2),Edges)./length(NonMatchIDx);

        subplot(nRows,nCols,vid)
        plot(Vect,hNM,'b-')
        hold on
        plot(Vect,hM,'r-')
        if vid==1
            title(['Full model: ' num2str(round(FoundAsMatchN(vid,id)*1000)./10) '%'])
            legend('False Negatives','Correct')
        else
            title(['wo ' VariableNames{vid-1} ': ' num2str(round(FoundAsMatchN(vid,id)*1000)./10) '%'])
        end
        ylabel('nSamples')
        xlabel('p|Match')
        makepretty

    end
    linkaxes


    %% We should also consider the false positives
    if UMparam.RunPyKSChronicStitched
        disp('Across recording cross-validation: PyKS was ran stitched')
        GTidx = find(MatchTable.ID1 ~= MatchTable.ID2 & MatchTable.RecSes1 ~= MatchTable.RecSes2);
    else
        disp('Within recording cross-validation: PyKS was ran on individual sessions')
        GTidx = find(MatchTable.ID1~=MatchTable.ID2 & MatchTable.RecSes1 == MatchTable.RecSes2);% Now we find the Ground truth NON-matches; within day NOT the same ID
    end
    % What is the match probability of
    MatchIdx = GTidx(find(MatchTable.UID1(GTidx)==MatchTable.UID2(GTidx)));
    NonMatchIDx = GTidx(find(MatchTable.UID1(GTidx)~=MatchTable.UID2(GTidx)));
    hM = histcounts(MatchTable.MatchProb(MatchIdx),Edges)./length(MatchIdx);
    hNM = histcounts(MatchTable.MatchProb(NonMatchIDx),Edges)./length(NonMatchIDx);
    if ShowScores
        figure('name',['Scores ' extraname ' False Positives']);
        subplot(nRows,nCols,1)
        plot(Vect,hNM,'b-')
        hold on
        plot(Vect,hM,'r-')
        legend('Correct','False Positives')
        title('Match probability of KS identified clusters')
        ylabel('nSamples')
        makepretty

        % Now show distributions for the other parameters
        for vid = 1:length(VariableNames)
            eval(['hM = histcounts(MatchTable.'  VariableNames{vid} '(MatchIdx),Edges)./length(MatchIdx);'])
            eval(['hNM = histcounts(MatchTable.'  VariableNames{vid} '(NonMatchIDx),Edges)./length(NonMatchIDx);'])
            subplot(nRows,nCols,1+vid)
            plot(Vect,hNM,'b-')
            hold on
            plot(Vect,hM,'r-')
            title(['Score Distribution ' VariableNames{vid}])
            ylabel('nSamples')
            makepretty
        end
    end
    % Leave-One-Out Analysis; would we have performed better when we left out one of the variables?
    figure('name',['Leave Out Analysis ' extraname ' False Positives'])
    if id==1
        FoundAsMatch = nan(length(VariableNames)+1,2);
        FalsePositiveChanges = nan(2,length(VariableNames)+1);
    end

    for vid = 1:length(VariableNames)+1
        VarNameSet = VariableNames;
        if vid>1
            VarNameSet(vid-1) = []; % Leave this one out
        end
        % Re-run model
        Tbl = MatchTable(:,ismember(MatchTable.Properties.VariableNames,VarNameSet));
        [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(BestMdl.VariableNames,VarNameSet),:),[0,1],BestMdl.Priors);

        % Were these neurons assigned as final matches?
        FoundAsMatch(vid,id) = sum(Fakelabel(GTidx))./length(GTidx);
        if vid==1
            FalsePositiveChanges(id,vid) = FoundAsMatch(vid,id) - FoundAsMatch(1,1);
        else
            FalsePositiveChanges(id,vid) = FoundAsMatch(vid,id) - FoundAsMatch(1,id);
        end
        %         if vid==1
        %             disp(['Full model UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         else
        %             disp(['Leave out ' VariableNames{vid-1} ' UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
        %         end

        % What is the match probability of
        MatchIdx = GTidx(find(Fakelabel(GTidx)));
        NonMatchIDx = GTidx(find(~Fakelabel(GTidx)));
        hM = histcounts(LeaveOutProb(MatchIdx,2),Edges)./length(MatchIdx);
        hNM = histcounts(LeaveOutProb(NonMatchIDx,2),Edges)./length(NonMatchIDx);

        subplot(nRows,nCols,vid)
        plot(Vect,hNM,'b-')
        hold on
        plot(Vect,hM,'r-')
        if vid==1
            title(['Full model: ' num2str(round(FoundAsMatch(vid,id)*1000)./10) '%'])
            legend('Correct','False Positives')
        else
            title(['wo ' VariableNames{vid-1} ': ' num2str(round(FoundAsMatch(vid,id)*1000)./10) '%'])
        end
        ylabel('nSamples')
        xlabel('p|Match')
        makepretty

    end
    linkaxes

    if 0 % not so useful
        %% Leave-One-In Analysis; would we have performed better when we left out one of the variables?
        figure('name',['Leave In Analysis ' extraname])
        FoundAsMatch = nan(length(VariableNames)+1,1);
        for vid = 1:length(VariableNames)+1
            if vid>1
                VarNameSet = VariableNames{vid-1};
            else
                VarNameSet = VariableNames;
            end
            % Re-run model
            Tbl = MatchTable(:,ismember(MatchTable.Properties.VariableNames,VarNameSet));
            [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(BestMdl.VariableNames,VarNameSet),:),[0,1],BestMdl.Priors);

            % Were these neurons assigned as final matches?
            FoundAsMatch(vid) = sum(Fakelabel(GTidx))./length(GTidx);
            if vid==1
                disp(['Full model UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
            else
                disp(['Only ' VariableNames{vid-1} ' UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
            end

            % What is the match probability of
            MatchIdx = GTidx(find(Fakelabel(GTidx)));
            NonMatchIDx = GTidx(find(~Fakelabel(GTidx)));
            hM = histcounts(LeaveOutProb(MatchIdx,2),Edges)./length(MatchIdx);
            hNM = histcounts(LeaveOutProb(NonMatchIDx,2),Edges)./length(NonMatchIDx);

            subplot(nRows,nCols,vid)
            plot(Vect,hNM,'b-')
            hold on
            plot(Vect,hM,'r-')
            if vid==1
                title(['Full model: ' num2str(round(FoundAsMatch(vid)*1000)./10) '%'])
                legend('KSMatch, Non UnitMatch','KSMatch and UnitMatch')
            else
                title(['Only ' VariableNames{vid-1} ': ' num2str(round(FoundAsMatch(vid)*1000)./10) '%'])
            end
            ylabel('nSamples')
            xlabel('p|Match')
            makepretty

        end
        linkaxes
    end
end
%% Results:
disp(['Detection rate: '  num2str(round(FoundAsMatchN(1,1)*1000)/10)])
disp(['False positive rate: '  num2str(round(FoundAsMatch(1,1)*1000)/10)])

%% Some Advise:
if FalseNegativeChanges(2,1)>FalseNegativeChanges(1,1) & FalsePositiveChanges(2,1)<FalsePositiveChanges(1,1)
    disp('We advise to use the standard model on this data:')
    disp('PrepareClusInfoparams.ApplyExistingBayesModel = 1');
    id2take = 2;
else
    disp('We advise to use the model you just used on this data:')
    id2take = 1;
end
if any((FalseNegativeChanges(id2take,2:end)>0 & FalsePositiveChanges(id2take,2:end)<0))
    disp(['To optimize detection, consider taking out:'])
    VariableNames((FalseNegativeChanges(id2take,2:end)>0 & FalsePositiveChanges(id2take,2:end)<0))

    disp(['So please try running the model with: PrepareClusInfoparams.Scores2Include = '])
    VariableNames(~(FalseNegativeChanges(id2take,2:end)>0 & FalsePositiveChanges(id2take,2:end)<0))
else
    disp('This parameter set is probably the best')
end

return
