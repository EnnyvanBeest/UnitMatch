function EvaluatingUnitMatch(DirToEvaluate)
%% Evaluating UnitMatch
% DirToEvaluate = 'H:\MatchingUnits\Output\AL032\UnitMatch';
stepsize = 0.01;
% Load UnitMatch Output
Output2Evaluate = dir(fullfile(DirToEvaluate,'UnitMatch','UnitMatch.mat'));
Output2Evaluate = matfile(fullfile(Output2Evaluate.folder,Output2Evaluate.name));

Model2Evaluate = dir(fullfile(DirToEvaluate,'UnitMatch','UnitMatchModel.mat'));
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
VariableNames = BestMdl.VariableNames;
nRows = ceil(sqrt(length(VariableNames)+1));
nCols = round(sqrt(length(VariableNames)+1));
%% 'Ground truth' (or as best as we can): Take the set where ID1 == ID2
GTidx = find(MatchTable.ID1 == MatchTable.ID2);

% Were these neurons assigned as final matches?
FoundAsMatch = sum(MatchTable.UID1(GTidx)==MatchTable.UID2(GTidx))./length(GTidx);
disp(['UnitMatch identified ' num2str(round(FoundAsMatch*1000)./10) '% of Kilosort single units as such'])

% What is the match probability of 
MatchIdx = GTidx(find(MatchTable.UID1(GTidx)==MatchTable.UID2(GTidx)));
NonMatchIDx = GTidx(find(MatchTable.UID1(GTidx)~=MatchTable.UID2(GTidx)));
Edges = [0:stepsize:1];
Vect = stepsize/2:stepsize:1-stepsize/2;
hM = histcounts(MatchTable.MatchProb(MatchIdx),Edges)./length(MatchIdx);
hNM = histcounts(MatchTable.MatchProb(NonMatchIDx),Edges)./length(NonMatchIDx);

figure; 
subplot(nRows,nCols,1)
plot(Vect,hNM,'b-')
hold on
plot(Vect,hM,'r-')
legend('KSMatch, Non UnitMatch','KSMatch and UnitMatch')
title('Match probability of KS identified clusters')
ylabel('nSamples')
makepretty

%% Now show distributions for the other parameters
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
linkaxes
%% Leave-One-Out Analysis; would we have performed better when we left out one of the variables?
figure('name','Leave Out Analysis')
FoundAsMatch = nan(length(VariableNames)+1,1);
for vid = 1:length(VariableNames)+1
    VarNameSet = VariableNames;
    if vid>1
        VarNameSet(vid-1) = []; % Leave this one out
    end
    % Re-run model
    Tbl = MatchTable(:,ismember(MatchTable.Properties.VariableNames,VarNameSet));
    [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(VariableNames,VarNameSet),:),[0,1],BestMdl.Priors);

    % Were these neurons assigned as final matches?
    FoundAsMatch(vid) = sum(Fakelabel(GTidx))./length(GTidx);
    if vid==1
        disp(['Full model UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
    else
        disp(['Leave out ' VariableNames{vid-1} ' UnitMatch identified ' num2str(round(FoundAsMatch(vid)*1000)./10) '% of Kilosort single units as such'])
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
        title(['wo ' VariableNames{vid-1} ': ' num2str(round(FoundAsMatch(vid)*1000)./10) '%'])
    end
    ylabel('nSamples')
    xlabel('p|Match')
    makepretty

end
linkaxes
%% Leave-One-In Analysis; would we have performed better when we left out one of the variables?
figure('name','Leave In Analysis')
FoundAsMatch = nan(length(VariableNames)+1,1);
for vid = 1:length(VariableNames)+1
    if vid>1
        VarNameSet = VariableNames{vid-1};
    else
        VarNameSet = VariableNames;
    end
    % Re-run model
    Tbl = MatchTable(:,ismember(MatchTable.Properties.VariableNames,VarNameSet));
    [Fakelabel, LeaveOutProb, performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(VariableNames,VarNameSet),:),[0,1],BestMdl.Priors);

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

return
