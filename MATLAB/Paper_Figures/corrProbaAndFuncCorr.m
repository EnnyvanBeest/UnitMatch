UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap\AL032\19011111882\2\UnitMatch\UnitMatch.mat"};
% UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap\AL036\19011116882\3\UnitMatch\UnitMatch.mat"};
% UMFile = {"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\FullAnimal_KSChanMap\CB016\Probe0\IMRO_1\UnitMatch\UnitMatch.mat"};
load(UMFile{1})


%%
MatchTable.FlipMatchProb = mat2vec(reshape(MatchTable.MatchProb,[sqrt(size(MatchTable,1)), sqrt(size(MatchTable,1))])');
tmp = MatchTable(1:sqrt(size(MatchTable,1))+1:size(MatchTable),:).MatchProb;
MatchTable.MatchProbID1 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))])');
MatchTable.MatchProbID2 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))]));
tmp = MatchTable(1:sqrt(size(MatchTable,1))+1:size(MatchTable),:).ISICorr;
MatchTable.ISICorrID1 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))])');
MatchTable.ISICorrID2 = mat2vec(repmat(tmp,[1,sqrt(size(MatchTable,1))]));

%%

FPNames = {'ISICorr','refPopCorr','natImRespCorr'};

probaBins = [0:0.05:1]; %[0 10.^(-7:0.2:0)]; %
probaBins(end) = 1.01; % to include p=1 in last bin
x = .5*probaBins(1:end-1)+.5*probaBins(2:end);

whichID = '';
% whichID = 'Liberal';
% whichID = 'Conservative';

matchIdx = MatchTable.(['UID1' whichID]) == MatchTable.(['UID2' whichID]);
% matchIdx = MatchTable.MatchProb*0.5+MatchTable.FlipMatchProb*0.5 > 0.5; 

splitUnitsUIDs = unique(MatchTable(MatchTable.MatchProb > 0.5 & MatchTable.FlipMatchProb > 0.5 ...
    & (MatchTable.ID1 ~= MatchTable.ID2) ...
    & (MatchTable.RecSes1 == MatchTable.RecSes2),:).(['UID1' whichID]));
splitUnitsIdx = ismember(MatchTable.(['UID1' whichID]),splitUnitsUIDs);

figure('Position',[809   595   700   420]);
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};
    y = nan(1, numel(probaBins)-1);
    err = nan(1, numel(probaBins)-1);
    n = nan(1, numel(probaBins)-1);
    y_match = nan(1, numel(probaBins)-1);
    err_match = nan(1, numel(probaBins)-1);
    n_match = nan(1, numel(probaBins)-1);
    y_nonmatch = nan(1, numel(probaBins)-1);
    err_nonmatch = nan(1, numel(probaBins)-1);
    n_nonmatch = nan(1, numel(probaBins)-1);
    for pb = 1:numel(probaBins)-1
        idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ... 
            splitUnitsIdx;% & ...
            % MatchTable.EucledianDistance < 10;
        % all
        n(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx)));
        y(pb) = nanmean(MatchTable.(FPNameCurr)(idx));
        err(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx))/sqrt(n(pb));
        % matches only
        n_match(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx & matchIdx)));
        y_match(pb) = nanmean(MatchTable.(FPNameCurr)(idx & matchIdx));
        err_match(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx & matchIdx))/sqrt(n(pb));
        % non-matches only
        n_nonmatch(pb) = sum(~isnan(MatchTable.(FPNameCurr)(idx & ~matchIdx)));
        y_nonmatch(pb) = nanmean(MatchTable.(FPNameCurr)(idx & ~matchIdx));
        err_nonmatch(pb) = 2*nanstd(MatchTable.(FPNameCurr)(idx & ~matchIdx))/sqrt(n(pb));
    end
    subplot(2,numel(FPNames),fpIdx); hold all
%     errorbar(x,y,err,'k')
%     errorbar(x,y_match,err_match,'r')
%     errorbar(x,y_nonmatch,err_nonmatch,'b')
    plot(x,y,'k.-')
    plot(x,y_match,'r.-')
    plot(x,y_nonmatch,'b.-')
    ylim([0 1])
    xlabel('Matching probability')
    ylabel('Correlation')
    title(FPNameCurr)
%     set(gca,'XScale','log')
    subplot(2,numel(FPNames),numel(FPNames)+fpIdx); hold all
    plot(x,n,'k.-')
    plot(x,n_match,'r.-')
    plot(x,n_nonmatch,'b.-')
    ylim([1 1e9])
    xlabel('Matching probability')
    ylabel('Number of pairs')
    title(FPNameCurr)
    set(gca,'YScale','log')
    yticks([1e2 1e4 1e6 1e8])
end

%% Compute AUCs for each bin

clear AUC
for fpIdx = 1:numel(FPNames)
    FPNameCurr = FPNames{fpIdx};
    for pb = 1:numel(probaBins)-1
        idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2; % & ... 
            % splitUnitsIdx;% & ...
            % MatchTable.EucledianDistance < 10;
        
            if any(idx & ~matchIdx) && any(idx & matchIdx)
                [~, AUC(pb,fpIdx)] = getAUC(MatchTable.(FPNameCurr),find(idx & ~matchIdx),find(idx & matchIdx),0:0.01:1,0);
            end
    end
    [~, AUCall(fpIdx)] = getAUC(MatchTable.(FPNameCurr),find(~matchIdx),find(matchIdx),0:0.01:1,0);
end

figure; hold all
for fpIdx = 1:numel(FPNames)
    plot(x,AUC(:,fpIdx))
    plot([0 1], [AUCall(fpIdx), AUCall(fpIdx)], 'k--')
end
ylim([0 1])
xlabel('Matching probability')
ylabel('AUC')

%%
figure; hold all
% histogram(MatchTable(matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')
pb = 1;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')
pb = numel(probaBins)-1;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).ISICorr,-1:0.05:1,'Normalization','probability')

%%
figure; hold all
pb = 1;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).FlipMatchProb,0:0.03:1,'Normalization','probability')
pb = 2;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
            MatchTable.EucledianDistance < 10;
histogram(MatchTable(idx & ~matchIdx & ~splitUnitsIdx,:).FlipMatchProb,0:0.03:1,'Normalization','probability')


%% Estimate number of false positives and false negatives

%% Get "true" distributions

fpIdx = 1;
FPNameCurr = FPNames{fpIdx};
corrBins = -1:0.1:5;
trackedTrueIdx = matchIdx & ~splitUnitsIdx & MatchTable.MatchProb > 0.99 & MatchTable.RecSes1 ~= MatchTable.RecSes2;
trackedTrueDistr = histcounts(atanh(MatchTable(trackedTrueIdx,:).(FPNameCurr)),corrBins,'Normalization','probability');
differentTrueIdx = ~matchIdx & ~splitUnitsIdx & MatchTable.MatchProb < 1e-5 & MatchTable.RecSes1 ~= MatchTable.RecSes2;
differentTrueDistr = histcounts(atanh(MatchTable(differentTrueIdx,:).(FPNameCurr)),corrBins,'Normalization','probability');

X = [trackedTrueDistr; differentTrueDistr];

%% average case 1 --  all tracked

actualDistr = histcounts(atanh(MatchTable(matchIdx & ~splitUnitsIdx & MatchTable.RecSes1 ~= MatchTable.RecSes2,:).(FPNameCurr)),corrBins,'Normalization','probability');

b = X'\actualDistr';

figure; hold all
stairs(corrBins(1:end-1), trackedTrueDistr, 'r')
stairs(corrBins(1:end-1), differentTrueDistr, 'b')
stairs(corrBins(1:end-1), actualDistr, 'k')
stairs(corrBins(1:end-1), X'*b, 'Color', [.5 .5 .5])

%% average case 2 --  all different

actualDistr = histcounts(atanh(MatchTable(~matchIdx & ~splitUnitsIdx & MatchTable.RecSes1 ~= MatchTable.RecSes2,:).(FPNameCurr)),corrBins,'Normalization','probability');

b = X'\actualDistr';

figure; hold all
stairs(corrBins(1:end-1), trackedTrueDistr, 'r')
stairs(corrBins(1:end-1), differentTrueDistr, 'b')
stairs(corrBins(1:end-1), actualDistr, 'k')
stairs(corrBins(1:end-1), X'*b, 'Color', [.5 .5 .5])

%% extreme case 1 -- prob = 0, tracked

pb = 1;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & matchIdx & ~splitUnitsIdx ;
actualDistr = histcounts(atanh(MatchTable(idx,:).(FPNameCurr)),corrBins,'Normalization','probability');

b = X'\actualDistr';

figure; hold all
stairs(corrBins(1:end-1), trackedTrueDistr, 'r')
stairs(corrBins(1:end-1), differentTrueDistr, 'b')
stairs(corrBins(1:end-1), actualDistr, 'k')
stairs(corrBins(1:end-1), X'*b, 'Color', [.5 .5 .5])

figure;
histogram(abs(MatchTable(differentTrueIdx,:).RecSes2 - MatchTable(differentTrueIdx,:).RecSes1),'Normalization','probability')
hold all
histogram(abs(MatchTable(idx,:).RecSes2 - MatchTable(idx,:).RecSes1),'Normalization','probability')

%% extreme case 2 -- prob = 1, different

pb = numel(probaBins)-1;
idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ~matchIdx & ~splitUnitsIdx;
actualDistr = histcounts(atanh(MatchTable(idx,:).(FPNameCurr)),corrBins,'Normalization','probability');

b = X'\actualDistr';

figure; hold all
stairs(corrBins(1:end-1), trackedTrueDistr, 'r')
stairs(corrBins(1:end-1), differentTrueDistr, 'b')
stairs(corrBins(1:end-1), actualDistr, 'k')
stairs(corrBins(1:end-1), X'*b, 'Color', [.5 .5 .5])

figure;
histogram(abs(MatchTable(trackedTrueIdx,:).RecSes2 - MatchTable(trackedTrueIdx,:).RecSes1),'Normalization','probability')
hold all
histogram(abs(MatchTable(idx,:).RecSes2 - MatchTable(idx,:).RecSes1),'Normalization','probability')

%% across all proba bins

fpIdx = 1;
for pb = 1:numel(probaBins)-1
    idx = MatchTable.MatchProb >= probaBins(pb) & MatchTable.MatchProb < probaBins(pb+1) & MatchTable.RecSes1 ~= MatchTable.RecSes2 & ...
    splitUnitsIdx;% & ...
    % MatchTable.EucledianDistance < 10;
    % all
    actualDistr = histcounts(atanh(MatchTable(idx,:).(FPNameCurr)),corrBins,'Normalization','probability');
    b_all(:,pb) = X'\actualDistr';
    % matches only
    actualDistr = histcounts(atanh(MatchTable(idx & matchIdx,:).(FPNameCurr)),corrBins,'Normalization','probability');
    b_match(:,pb) = X'\actualDistr';
    n_match_TP(pb) = sum(idx & matchIdx) * b_match(1,pb);
    n_match_FP(pb) = sum(idx & matchIdx) * b_match(2,pb);
    % non-matches only
    actualDistr = histcounts(atanh(MatchTable(idx & ~matchIdx,:).(FPNameCurr)),corrBins,'Normalization','probability');
    b_nonmatch(:,pb) = X'\actualDistr';
    n_nonmatch_TN(pb) = sum(idx & ~matchIdx) * b_nonmatch(2,pb);
    n_nonmatch_FN(pb) = sum(idx & ~matchIdx) * b_nonmatch(1,pb);
end

figure;
subplot(131); hold all
stairs(probaBins(1:end-1), b_all(1,:),'r')
stairs(probaBins(1:end-1), b_all(2,:),'b')
legend({'same','different'})
xlabel('Match probability')
ylabel('%')
title('all')
ylim([0 1])
subplot(132); hold all
stairs(probaBins(1:end-1), b_match(1,:),'r')
stairs(probaBins(1:end-1), b_match(2,:),'b')
title('tracked')
ylim([0 1])
subplot(133); hold all
stairs(probaBins(1:end-1), b_nonmatch(1,:),'r')
stairs(probaBins(1:end-1), b_nonmatch(2,:),'b')
title('different')
ylim([0 1])