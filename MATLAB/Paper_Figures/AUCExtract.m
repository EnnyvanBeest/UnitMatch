function AUCExtract(UMfiles)
TrajectoryScores = nan(2,2,length(UMfiles));
for umid = 1:length(UMfiles)
    load(strrep(UMfiles{umid},'UnitMatch.mat','AUC.mat'),'AUCStruct')
    % TrajectoryScores(:,:,umid) = AUCStruct.TrajectoryScore;

    if umid == 1
        ScoreItems = AUCStruct.ParamNames;
        AUCScores = nan(numel(ScoreItems),numel(UMfiles));
    end
    AUCScores(ismember(ScoreItems,AUCStruct.ParamNames),umid) = AUCStruct.AUC(ismember(AUCStruct.ParamNames,ScoreItems));

end
[~,sortidx] = sort(nanmean(AUCScores,2),'descend');
ExcludeNans = all(isnan(AUCScores),2);
sortidx(ExcludeNans) = [];

% figure('name','TrajectoryScoreOptions')
% barwitherr(nanstd(TrajectoryScores,[],3),nanmean(TrajectoryScores,3))
% ylabel('AUC value within recording')
% set(gca,'XTickLabel',{'AngleOld','DotProductAngle'})
% legend({'Subtraction','Correlation'})
% makepretty
% offsetAxes

figure('name','AUC per Score')
violinplot(AUCScores(sortidx,:)')
ylabel('AUC value within recording')
set(gca,'XTickLabel',ScoreItems(sortidx))
makepretty
offsetAxes



return