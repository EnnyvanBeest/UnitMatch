function AUCExtract(UMfiles)
TrajectoryScores = nan(2,2,length(UMfiles));
for umid = 1:length(UMfiles)
    load(strrep(UMfiles{umid},'UnitMatch.mat','AUC.mat'),'AUCStruct')
    TrajectoryScores(:,:,umid) = AUCStruct.TrajectoryScore;

end

figure('name','TrajectoryScoreOptions')
barwitherr(nanstd(TrajectoryScores,[],3),nanmean(TrajectoryScores,3))
ylabel('AUC value within recording')
set(gca,'XTickLabel',{'AngleOld','DotProductAngle'})
legend({'Subtraction','Correlation'})
makepretty
offsetAxes




return