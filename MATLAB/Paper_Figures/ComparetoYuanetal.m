if 0
%% Within day comparison
datapath = 'H:\MatchingUnits\Yuan\WithinDay';
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.9;

FPFNYuan = nan(length(miceopt),2); %FP/FN
FPFNUM = nan(length(miceopt),2); %FP/FN

% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));

    ResultTable = tmpYuan.output.results_wth(:,[2,3]); % This is with their 10 micron cut-off
    % ResultTable = tmpYuan.output.all_results_post(:,[2,3]);

    % Proportion - only for within day
    FPFNYuan(midx,2) = (tmpYuan.output.KSgood_f1 - sum(ResultTable(:,1)==ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were no found?
    FPFNYuan(midx,1) = (sum(ResultTable(:,1)~=ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were found extra?

    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapath,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    if nclus ~= tmpYuan.output.KSgood_f1
        disp('Not same amount of units.. Unfair comparison??')
        keyboard
    end
    MatchProbability = reshape(tmpUM.MatchTable.MatchProb,nclus,nclus);
    EuclDist = reshape(tmpUM.MatchTable.EucledianDistance,nclus,nclus);

    % Also use distance threshold that was used in running UnitMatch
    MatchProbability(EuclDist>tmpUM.UMparam.maxdist) = 0;
    AverageMP = nanmean(cat(3,MatchProbability,MatchProbability'),3);

    % proportion
    FPFNUM(midx,2) = (nclus - sum(AverageMP(logical(eye(nclus)))>ProbThresh))./nclus;
    FPFNUM(midx,1) = sum(AverageMP(~logical(eye(nclus)))>ProbThresh)./nclus;
end

% Histogram
figure; barwitherr(squeeze(nanstd(cat(3,FPFNYuan,FPFNUM),[],1))./sqrt(length(miceopt)-1),squeeze(nanmean(cat(3,FPFNYuan,FPFNUM),1)));
set(gca,'XTickLabel',{'False Positives','False Negatives'})
ylabel('Proportion (of total neurons)')
legend({'Yuan et al.',['UnitMatch, p=' num2str(ProbThresh)]})
title('Within-day comparison')

makepretty
offsetAxes
end
%% Across days
datapath = 'H:\MatchingUnits\Yuan\AcrossDays';
miceopt = {'AL032','AV008','CB016','EB019','JF067'};
ProbThresh = 0.5;

FPFNYuan = nan(length(miceopt),2); %FP/FN
FPFNUM = nan(length(miceopt),2); %FP/FN

% Data was split up first versus second half of recording. How many of good
% units (Bombcell isolated) are found back by the two algorithms?
for midx = 1:length(miceopt)
    % What about UnitMatch?
    UMOutputFile = dir(fullfile(datapath,miceopt{midx},'**','UnitMatch.mat')); % Find output file
    tmpUM = load(fullfile(UMOutputFile.folder,UMOutputFile.name));

    nclus = size(tmpUM.WaveformInfo.MaxChannel,1);
    if nclus ~= tmpYuan.output.KSgood_f1
        disp('Not same amount of units.. Unfair comparison??')
        keyboard
    end
    MatchProbability = reshape(tmpUM.MatchTable.MatchProb,nclus,nclus);
    EuclDist = reshape(tmpUM.MatchTable.EucledianDistance,nclus,nclus);

    % Also use distance threshold that was used in running UnitMatch
    MatchProbability(EuclDist>tmpUM.UMparam.maxdist) = 0;
    AverageMP = nanmean(cat(3,MatchProbability,MatchProbability'),3);

    % proportion
    FPFNUM(midx,2) = (nclus - sum(AverageMP(logical(eye(nclus)))>ProbThresh))./nclus;
    FPFNUM(midx,1) = sum(AverageMP(~logical(eye(nclus)))>ProbThresh)./nclus;
    
    YuanOutputFile = dir(fullfile(datapath,miceopt{midx},'**','Output.mat')); % Find output file
    tmpYuan = load(fullfile(YuanOutputFile.folder,YuanOutputFile.name));

    % ResultTable = tmpYuan.output.results_wth(:,[2,3]);
    ResultTable = tmpYuan.output.all_results_post(:,[2,3]);

    % Proportion - only for within day
    FPFNYuan(midx,2) = (tmpYuan.output.KSgood_f1 - sum(ResultTable(:,1)==ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were no found?
    FPFNYuan(midx,1) = (sum(ResultTable(:,1)~=ResultTable(:,2)))./tmpYuan.output.KSgood_f1; % How many of the neurons were found extra?


end
