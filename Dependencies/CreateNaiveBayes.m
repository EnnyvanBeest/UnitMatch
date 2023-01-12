function [Parameterkernels, Performance] = CreateNaiveBayes(Tbl,label,Prior)
% Table with a row per parameter that will be used as a predictor
% label attribute each to the correct group (0 or 1)
%% better not to change these, unless you change them everywhere
ncross = 5;
Scores2Include = Tbl.Properties.VariableNames;
SmoothProbability = 5; % smoothing, since we need to make sure none of the scores have a probability of 0
SmoothSpatialFactor = 5; % Apply another smoothing factor for spatial decay and space
addone=1; %Add one to the distribution to not have a probability of 0
if nargin<3
    Prior = [0.5 0.5];
end
global stepsize

%% inatialization
Cond = unique(label); %Number of conditions
Edges = [0:stepsize:1];
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
Parameterkernels = nan(length(ScoreVector),length(Scores2Include),length(Cond)); %These can be used as kernels in the Naive Bayes

%% Final parameters 
figure('name','UsedKernels')
for scid=1:length(Scores2Include)
    % Get data
    eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])

    if ismember(Scores2Include{scid},{'LocDistSim','spatialdecaySim'})
        Smoothtmp = SmoothSpatialFactor;
    else
        Smoothtmp = SmoothProbability;
    end

    % Prepare plot
    subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),scid)
    hold on

    % Scores irrespective of label:
%     dlim = quantile(ScoresTmp,0.025);
%     ulim = quantile(ScoresTmp,0.975);
    % Save out probability distribution
    for Ck = 1:length(Cond)
        Parameterkernels(:,scid,Ck) = smooth(histcounts(ScoresTmp(label==Cond(Ck)),Edges)+addone,Smoothtmp)./(sum(label==Cond(Ck))+addone);
     
%         % modify the distribution for not frequently observed values
%         Parameterkernels(1:find(dlim-Edges<0,1,'first'),scid,Ck) = Parameterkernels(find(dlim-Edges<0,1,'first'),scid,Ck);
%         Parameterkernels(find(ulim-Edges<0,1,'first'):end,scid,Ck) = Parameterkernels(find(ulim-Edges<0,1,'first'),scid,Ck);

        %Sum should still be 1
        plot(ScoreVector,Parameterkernels(:,scid,Ck));
    end

    % Make plot prettier
    xlabel('Score')
    ylabel('p|observation')
    title(Scores2Include{scid})
    makepretty
end

%% Cross validate
shuffleddata = datasample(1:length(label),length(label),'Replace',false);
sampleperfold = floor(length(shuffleddata)./ncross);
Performance = nan(2,ncross);
ParameterkernelsCV = nan(size(Parameterkernels));
for cv = 1:ncross
    % Get test and train data
    testidx = shuffleddata((cv-1)*sampleperfold+1:cv*sampleperfold);
    trainidx = shuffleddata;
    trainidx((cv-1)*sampleperfold+1:cv*sampleperfold) = [];

    % Create probability distributions based on train data
    for scid=1:length(Scores2Include)
        eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
        if strcmp(Scores2Include,'LocDistSim')
            Smoothtmp = SmoothSpatialFactor;
        else
            Smoothtmp = SmoothProbability;
        end

        ScoresTmp = ScoresTmp(trainidx);
        % Scores irrespective of label:
%         dlim = quantile(ScoresTmp,0.025);
%         ulim = quantile(ScoresTmp,0.975);

        % Save out probability distribution
        for Ck = 1:length(Cond)
            % Always add 1 count!!
            ParameterkernelsCV(:,scid,Ck) = smooth(histcounts(ScoresTmp(label(trainidx)==Cond(Ck)),Edges)+addone,Smoothtmp)./(sum(label(trainidx)==Cond(Ck))+addone);

            % modify the distribution for not frequently observed values
%             ParameterkernelsCV(1:find(dlim-Edges<0,1,'first'),scid,Ck) = ParameterkernelsCV(find(dlim-Edges<0,1,'first'),scid,Ck);
%             ParameterkernelsCV(find(ulim-Edges<0,1,'first'):end,scid,Ck) = ParameterkernelsCV(find(ulim-Edges<0,1,'first'),scid,Ck);

        end
    end
    % Apply
    [~, ~, Performance(:,cv)] = ApplyNaiveBayes(Tbl(testidx,:),ParameterkernelsCV,label(testidx),Prior);
    
end
Performance = round(nanmean(Performance,2)*1000)/10;
for Ck=1:length(Cond)
    disp(['Average cross-validation match = ' num2str(Cond(Ck)) ' performance: ' num2str(Performance(Ck)) '%'])
end
return

