function [Parameterkernels, Performance] = CreateNaiveBayes(Tbl,label,Prior)
% Table with a row per parameter that will be used as a predictor
% label attribute each to the correct group (0 or 1)
%% better not to change these, unless you change them everywhere
ncross = 5;
Scores2Include = Tbl.Properties.VariableNames;
SmoothProbability = 10; % smoothing, since we need to make sure none of the scores have a probability of 0
SmoothSpatialFactor = 10; % Apply another smoothing factor for spatial decay and space
global stepsize

addone=1; %Add min probability to the distribution to not have a probability of 0 with small sample size(default 0)
if nargin<3
    Prior = [0.5 0.5];
end

%% inatialization
Cond = unique(label); %Number of conditions
Edges = 0:stepsize:1;
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
Parameterkernels = nan(length(ScoreVector),length(Scores2Include),length(Cond)); %These can be used as kernels in the Naive Bayes

%% Final parameters 
figure('name','UsedKernels')
for scid=1:length(Scores2Include)
    % Get data
    ScoresTmp = Tbl.(Scores2Include{scid});

    if ismember(Scores2Include{scid},{'LocDistSim','spatialdecaySim'})
        Smoothtmp = SmoothSpatialFactor;
    else
        Smoothtmp = SmoothProbability;
    end

    % Prepare plot
    subplot(1,length(Scores2Include),scid)
    hold on

    % Scores irrespective of label:
%     dlim = quantile(ScoresTmp,0.025);
%     ulim = quantile(ScoresTmp,0.975);
    % Save out probability distribution
    for Ck = 1:length(Cond)
        % Probability Distributions
        Parameterkernels(:,scid,Ck) = smooth(histcounts(ScoresTmp(label==Cond(Ck)),Edges),Smoothtmp);
        Parameterkernels(:,scid,Ck) = Parameterkernels(:,scid,Ck)/sum(Parameterkernels(:,scid,Ck)); % normalize to 1
     
        % Avoid having 0 probability
        Parameterkernels(:,scid,Ck) = Parameterkernels(:,scid,Ck)+addone*min(Parameterkernels(Parameterkernels(:,scid,Ck)~=0,scid,Ck),[],1);
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
        ScoresTmp = Tbl.(Scores2Include{scid});
        
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
            % Probability Distributions
            ParameterkernelsCV(:,scid,Ck) = smooth(histcounts(ScoresTmp(label(trainidx)==Cond(Ck)),Edges),Smoothtmp);
            ParameterkernelsCV(:,scid,Ck) = ParameterkernelsCV(:,scid,Ck)/sum(ParameterkernelsCV(:,scid,Ck)); % normalize to 1

            % Avoid having 0 probability
            ParameterkernelsCV(:,scid,Ck) = ParameterkernelsCV(:,scid,Ck)+addone*min(ParameterkernelsCV(Parameterkernels(:,scid,Ck)~=0,scid,Ck),[],1);
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

