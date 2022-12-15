function [Parameterkernels, Performance] = CreateNaiveBayes(Tbl,label)
% Table with a row per parameter that will be used as a predictor
% label attribute each to the correct group (0 or 1)
ncross = 5;

Scores2Include = Tbl.Properties.VariableNames;
stepsize = 0.01;
Edges = [0:stepsize:1];
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
Parameterkernels = nan(length(ScoreVector),length(Scores2Include),2); %These can be used as kernels in the Naive Bayes

%% Cross validate
shuffleddata = datasample(1:length(label),length(label),'Replace',false);
sampleperfold = floor(length(shuffleddata)./ncross);
% figure('name','UsedKernels')
Performance = nan(2,ncross);
for cv = 1:ncross
    % Get test and train data
    testidx = shuffleddata((cv-1)*sampleperfold+1:cv*sampleperfold);
    trainidx = shuffleddata;
    trainidx((cv-1)*sampleperfold+1:cv*sampleperfold) = [];
    for scid=1:length(Scores2Include)
        eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
        ScoresTmp = ScoresTmp(trainidx);
        % Save out probability distribution
        Parameterkernels(:,scid,1) = histcounts(ScoresTmp(~label(trainidx)),Edges)./sum(~label(trainidx));
        Parameterkernels(:,scid,2) = histcounts(ScoresTmp(label(trainidx)),Edges)./sum(label(trainidx));

%         subplot(5,length(Scores2Include),(cv-1)*length(Scores2Include)+scid)
%         plot(ScoreVector,Parameterkernels(:,scid,1)); hold on
%         plot(ScoreVector,Parameterkernels(:,scid,2));
%         xlabel('Score')
%         ylabel('p|observation')
%         title(Scores2Include{scid})
%         makepretty
    end

    % Apply
    [labelcv, probability, Performance(:,cv)] = ApplyNaiveBayes(Tbl(testidx,:),Parameterkernels,label(testidx));
    
end
Performance = round(nanmean(Performance,2)*1000)/10;
disp(['Average cross-validation match = 1 performance: ' num2str(Performance(1)) '%, match = 0 performance: ' num2str(Performance(2))])
%% Final parameters 
figure('name','UsedKernels')
for scid=1:length(Scores2Include)
    eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
    % Save out probability distribution
    Parameterkernels(:,scid,1) = histcounts(ScoresTmp(~label),Edges)./sum(~label);    
    Parameterkernels(:,scid,2) = histcounts(ScoresTmp(label),Edges)./sum(label);
    
    subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),scid)
    plot(ScoreVector,Parameterkernels(:,scid,1)); hold on
    plot(ScoreVector,Parameterkernels(:,scid,2));
    xlabel('Score')
    ylabel('p|observation')
    title(Scores2Include{scid})
    makepretty
end
return

