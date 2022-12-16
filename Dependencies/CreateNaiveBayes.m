function [Parameterkernels, Performance] = CreateNaiveBayes(Tbl,label)
% Table with a row per parameter that will be used as a predictor
% label attribute each to the correct group (0 or 1)
%% better not to change these, unless you change them everywhere
ncross = 5;
Scores2Include = Tbl.Properties.VariableNames;
stepsize = 0.01;


%% inatialization
Cond = unique(label); %Number of conditions
Edges = [0:stepsize:1];
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
Parameterkernels = nan(length(ScoreVector),length(Scores2Include),length(Cond)); %These can be used as kernels in the Naive Bayes

%% Cross validate
shuffleddata = datasample(1:length(label),length(label),'Replace',false);
sampleperfold = floor(length(shuffleddata)./ncross);
Performance = nan(2,ncross);
for cv = 1:ncross
    % Get test and train data
    testidx = shuffleddata((cv-1)*sampleperfold+1:cv*sampleperfold);
    trainidx = shuffleddata;
    trainidx((cv-1)*sampleperfold+1:cv*sampleperfold) = [];

    % Create probability distributions based on train data
    for scid=1:length(Scores2Include)
        eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
        ScoresTmp = ScoresTmp(trainidx);
        % Save out probability distribution
        for Ck = 1:length(Cond)
            Parameterkernels(:,scid,Ck) = smooth(histcounts(ScoresTmp(label(trainidx)==Cond(Ck)),Edges),10)./sum(label(trainidx)==Cond(Ck));
        end
    end
    % Apply
    [~, ~, Performance(:,cv)] = ApplyNaiveBayes(Tbl(testidx,:),Parameterkernels,label(testidx));
    
end
Performance = round(nanmean(Performance,2)*1000)/10;
for Ck=1:length(Cond)
    disp(['Average cross-validation match = ' num2str(Cond(Ck)) ' performance: ' num2str(Performance(Ck)) '%'])
end
%% Final parameters 
figure('name','UsedKernels')
for scid=1:length(Scores2Include)
    % Get data
    eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])

    % Prepare plot
    subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),scid)
    hold on

    % Save out probability distribution
    for Ck = 1:length(Cond)
        Parameterkernels(:,scid,Ck) = smooth(histcounts(ScoresTmp(label==Cond(Ck)),Edges),10)./sum(label==Cond(Ck));
        plot(ScoreVector,Parameterkernels(:,scid,Ck));
    end

    % Make plot prettier
    xlabel('Score')
    ylabel('p|observation')
    title(Scores2Include{scid})
    makepretty
end
return

