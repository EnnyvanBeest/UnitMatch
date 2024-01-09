function [MatchProbability,label,Tbl,BestMdl] = RunNaiveBayes(Predictors,TotalScore,Scores2Include,clusinfo,param,SortingOrder,EuclDist)

%% Extract parameters
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')

end
GoodRecSesID = clusinfo.RecSesID(Good_Idx);
recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
ndays = length(unique(recsesGood));
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),unique(recsesGood),'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch); nclus+1];
nCellsPerSession = diff(SessionSwitch);

%% Prepare naive bayes - inspect probability distributions
IncludeThesePairs = find(EuclDist<param.maxdist);
% Prepare a set INCLUDING the cross-validated self-scores, otherwise the probability
% distributions are just weird
priorMatch = 1-(param.nExpectedMatches./length(IncludeThesePairs)); %

% priorMatch = 1-((nclus+nclus.*sqrt(ndays-1)*2*param.ExpectMatches)./length(IncludeThesePairs)); %Punish multiple days (unlikely to find as many matches after a few days) *2 for symmetry
% priorMatch = 1-(nclus*ndays)./(nclus*nclus); %Now use the actual expected prior for bayes'
ThrsOpt = quantile(TotalScore(IncludeThesePairs),priorMatch);
CandidatePairs = TotalScore>ThrsOpt;%

figure('name','Potential Matches')
imagesc(CandidatePairs(SortingOrder,SortingOrder))
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
%     xlim([SessionSwitch nclus])
%     ylim([1 SessionSwitch-1])
xlabel('Unit_i')
ylabel('Unit_j')
title('Potential Matches')
makepretty
[uid,uid2] = find(CandidatePairs);
Pairs = cat(2,uid,uid2);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs,'rows');
global stepsize

%% Initialize
% Usually this means there's no variance in the match distribution
% (which in a way is great). Create some small variance
flag = 0;
MinLoss=1;
MaxPerf = [0 0];
runid = 0;
% Priors = [0.5 0.5];

priorMatch = 1-(param.nExpectedMatches./nclus.^2) %
Priors = [priorMatch 1-priorMatch];
BestMdl = [];
Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
while flag<2 && runid<param.maxrun
    timercounter = tic;

    flag = 0;
    runid = runid+1;
    filedir = which('UnitMatch');
    filedir = dir(filedir);
    if param.ApplyExistingBayesModel && exist(fullfile(filedir.folder,'UnitMatchModel.mat'))
        disp('Loading the existing Naive Bayes model...')

        load(fullfile(filedir.folder,'UnitMatchModel.mat'),'BestMdl')
        % Apply naive bays classifier
        if isfield(BestMdl,'Parameterkernels')
            [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels(:,ismember(BestMdl.VariableNames,Scores2Include),:),[0 1],Priors);
        else
            [label, posterior, cost] = predict(BestMdl,Tbl);
        end
    else
        tmp = reshape(Predictors,[],length(Scores2Include));
        Tbltmp = array2table(reshape(tmp,[],size(tmp,2)),'VariableNames',Scores2Include); %All parameters
        label = CandidatePairs(:);

        if param.MakeOwnNaiveBayes
            % Work in progress
            fprintf('Creating the Naive Bayes model...\n')
            timercounter2 = tic;
            [Parameterkernels,Performance] = CreateNaiveBayes(Tbltmp,label,Priors);
            fprintf('Done in %ds.\n', round(toc(timercounter2)))
            if any(Performance'<MaxPerf)
                flag = flag+1;
            else
                BestMdl.Parameterkernels = Parameterkernels;
            end
            % Apply naive bays classifier
            [label, posterior] = ApplyNaiveBayes(Tbl,Parameterkernels,[0 1],Priors);
            saveas(gcf,fullfile(param.SaveDir,'ProbabilityDistribution.fig'))
            saveas(gcf,fullfile(param.SaveDir,'ProbabilityDistribution.bmp'))
        else % This uses matlab package. Warning: normal distributions assumed?
            try
                Mdl = fitcnb(Tbltmp,label);
            catch ME
                disp(ME)
                keyboard
                for id = 1:size(Predictors,3)
                    tmp = Predictors(:,:,id);
                    if nanvar(tmp(CandidatePairs(:)==1)) == 0
                        %Add some noise
                        tmp(CandidatePairs(:)==1) = tmp(CandidatePairs(:)==1)+(rand(sum(CandidatePairs(:)==1),1)-0.5)./2;
                        tmp(tmp>1)=1;
                        Predictors(:,:,id)=tmp;
                    end
                end
            end
            % Cross validate on model that uses only prior
            DefaultPriorMdl = Mdl;
            FreqDist = cell2table(tabulate(label==1));
            DefaultPriorMdl.Prior = FreqDist{:,3};
            rng(1);%
            defaultCVMdl = crossval(DefaultPriorMdl);
            defaultLoss = kfoldLoss(defaultCVMdl);

            CVMdl = crossval(Mdl);
            Loss = kfoldLoss(CVMdl);

            if Loss>defaultLoss
                warning('Model doesn''t perform better than chance')
            end
            if round(Loss*10000) >= round(MinLoss*10000)
                flag = flag+1;
            elseif Loss<MinLoss
                MinLoss=Loss;
                BestMdl = Mdl;
            end
            disp(['Loss = ' num2str(round(Loss*10000)/10000)])

            % Apply naive bays classifier
            [label, posterior, cost] = predict(Mdl,Tbl);

            %% Evaluate Model:
            figure('name','NaiveBayesEstimates')
            for parid=1:size(Predictors,3)

                subplot(size(Predictors,3),1,parid)
                mu = BestMdl.DistributionParameters{1,parid}(1);
                sigma = BestMdl.DistributionParameters{1,parid}(2);
                x = (-5 * sigma:0.01:5*sigma)+mu;
                plot(x,normpdf(x,mu,sigma),'b-')
                hold on
                mu = BestMdl.DistributionParameters{2,parid}(1);
                sigma = BestMdl.DistributionParameters{2,parid}(2);
                x = (-5 * sigma:0.01:5*sigma)+mu;
                plot(x,normpdf(x,mu,sigma),'r-')
                title(Scores2Include{parid})

                makepretty
                xlim([0 1])
            end
            saveas(gcf,fullfile(param.SaveDir,'ProbabilityDistribution.fig'))
            saveas(gcf,fullfile(param.SaveDir,'ProbabilityDistribution.bmp'))

        end
    end
    drawnow

    disp(['Getting the Naive Bayes model took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
end

%% Extract final pairs:
disp('Extracting final pairs of units...')
timercounter = tic;

if isfield(BestMdl,'Parameterkernels')
    BestMdl.VariableNames = Scores2Include;
    BestMdl.Parameterkernels = BestMdl.Parameterkernels(:,ismember(BestMdl.VariableNames,Scores2Include),:);
    Edges = [0:stepsize:1];
    ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
    BestMdl.ScoreVector = ScoreVector;
    [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,[0 1],Priors);
else
    [label, posterior, cost] = predict(BestMdl,Tbl);
end
MatchProbability = reshape(posterior(:,2),size(Predictors,1),size(Predictors,2));
label = (MatchProbability>=param.ProbabilityThreshold);% | (MatchProbability>0.05 & RankScoreAll==1 & SigMask==1);
% label = reshape(label,nclus,nclus);
[r, c] = find(label==1); %Find matches across 2 days
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs,'rows');

figure;
subplot(2,2,1)
imagesc(MatchProbability(SortingOrder,SortingOrder),[0 1])
colormap(flipud(gray))
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
title('match probability')
makepretty

subplot(2,2,2)
imagesc(MatchProbability(SortingOrder,SortingOrder)>0.5)
colormap(flipud(gray))
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
title('thresholded')
makepretty


subplot(2,2,3)
tmp = MatchProbability;
% Take centroid dist > maxdist out
tmp(EuclDist>param.NeighbourDist)=nan;
% Take between session out
for did = 1:ndays
    for did2 = 1:ndays
        if did==did2
            continue
        end
        tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did2):SessionSwitch(did2+1)-1)=nan;
    end
end
stepsz = 0.01;

hd = histcounts(diag(tmp),0:stepsz:1);%./nclus + 0.0001;
hnd = histcounts(tmp(~eye(size(tmp))),0:stepsz:1);%./sum(~isnan(tmp(~eye(size(tmp))))) +0.0001;
plot(stepsz./2:stepsz:1-stepsz./2,hd,'-','color',[0 0.7 0]); hold on; plot(stepsz./2:stepsz:1-stepsz./2,hnd,'b-')
tmp = MatchProbability;
% Take centroid dist > maxdist out
tmp(EuclDist>param.NeighbourDist)=nan;
% Take within session out
for did = 1:ndays
    tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)=nan;
end
ha = histcounts(tmp(:),0:stepsz:1);%./sum(~isnan(tmp(:))) + 0.0001;
plot(stepsz./2:stepsz:1-stepsz./2,ha,'-','color',[1 0 0]);
set(gca,'yscale','linear')
line([param.ProbabilityThreshold param.ProbabilityThreshold ],get(gca,'ylim'),'LineStyle','--','color',[0 0 0])
xlabel('MatchProbability')
ylabel('number of pairs')
makepretty
title('P(match) distributions')
legend('Same Unit','Neighbors','Across','Threshold','Location', 'best')

%Cumulative density
subplot(2,2,4)
hold on
tmp = MatchProbability;
% Take centroid dist > maxdist out
tmp(EuclDist>param.NeighbourDist)=nan;
% Take between session out
for did = 1:ndays
    for did2 = 1:ndays
        if did==did2
            continue
        end
        tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did2):SessionSwitch(did2+1)-1)=nan;
    end
end
[h,stats] = cdfplot(diag(tmp)); %same
h.Color = [0 0.5 0];

tmp(logical(eye(size(tmp)))) = nan;

[h,stats] = cdfplot(tmp(~eye(size(tmp)))); %Neighbours
h.Color = [0 0 0.5];

tmp = MatchProbability;
% Take centroid dist > maxdist out
tmp(EuclDist>param.NeighbourDist)=nan;
% Take within session out
for did = 1:ndays
    tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)=nan;
end
if any(~isnan(tmp(:)))
[h,stats] = cdfplot(tmp(:)); %Across days
h.Color = [1 0 0];
end


xlabel('MatchProbability')
ylabel('Cumulative density')
line([param.ProbabilityThreshold,param.ProbabilityThreshold],[0 1],'color',[1 0 0],'LineStyle','--')
makepretty

legend('Same unit','Neighbors','Across','threshold')

saveas(gcf,fullfile(param.SaveDir,'MatchProbability.fig'))
saveas(gcf,fullfile(param.SaveDir,'MatchProbability.bmp'))

%%
BestMdl.Priors = Priors;
disp(['Extracting final pair of units took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
