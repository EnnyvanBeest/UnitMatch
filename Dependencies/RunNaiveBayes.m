 function [MatchProbability,label,Pairs,Tbl,BestMdl] = RunNaiveBayes(Predictors,TotalScore,Scores2Include,RankScoreAll,SigMask,clusinfo,param)
%% Extract parameters
Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
GoodRecSesID = clusinfo.RecSesID(Good_Idx);
OriginalClusterIDs = clusinfo.cluster_id;

recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
ndays = length(unique(recsesAll));
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),1:ndays,'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch) nclus+1];

%% Prepare naive bayes - inspect probability distributions
% Prepare a set INCLUDING the cross-validated self-scores, otherwise the probability
% distributions are just weird
priorMatch = 1-(nclus*ndays)./(nclus*nclus); %Now use a slightly more lenient prior
ThrsOpt = quantile(TotalScore(:),priorMatch);
CandidatePairs = TotalScore>ThrsOpt & RankScoreAll == 1 & SigMask == 1; %Best results if all 3 are true
% Also assume kilosort does a good job ?
CandidatePairs(logical(eye(size(CandidatePairs))))=1;

figure('name','Potential Matches')
imagesc(CandidatePairs)
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
npairs = 0;
MinLoss=1;
MaxPerf = [0 0];
npairslatest = 0;
runid = 0;
% Priors = [0.5 0.5];
Priors = [priorMatch 1-priorMatch];
BestMdl = [];
Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters

while flag<2 && runid<param.maxrun

    timercounter = tic;
    disp('Getting the Naive Bayes model...')

    flag = 0;
    runid = runid+1;
    filedir = which('UnitMatch');
    filedir = dir(filedir);
    if param.ApplyExistingBayesModel && exist(fullfile(filedir.folder,'UnitMatchModel.mat'))
        load(fullfile(param.SaveDir,'UnitMatchModel.mat'),'BestMdl')
        % Apply naive bays classifier

        if isfield(BestMdl,'Parameterkernels')
            [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,[0 1],Priors);
        else
            [label, posterior, cost] = predict(BestMdl,Tbl);
        end
    else
        tmp= reshape(Predictors(Pairs(:,1),Pairs(:,2),:),[],length(Scores2Include));
        Tbltmp = array2table(reshape(tmp,[],size(tmp,2)),'VariableNames',Scores2Include); %All parameters
        % Use Rank as 'correct' label
        label = reshape(CandidatePairs(Pairs(:,1),Pairs(:,2)),1,[])';
       
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

    if runid<param.maxrun % Otherwise a waste of time!
        label = reshape(label,size(Predictors,1),size(Predictors,2));
        [r, c] = find(triu(label)==1); %Find matches

        Pairs = cat(2,r,c);
        Pairs = sortrows(Pairs);
        Pairs = unique(Pairs,'rows');
        %     Pairs(Pairs(:,1)==Pairs(:,2),:)=[];
        MatchProbability = reshape(posterior(:,2),size(Predictors,1),size(Predictors,2));
        %     figure; imagesc(label)

        % Functional score for optimization: compute Fingerprint for the matched units - based on Célian Bimbard's noise-correlation finger print method but applied to across session correlations
        % Not every recording day will have the same units. Therefore we will
        % correlate each unit's activity with average activity across different
        % depths
        disp('Recalculate activity correlations')

        % Use a bunch of units with high total scores as reference population
        [PairScore,sortid] = sort(cell2mat(arrayfun(@(X) MatchProbability(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
        Pairs = Pairs(sortid,:);
        [FingerprintR,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(sessionCorrelationsAll,Pairs,Unit2Take,recsesGood);
%         CrossCorrelationFingerPrint_BU

        tmpf = triu(FingerprintR,1);
        tmpm = triu(MatchProbability,1);
        tmpm = tmpm(tmpf~=0);
        tmpf = tmpf(tmpf~=0);
        tmpr = triu(RankScoreAll,1);
        tmpr = tmpr(tmpr~=0);

        figure;
        scatter(tmpm,tmpf,14,tmpr,'filled')
        colormap(cat(1,[0 0 0],winter))
        xlabel('Match Probability')
        ylabel('Cross-correlation fingerprint')
        makepretty
        drawnow

        % New Pairs for new round
        CandidatePairs = label==1 & RankScoreAll==1& SigMask==1;
        CandidatePairs(tril(true(size(CandidatePairs)),-1))=0;
        [uid,uid2] = find(CandidatePairs);
        Pairs = cat(2,uid,uid2);
        Pairs = sortrows(Pairs);
        Pairs=unique(Pairs,'rows');
    end
    disp(['Getting the Naive Bayes model took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
end

%% Extract final pairs:
disp('Extracting final pairs of units...')
timercounter = tic;
   
if isfield(BestMdl,'Parameterkernels')
    BestMdl.VariableNames = Scores2Include;
    Edges = [0:stepsize:1];
    ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
    BestMdl.ScoreVector = ScoreVector;
    if param.RunPyKSChronicStitched
        [label, posterior,performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,PyKSLabel(:),Priors);
        disp(['Correctly labelled ' num2str(round(performance(2)*1000)/10) '% of PyKS Matches and ' num2str(round(performance(1)*1000)/10) '% of PyKS non matches'])

        disp('Results if training would be done with PyKs stitched')
        [ParameterkernelsPyKS,Performance] = CreateNaiveBayes(Tbl,PyKSLabel(:),Priors);
        [Fakelabel, Fakeposterior,performance] = ApplyNaiveBayes(Tbl,ParameterkernelsPyKS,PyKSLabel(:),Priors);
        disp(['Correctly labelled ' num2str(round(performance(2)*1000)/10) '% of PyKS Matches and ' num2str(round(performance(1)*1000)/10) '% of PyKS non matches'])
    else
        [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,[0 1],Priors);
    end
else
    [label, posterior, cost] = predict(BestMdl,Tbl);
end
MatchProbability = reshape(posterior(:,2),size(Predictors,1),size(Predictors,2));
label = (MatchProbability>=param.ProbabilityThreshold);% | (MatchProbability>0.05 & RankScoreAll==1 & SigMask==1);
% label = reshape(label,nclus,nclus);
[r, c] = find(triu(label)==1); %Find matches across 2 days
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs,'rows');
Pairs(Pairs(:,1)==Pairs(:,2),:) = [];

figure; imagesc(label)
colormap(flipud(gray))
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
title('Identified matches')
makepretty
saveas(gcf,fullfile(param.SaveDir,'IdentifiedMatches.fig'))
saveas(gcf,fullfile(param.SaveDir,'IdentifiedMatches.bmp'))

disp(['Extracting final pair of units took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
