function  [UniqueIDConversion, MatchTable, WaveformInfo, AllSessionCorrelations, param] = UnitMatch(clusinfo,param)
%% Match units on neurophysiological evidence
% Input:
% - clusinfo (this is phy output, see also prepareinfo/spikes toolbox)
% - param: parameters for ephys extraction

% Output:
% - UniqueIDConversion (Units that are found to be a match are given the same
% UniqueID. This UniqueID can be used for further analysis
% - MatchTable: Probability, rank score and cross-correlation correlation
% (Fingerprint correlation) of all possible unit pairs
% - WaveformInfo: most important waveform information that is useful to
% keep
% - AllSessionCorrelations: Cross correlations between sessions, useful for
% functional confirmation of matches
% - param (additional parameters used in UnitMatch)

% Matching occurs on:
% - Waveform Similarity: Correlation and errors
% - Projected location difference (Centroid): distance and direction, also
% per time point
% - Amplitude differences
% - Spatial decay (decrease in signal over space)

% fine tuning the initial training set for matching:
% - cross-correlation finger prints --> units that are the same are likely
% to correlate in a similar way with other units

% Contributions:
% Enny van Beest (2022-2023)
% Célian Bimbard (2022-2023)

%% Parameters - tested on these values, but feel free to try others
Scores2Include = param.Scores2Include % Good to show for failure prevention
TakeChannelRadius = 75; %in micron around max channel
maxdist = 200; % Maximum distance at which units are considered as potential matches
param.MakeOwnNaiveBayes = 1; % if 0, use standard matlab version, which assumes normal distributions --> not recommended
SaveScoresAsProbability = 0; %If 1, the individual scores are converted to probabiliti
param.maxrun = 1; % This is whether you want to use Bayes' output to create a new potential candidate set to optimize the probability distributions. Probably we don't want to keep optimizing?, as this can be a bit circular (?)
drawmax = inf; % Maximum number of drawed matches (otherwise it takes forever!)
VisibleSetting = 'on'; %Do we want to see the figures being plot online?
Draw2DMatrixes = 0; % If you find this useful
global stepsize
stepsize = 0.01;
%% Read in from param
Allchannelpos = param.channelpos;
if ~iscell(Allchannelpos)
    Allchannelpos = {Allchannelpos};
end
RunPyKSChronicStitched = param.RunPyKSChronicStitched;
SaveDir = param.SaveDir;
% AllDecompPaths = param.AllDecompPaths;
% AllRawPaths = param.AllRawPaths;
param.nChannels = length(Allchannelpos{1})+1; %First assume there's a sync channel as well.
% sampleamount = param.sampleamount; %500; % Nr. waveforms to include
spikeWidth = param.spikeWidth; %83; % in sample space (time)
% UseBombCelRawWav = param.UseBombCelRawWav; % If Bombcell was also applied on this dataset, it's faster to read in the raw waveforms extracted by Bombcell
NewPeakLoc = floor(spikeWidth./2); % This is where all peaks will be aligned to!
% waveidx = NewPeakLoc-7:NewPeakLoc+16; % Force this analysis window
waveidx = NewPeakLoc-7:NewPeakLoc+15; % Force this analysis window So far
% best option
param.TakeChannelRadius = TakeChannelRadius;
param.waveidx = waveidx;
param.SaveScoresAsProbability = SaveScoresAsProbability;
param.NewPeakLoc = NewPeakLoc;
param.maxdist = maxdist;

if ~exist(SaveDir)
    mkdir(SaveDir)
end

%% Extract all cluster info
OriginalClusterIDs = clusinfo.cluster_id;
% nses = length(AllDecompPaths);
% OriginalClusID = AllClusterIDs; % Original cluster ID assigned by KS
UniqueID = 1:length(OriginalClusterIDs); % Initial assumption: All clusters are unique
Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
GoodRecSesID = clusinfo.RecSesID(Good_Idx);

% Define day stucture
recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
ndays = length(unique(recsesAll));
% x = repmat(GoodRecSesID,[1 numel(GoodRecSesID)]);
% SameSesMat = x == x';
% OriSessionSwitch = cell2mat(arrayfun(@(X) find(recsesAll==X,1,'first'),1:ndays,'Uni',0));
% OriSessionSwitch = [OriSessionSwitch nclus+1];
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),1:ndays,'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch) nclus+1];

%% Extract raw waveforms
% This script does the actual extraction (if necessary) and saves out paths
% to NPY for individual unit data
Path4UnitNPY = ExtractAndSaveAverageWaveforms(clusinfo,param);

%% Extract parameters used in UnitMatch
AllWVBParameters = ExtractParameters(Path4UnitNPY,clusinfo,param);

%% Metrics
ExtractSimilarityMetrics(Scores2Include,AllWVBParameters,clusinfo,param)% All Scores2Include are pushed to the workspace

%% Naive bayes classifier
[MatchProbability,label,Pairs,Tbl,BestMdl] = RunNaiveBayes(Predictors,TotalScore,Scores2Include,clusinfo,param);

%% Some evaluation:
% Units on the diagonal are matched by (Py)KS within a day. Very likely to
% be correct:
disp(['Evaluating Naive Bayes...'])
FNEst = (1-(sum(diag(MatchProbability)>param.ProbabilityThreshold)./nclus))*100;
disp([num2str(round(sum(diag(MatchProbability)>param.ProbabilityThreshold)./nclus*100)) '% of units were matched with itself'])
disp(['False negative estimate: ' num2str(round(FNEst*100)/100) '%'])
if FNEst>10
    warning('Warning, false negatives very high!')
end
lowselfscores = find(diag(MatchProbability<param.ProbabilityThreshold));

% Units off diagonal within a day are not matched by (Py)KS within a day. 
FPEst=nan(1,ndays);
for did = 1:ndays
    tmpprob = double(MatchProbability(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)>param.ProbabilityThreshold);
    tmpprob(logical(eye(size(tmpprob)))) = nan;
    FPEst(did) = sum(tmpprob(:)==1)./sum(~isnan(tmpprob(:)))*100;
    disp(['False positive estimate recording ' num2str(did) ': ' num2str(round(FPEst(did)*100)/100) '%'])
end
BestMdl.FalsePositiveEstimate = FPEst;
BestMdl.FalseNegativeEstimate = FNEst;

%% Save model and relevant information
save(fullfile(SaveDir,'MatchingScores.mat'),'BestMdl','SessionSwitch','GoodRecSesID','OriginalClusterIDs','Good_Idx','Tbl','TotalScore','label','MatchProbability','param')
save(fullfile(SaveDir,'UnitMatchModel.mat'),'BestMdl')

%% Change these parameters to probabilities of being a match
if SaveScoresAsProbability
    for pidx = 1:length(Scores2Include)
        [~,IndividualScoreprobability, ~]=ApplyNaiveBayes(Tbl(:,pidx),Parameterkernels(:,pidx,:),BestMdl.Priors);
        eval([Scores2Include{pidx} ' = reshape(IndividualScoreprobability(:,2),nclus,nclus).*100;'])
    end
end

%% Compute functional score (cross-correlation fingerprint)
if ndays<5
    drawdrosscorr = 1;
else
    drawdrosscorr = 0;
end
[r, c] = find(MatchProbability>param.ProbabilityThreshold); %Find matches
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs,'rows');
[FingerprintR,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(sessionCorrelationsAll,Pairs,Unit2Take,recsesGood,drawdrosscorr);


%% Assign same Unique ID
OriUniqueID = UniqueID; %need for plotting
[PairID1,PairID2]=meshgrid(OriginalClusterIDs(Good_Idx));
[recses1,recses2] = meshgrid(recsesAll(Good_Idx));
[PairID3,PairID4]=meshgrid(OriUniqueID(Good_Idx));

MatchTable = table(PairID1(:),PairID2(:),recses1(:),recses2(:),PairID3(:),PairID4(:),MatchProbability(:),RankScoreAll(:),FingerprintR(:),TotalScore(:),EuclDist(:),'VariableNames',{'ID1','ID2','RecSes1','RecSes2','UID1','UID2','MatchProb','RankScore','FingerprintCor','TotalScore','EucledianDistance'});
if param.AssignUniqueID
    [UniqueID, MatchTable] = AssignUniqueID(MatchTable,clusinfo,Path4UnitNPY,param);
end
% Add Scores2Include to MatchTable
for sid = 1:length(Scores2Include)
    eval(['MatchTable.' Scores2Include{sid} ' = Tbl.' Scores2Include{sid} '(:);'])
end

%% Other useful parameters to keep:
WaveformInfo.MaxChannel = AllWVBParameters.MaxChannel;
WaveformInfo.ProjectedLocation = ProjectedLocation;
WaveformInfo.ProjectedWaveform = AllWVBParameters.ProjectedWaveform;
WaveformInfo.ProjectedLocationPerTP = ProjectedLocationPerTP;

%% Save out UniqueID conversion
UniqueIDConversion.UniqueID = UniqueID;
UniqueIDConversion.OriginalClusID = OriginalClusterIDs;
UniqueIDConversion.recsesAll = recsesAll;
UniqueIDConversion.GoodID = clusinfo.Good_ID;
UniqueIDConversion.Path4UnitNPY = Path4UnitNPY;

%% From here on we're just evaluating the model:
%% If this was stitched pykilosort, we know what pykilosort thought about the matches
PyKSLabel = [];
PairsPyKS = [];
if RunPyKSChronicStitched
    for uid = 1:nclus
        pairstmp = find(OriginalClusterIDs(Good_Idx)==OriginalClusterIDs(Good_Idx(uid)))';
        if length(pairstmp)>1
            PairsPyKS = cat(1,PairsPyKS,pairstmp);
        end
    end

    PyKSLabel = false(nclus,nclus);
    for pid = 1:size(PairsPyKS,1)
        PyKSLabel(PairsPyKS(pid,1),PairsPyKS(pid,2)) = true;
        PyKSLabel(PairsPyKS(pid,2),PairsPyKS(pid,1)) = true;
    end
    PyKSLabel(logical(eye(size(PyKSLabel)))) = true;
    PairsPyKS=unique(PairsPyKS,'rows');

    figure('name','Parameter Scores');
    Edges = [0:stepsize:1];
    ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;

    for scid=1:length(Scores2Include)
        eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
        ScoresTmp = reshape(ScoresTmp,nclus,nclus);
        ScoresTmp(tril(true(size(ScoresTmp))))=nan;
        subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),scid)
        hc = histcounts(ScoresTmp(~PyKSLabel),Edges)./sum(~PyKSLabel(:));
        plot(ScoreVector,hc,'b-')
        hold on

        hc = histcounts(ScoresTmp(PyKSLabel),Edges)./sum(PyKSLabel(:));
        plot(ScoreVector,hc,'r-')


        title(Scores2Include{scid})
        makepretty
    end
    legend('Non-matches','Matches')

    [~, ~,performance] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,PyKSLabel(:),BestMdl.Priors);
    disp(['Correctly labelled ' num2str(round(performance(2)*1000)/10) '% of PyKS Matches and ' num2str(round(performance(1)*1000)/10) '% of PyKS non matches'])

    disp('Results if training would be done with PyKs stitched')
    [ParameterkernelsPyKS,~] = CreateNaiveBayes(Tbl,PyKSLabel(:),BestMdl.Priors);
    [Fakelabel, Fakeposterior,performance] = ApplyNaiveBayes(Tbl,ParameterkernelsPyKS,PyKSLabel(:),BestMdl.Priors);
    disp(['Correctly labelled ' num2str(round(performance(2)*1000)/10) '% of PyKS Matches and ' num2str(round(performance(1)*1000)/10) '% of PyKS non matches'])

    Fakeposterior = reshape(Fakeposterior(:,2),nclus,nclus);
    FNEst = (1-(sum(diag(Fakeposterior)>param.ProbabilityThreshold)./nclus))*100;
    disp([num2str(round(sum(diag(Fakeposterior)>param.ProbabilityThreshold)./nclus*100)) '% of units were matched with itself'])
    disp(['False negative estimate: ' num2str(round(FNEst*100)/100) '%'])
    if FNEst>10
        warning('Warning, false negatives very high!')
    end
    lowselfscores = find(diag(Fakeposterior<param.ProbabilityThreshold));

    % Units off diagonal within a day are not matched by (Py)KS within a day.
    FPEst=nan(1,ndays);
    for did = 1:ndays
        tmpprob = double(Fakeposterior(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)>param.ProbabilityThreshold);
        tmpprob(logical(eye(size(tmpprob)))) = nan;
        FPEst(did) = sum(tmpprob(:)==1)./sum(~isnan(tmpprob(:)))*100;
        disp(['False positive estimate recording ' num2str(did) ': ' num2str(round(FPEst(did)*100)/100) '%'])
    end
end


%% Check different probabilities, what does the match graph look like?
figure;
takethisprob = [0.5 0.75 0.95 0.99];
for pid = 1:4
    subplot(2,2,pid)
    h = imagesc(MatchProbability>takethisprob(pid));
    colormap(flipud(gray))
    makepretty
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    title(['p>' num2str(takethisprob(pid))])
end
saveas(gcf,fullfile(SaveDir,'ProbabilitiesMatches.fig'))
saveas(gcf,fullfile(SaveDir,'ProbabilitiesMatches.bmp'))

%% Compare to functional scores
figure;
subplot(1,3,1)
imagesc(RankScoreAll==1 & SigMask==1)
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
title('Rankscore == 1*')
makepretty

subplot(1,3,2)
imagesc(MatchProbability>param.ProbabilityThreshold)
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
title(['Match Probability>' num2str(param.ProbabilityThreshold)])
makepretty

subplot(1,3,3)
imagesc(MatchProbability>=param.ProbabilityThreshold | (MatchProbability>0.05 & RankScoreAll==1 & SigMask==1));
% imagesc(MatchProbability>=0.99 | (MatchProbability>=0.05 & RankScoreAll==1 & SigMask==1))
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
title('Matching probability + rank')
makepretty
saveas(gcf,fullfile(SaveDir,'RankScoreVSProbability.fig'))
saveas(gcf,fullfile(SaveDir,'RankScoreVSProbability.bmp'))

tmpf = triu(FingerprintR);
tmpm = triu(MatchProbability);
tmpr = triu(RankScoreAll);
tmpr = tmpr(tmpf~=0);
tmpm = tmpm(tmpf~=0);
tmpf = tmpf(tmpf~=0);
figure;
scatter(tmpm,tmpf,14,tmpr,'filled')
colormap(cat(1,[0 0 0],winter))
xlabel('Match Probability')
ylabel('Cross-correlation fingerprint')
makepretty
saveas(gcf,fullfile(SaveDir,'RankScoreVSProbabilityScatter.fig'))
saveas(gcf,fullfile(SaveDir,'RankScoreVSProbabilityScatter.bmp'))


if RunPyKSChronicStitched

    [Int,A,B] = intersect(Pairs,PairsPyKS,'rows');
    PercDetected = size(Int,1)./size(PairsPyKS,1).*100;
    disp(['Detected ' num2str(PercDetected) '% of PyKS matched units'])

    PercOver = (size(Pairs,1)-size(Int,1))./size(PairsPyKS,1)*100;
    disp(['Detected ' num2str(PercOver) '% more units than just PyKS matched units'])

    % interesting: Not detected
    NotB = 1:size(PairsPyKS,1);
    NotB(B) = [];
    OnlyDetectedByPyKS = PairsPyKS(NotB,:);

    % Figure out what hasn't been detected:
    % Extract individual parameter scores for these specific 'pairs'
    TmpSc = cell2mat(arrayfun(@(X) squeeze(Predictors(OnlyDetectedByPyKS(X,1),OnlyDetectedByPyKS(X,2),:)),1:size(OnlyDetectedByPyKS,1),'Uni',0));

    % p(X|Match)
    [~,minidx] = arrayfun(@(X) min(abs(X-ScoreVector)),TmpSc,'Uni',0); % Find index for observation in score vector
    minidx=cell2mat(minidx);
    % Only interested in matches:
    MatchLkh = cell2mat(arrayfun(@(Y) BestMdl.Parameterkernels(minidx(Y,:),Y,2),1:size(minidx,1),'Uni',0));

    figure;
    subplot(2,2,1)
    imagesc(TmpSc')
    colormap(flipud(gray))
    title('Scores')

    subplot(2,2,2)
    imagesc(MatchLkh)
    colormap(flipud(gray))
    title('likelihood_Match')

    TmpMP = cell2mat(arrayfun(@(X) squeeze(MatchProbability(OnlyDetectedByPyKS(X,1),OnlyDetectedByPyKS(X,2))),1:size(OnlyDetectedByPyKS,1),'Uni',0));

    figure('name','OnlyDetecedPyKS');
    for scid=1:length(Scores2Include)
        subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),scid)
        histogram(TmpSc(scid,:),Edges)
        title(Scores2Include{scid})
    end

    % Too much detected
    NotA = 1:size(Pairs,1);
    NotA(A) = [];
    NotdetectedByPyKS = Pairs(NotA,:);

end


%% TotalScore Pair versus no pair
NeighbourDist = 30; % In micron
SelfScore = MatchProbability(logical(eye(size(MatchProbability))));
OtherScores = MatchProbability; %First being TotalScore, second being TemplateMatch
OtherScores(logical(eye(size(MatchProbability)))) = nan; %Get rid of diagonal
OtherScores(EuclDist>NeighbourDist) = nan;%Remove units that were too far away
AcrossScores = nan(size(EuclDist));
WithinScores = nan(size(EuclDist));
% Divide in within and across scores
for did = 1:ndays
    WithinScores(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1) = OtherScores(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1);
    for did2 = 1:ndays
        if did==did2
            continue
        end
        AcrossScores(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did2):SessionSwitch(did2+1)-1) = OtherScores(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did2):SessionSwitch(did2+1)-1);
    end
end

ThrsScore = min(MatchProbability(label==1));
figure;
hs = histcounts(SelfScore(:),[0:0.1:1]); 
hw = histcounts(WithinScores(~isnan(WithinScores)),[0:0.1:1]);
ha = histcounts(AcrossScores(~isnan(AcrossScores)),[0:0.1:1]);

plot([0.05:0.1:1-0.05],hs./sum(hs),'-','color',[0.5 0.5 0.5])
hold on
plot([0.05:0.1:1-.005],ha./sum(ha),'g-')
plot([0.05:0.1:1-.005],hw./sum(hw),'r-')

line([ThrsScore ThrsScore],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')

% histogram(scorematches(:,1),[0:0.02:6])
xlabel('Matching Probability')
ylabel('Proportion|Group')
legend('Self Score',['Across C_i_j<' num2str(NeighbourDist)],['Within C_i_j<' num2str(NeighbourDist)],'Threshold','Location','best')
makepretty
saveas(gcf,fullfile(SaveDir,'ScoresSelfvsMatch.fig'))
saveas(gcf,fullfile(SaveDir,'ScoresSelfvsMatch.bmp'))

%% inspect probability distributions
figure('name','Parameter Scores');
Edges = [0:0.01:1];
for scid=1:length(Scores2Include)
    eval(['ScoresTmp = Tbl.' Scores2Include{scid} ';'])
    ScoresTmp = reshape(ScoresTmp,nclus,nclus);

    ScoresTmp(tril(true(size(ScoresTmp))))=nan;
    subplot(length(Scores2Include),2,(scid-1)*2+1)
    histogram(ScoresTmp(~label),Edges)
    if scid==1
        title('Identified non-Matches')
    end
    ylabel(Scores2Include{scid})
    makepretty

    subplot(length(Scores2Include),2,scid*2)
    histogram(ScoresTmp(label),Edges)
    if scid==1
        title('Identified Matches')
    end
    makepretty
end

if Draw2DMatrixes 
figure('name','Projected Location Distance to [0 0]')
Dist2Tip = sqrt(nansum(ProjectedLocation.^2,1));
% Dist2TipMatrix = nan(size(CandidatePairs));

Dist2TipMatrix = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,Dist2Tip(X),Dist2Tip(Y)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
Dist2TipMatrix = cat(3,Dist2TipMatrix{:});
Dist2TipMatrix = reshape(Dist2TipMatrix,2,[]);
subplot(1,2,1)
[N,C] = hist3(Dist2TipMatrix(:,~label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(Dist2TipMatrix(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
zlabel('Counts')
title('Identified Matches')

% Waveform duration
figure('name','WaveDur')
waveformdurationMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,AllWVBParameters.waveformduration(X),AllWVBParameters.waveformduration(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
waveformdurationMat = cat(3,waveformdurationMat{:});
subplot(1,2,1)
[N,C] = hist3(waveformdurationMat(:,~label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(waveformdurationMat(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
zlabel('Counts')
title('Identified Matches')

% SpatialDecaySlope
figure('name','Spatial Decay Slope')
SpatDecMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,AllWVBParameters.spatialdecay(X),AllWVBParameters.spatialdecay(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
SpatDecMat = cat(3,SpatDecMat{:});
subplot(1,2,1)
[N,C] = hist3(SpatDecMat(:,~label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(SpatDecMat(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit_i')
ylabel('Unit_j')
zlabel('Counts')
title('Identified Matches')
end
%% Figures
if param.MakePlotsOfPairs
    % Pairs redefined:
    uId = unique(UniqueID(Good_Idx));
    Pairs = arrayfun(@(X) find(UniqueID(Good_Idx)==X),uId,'Uni',0);
    Pairs(cellfun(@length,Pairs)==1) = [];
    for id =1:length(lowselfscores) % Add these for plotting - inspection
        Pairs{end+1} = [lowselfscores(id) lowselfscores(id)];
    end

    if RunPyKSChronicStitched
        for id = 1:size(OnlyDetectedByPyKS,1) % Add these for plotting - inspection
            Pairs{end+1} = [OnlyDetectedByPyKS(id,1) OnlyDetectedByPyKS(id,2)];
        end
    end

    if size(Pairs,2)>drawmax
        DrawPairs = randsample(1:size(Pairs,2),drawmax,'false');
    else
        DrawPairs = 1:size(Pairs,2);
    end

    PlotTheseUnits_UM(Pairs(DrawPairs),MatchTable,UniqueIDConversion,WaveformInfo,AllSessionCorrelations,param,VisibleSetting)
end
return

