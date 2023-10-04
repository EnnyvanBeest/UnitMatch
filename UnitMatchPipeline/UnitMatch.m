function  [UniqueIDConversion, MatchTable, WaveformInfo, param] = UnitMatch(clusinfo,param)
%% Match units on neurophysiological evidence
% Input:
% - clusinfo (this is a struct that contains per unit the following information):
% * cluster_id (e.g. kilosort output clus_id)
% * Good_ID: ones for units that should be included in the analysis
% * RecSesID: Recording Session ID
% * Coordinates: Typically 3D coordinates per unit
% * Depth: depth on probe
% * Shank: Which shank 
% * Probe: Which probe
% - param: parameters for ephys extraction
% See UMDefaultParam.m for details

% Output:
% - UniqueIDConversion (Units that are found to be a match are given the same
% UniqueID. This UniqueID can be used for further analysis
% - MatchTable: Probability, rank score and cross-correlation correlation
% (Fingerprint correlation) of all possible unit pairs
% - WaveformInfo: most important waveform information that is useful to
% keep
% - param (additional parameters used in UnitMatch)

% Matching occurs on:
% - Waveform Similarity: Correlation and errors
% - Projected location difference (Centroid): distance and direction, also
% per time point
% - Amplitude differences
% - Spatial decay (decrease in signal over space)

% Contributions:
% Enny van Beest (2022-2023)
% Célian Bimbard (2022-2023)

%% Parameters - tested on these values, but feel free to try others
GlobalUnitMatchClock = tic;        
Scores2Include = param.Scores2Include % Good to show for failure prevention
TakeChannelRadius = 150; %in micron around max channel;
maxdist = 100; % Maximum distance at which units are considered as potential matches %best tested so far 100
param.removeoversplits = 0; % Remove oversplits based on ISI violations or not?
param.MakeOwnNaiveBayes = 1; % if 0, use standard matlab version, which assumes normal distributions --> not recommended
SaveScoresAsProbability = 0; %If 1, the individual scores are converted to probabiliti
param.maxrun = 1; % This is whether you want to use Bayes' output to create a new potential candidate set to optimize the probability distributions. Probably we don't want to keep optimizing?, as this can be a bit circular (?)
param.drawmax = inf; % Maximum number of drawed matches (otherwise it takes forever!)
param.VisibleSetting = 'off'; %Do we want to see the figures being plot online?
Draw2DMatrixes = 0; % If you find this useful
param.NeighbourDist = 50; % In micron

global stepsize
stepsize = 0.01;
%% Read in from param
Allchannelpos = param.AllChannelPos;
RunPyKSChronicStitched = param.RunPyKSChronicStitched;
SaveDir = param.SaveDir;
param.nChannels = length(Allchannelpos{1})+1; %First assume there's a sync channel as well.
spikeWidth = param.spikeWidth; %83; % in sample space (time)
NewPeakLoc = floor(spikeWidth./2); % This is where all peaks will be aligned to!
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
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
    if isempty(Good_Idx)
        error('No good units in these recordings (yet asking for it). check!')
    end
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')
end
GoodRecSesID = clusinfo.RecSesID(Good_Idx);

% Define day stucture
recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
ndays = length(unique(recsesGood));
% x = repmat(GoodRecSesID,[1 numel(GoodRecSesID)]);
% SameSesMat = x == x';
% OriSessionSwitch = cell2mat(arrayfun(@(X) find(recsesAll==X,1,'first'),1:ndays,'Uni',0));
% OriSessionSwitch = [OriSessionSwitch nclus+1];
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),unique(recsesGood),'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch); nclus+1];


%% Extract raw waveforms
% This script does the actual extraction (if necessary) and saves out paths
% to NPY for individual unit data
Path4UnitNPY = ExtractAndSaveAverageWaveforms(clusinfo,param);

%% Extract parameters used in UnitMatch
[AllWVBParameters,param] = ExtractParameters(Path4UnitNPY,clusinfo,param);

%% Metrics
param = ExtractSimilarityMetrics(Scores2Include,AllWVBParameters,clusinfo,param);% All Scores2Include are pushed to the workspace

%% Naive bayes classifier
[MatchProbability,label,Tbl,BestMdl] = RunNaiveBayes(Predictors,TotalScore,Scores2Include,clusinfo,param,SortingOrder,EuclDist);

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

%% Assign same Unique ID
OriUniqueID = UniqueID; %need for plotting
[PairID1,PairID2]=meshgrid(OriginalClusterIDs(Good_Idx));
[recses1,recses2] = meshgrid(recsesAll(Good_Idx));
[PairID3,PairID4]=meshgrid(OriUniqueID(Good_Idx));

MatchTable = table(PairID1(:),PairID2(:),recses1(:),recses2(:),PairID3(:),PairID4(:),MatchProbability(:),TotalScore(:),EuclDist(:),'VariableNames',{'ID1','ID2','RecSes1','RecSes2','UID1','UID2','MatchProb','TotalScore','EucledianDistance'});
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

%% Save
UMparam = param;
UMrunTime = toc(GlobalUnitMatchClock);

save(fullfile(UMparam.SaveDir, 'UnitMatch.mat'), 'UniqueIDConversion', 'MatchTable', 'WaveformInfo', 'UMparam','UMrunTime','-v7.3')

%% From here on we're just evaluating the model:
%% If this was stitched pykilosort, we know what pykilosort thought about the matches
PyKSLabel = [];
PairsPyKS = [];
if RunPyKSChronicStitched
    for uid = 1:nclus
        pairstmp = find(OriginalClusterIDs(Good_Idx)==OriginalClusterIDs(Good_Idx(uid)))';
        if length(pairstmp)>1
            for id1 = 1:length(pairstmp)
                for id2 = 1:length(pairstmp)
                    if id1 >= id2
                        continue
                    end
                    PairsPyKS = cat(1,PairsPyKS,[pairstmp(id1),pairstmp(id2)]);
                end
            end
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


    PercDetected = sum(PyKSLabel(:) == 1 & label(:) ==1)./sum(PyKSLabel(:)==1)*100;    
    disp(['Detected ' num2str(PercDetected) '% of PyKS matched units'])

    PercOver = sum(label(:)==1 & PyKSLabel(:)==0)./sum(PyKSLabel(:)==1)*100;
    disp(['Detected ' num2str(PercOver) '% more units than just PyKS matched units'])
end

%% Check different probabilities, what does the match graph look like?
figure('name','Different Posterior probability Thresholds');
takethisprob = [0.5 0.75 0.95 0.99];
for pid = 1:4
    subplot(2,2,pid)
    h = imagesc(MatchProbability(SortingOrder,SortingOrder)>takethisprob(pid));
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

%% TotalScore Pair versus no pair
figure('name','TotalScore vs Probability');

for id = 1:2
    if id == 2
        SelfScore = MatchProbability(logical(eye(size(MatchProbability))));
        OtherScores = MatchProbability; %First being TotalScore, second being TemplateMatch
        ThrsScore = min(MatchProbability(label==1));
        Edges = [0:0.1:1];
        Vector = [0.05:0.1:1-0.05];
    else % TotalScore
        SelfScore = TotalScore(logical(eye(size(MatchProbability))));
        OtherScores = TotalScore; %First being TotalScore, second being TemplateMatch
        ThrsScore = min(TotalScore(label==1));
         Edges = [0:0.1:1];
        Vector = [0.05:0.1:1-0.05];
    end
    OtherScores(logical(eye(size(MatchProbability)))) = nan; %Get rid of diagonal
    OtherScores(EuclDist>param.NeighbourDist) = nan;%Remove units that were too far away
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

    subplot(1,2,id)
    hs = histcounts(SelfScore(:),Edges);
    hw = histcounts(WithinScores(~isnan(WithinScores)),Edges);
    ha = histcounts(AcrossScores(~isnan(AcrossScores)),Edges);

    plot(Vector,hs./sum(hs),'-','color',[0.5 0.5 0.5])
    hold on
    plot(Vector,ha./sum(ha),'-','color',[0 0.5 0])
    plot(Vector,hw./sum(hw),'-','color',[0 0 0])

%     line([ThrsScore ThrsScore],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')

    % histogram(scorematches(:,1),[0:0.02:6])
    if id == 2
        xlabel('Matching Probability')
        title('After naive Bayes')
    else
        xlabel('Total Score')
        title('Using total score')
    end
    ylabel('Proportion|Group')
    legend('Self Score',['Across C_i_j<' num2str(param.NeighbourDist)],['Within C_i_j<' num2str(param.NeighbourDist)],'Location','best')
    makepretty
end
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


return

