%function  [UniqueID, Prob] = MatchUnitsAlgorithm(clusinfo,AllRawPaths)
%% Match units on neurophysiological evidence
% Input:
% - clusinfo (this is phy output, see also prepareinfo/spikes toolbox)
% - AllRawPaths: cell struct with paths for individual recording sessions

% Output:
% - UniqueID (Units with large overlap in QMs are likely the same unit, and
% will share a UniqueID)
% - Prob: Probability of all units to be the same as every other unit

% Matching occurs on:
% - Wavform Similarity
% - Location differences
% - Spatial decay differences
% - Waveform duration differences
% - peak time differences

% Testing the match:
% - cross-correlation finger prints --> units that are the same are likely
% to correlate in a similar way with other units

% Contributions: 
% Enny van Beest (Dec 2022)


%% Parameters
sampleamount = 500; % Nr. waveforms to include
spikeWidth = 83; % in sample space (time)
TakeChannelRadius = 50; %in micron around max channel
TakeChannelRadiusWaveform = 25; %Slightly more selective; to calculate average waveform shape
binsz = 0.01; % Binsize in time for the cross-correlation fingerprint. We recommend ~2-10ms time windows
RedoExtraction = 0; % Raw waveform and parameter extraction
Scores2Include = {'WavformSimilarity','LocationCombined','waveformdurationDiff','spatialdecayDiff','PeakTimeDiff'};%}
MakeOwnNaiveBayes = 1; % if 0, use standard matlab version, which assumes normal distributions
ApplyExistingBayesModel = 0; %If 1, look if a Bayes model already exists for this mouse and applies that

%% Extract all clusters
AllClusterIDs = clusinfo.cluster_id;
nses = length(AllQMsPaths);
OriginalClusID = AllClusterIDs;
UniqueID = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
GoodRecSesID = clusinfo.RecSesID(Good_Idx);
SessionSwitch = find(GoodRecSesID==2,1,'first');
% Define day stucture
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
SameSesMat = arrayfun(@(X) cell2mat(arrayfun(@(Y) GoodRecSesID(X)==GoodRecSesID(Y),1:nclus,'Uni',0)),1:nclus,'Uni',0);
SameSesMat = cat(1,SameSesMat{:});

%% Load raw waveforms and extract waveform parameters
halfWidth = floor(spikeWidth / 2);
dataTypeNBytes = numel(typecast(cast(0, 'uint16'), 'uint8'));
nChannels = param.nChannels;

% Initialize
ProjectedLocation = nan(2,nclus,2);
ProjectedLocationPerTP = nan(2,nclus,spikeWidth,2);
ProjectedWaveform = nan(spikeWidth,nclus,2); % Just take waveform on maximal channel
PeakTime = nan(nclus,2); % Peak time first versus second half
MaxChannel = nan(nclus,2); % Max channel first versus second half
waveformduration = nan(nclus,2); % Waveformduration first versus second half
spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second half
%Calculate how many channels are likely to be included
fakechannel = [channelpos(1,1) nanmean(channelpos(:,2))];
ChanIdx = find(cell2mat(arrayfun(@(Y) norm(fakechannel-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
% MultiDimMatrix = nan(spikeWidth,length(ChanIdx),nclus,2); % number time points (max), number of channels (max?), per cluster and cross-validated

% Take geographically close channels (within 50 microns!), not just index!
timercounter = tic;
fprintf(1,'Extracting raw waveforms. Progress: %3d%%',0)
pathparts = strsplit(AllRawPaths{GoodRecSesID(1)},'\');
rawdatapath = dir(fullfile('\\',pathparts{1:end-1}));
if isempty(rawdatapath)
    rawdatapath = dir(fullfile(pathparts{1:end-1}));
end
if exist(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(AllClusterIDs(Good_Idx(1))) '_RawSpikes.mat'])) && RedoExtraction
    delete(fullfile(rawdatapath(1).folder,'RawWaveforms','*'))
end

Currentlyloaded = 0;
for uid = 1:nclus
    pathparts = strsplit(AllRawPaths{GoodRecSesID(uid)},'\');
    rawdatapath = dir(fullfile('\\',pathparts{1:end-1}));
    if isempty(rawdatapath)
        rawdatapath = dir(fullfile(pathparts{1:end-1}));
    end

    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    if exist(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(AllClusterIDs(Good_Idx(uid))) '_RawSpikes.mat'])) && ~RedoExtraction
        load(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(AllClusterIDs(Good_Idx(uid))) '_RawSpikes.mat']))
    else
        if ~(GoodRecSesID(uid) == Currentlyloaded) % Only load new memmap if not already loaded
            % Map the data
            spikeFile = dir(AllDecompPaths{GoodRecSesID(uid)});
            try %hacky way of figuring out if sync channel present or not
                n_samples = spikeFile.bytes / (param.nChannels * dataTypeNBytes);
                nChannels = param.nChannels - 1; % Last channel is sync, ignore for now
                ap_data = memmapfile(AllDecompPaths{GoodRecSesID(uid)}, 'Format', {'int16', [param.nChannels, n_samples], 'data'});
            catch
                nChannels = param.nChannels - 1;
                n_samples = spikeFile.bytes / (nChannels * dataTypeNBytes);
                ap_data = memmapfile(AllDecompPaths{GoodRecSesID(uid)}, 'Format', {'int16', [nChannels, n_samples], 'data'});
            end
            memMapData = ap_data.Data.data;
            Currentlyloaded = GoodRecSesID(uid);
        end

        % Spike samples
        idx1=(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid)) & sp.RecSes == GoodRecSesID(uid)).*round(sp.sample_rate));  % Spike times in samples;

        %Extract raw waveforms on the fly - % Unit uid
        try
            spikeIndicestmp = sort(datasample(idx1,sampleamount,'replace',false));
        catch ME
            spikeIndicestmp = idx1;
        end
        spikeMap = nan(spikeWidth,nChannels,sampleamount);
        for iSpike = 1:length(spikeIndicestmp)
            thisSpikeIdx = int32(spikeIndicestmp(iSpike));
            if thisSpikeIdx > halfWidth && (thisSpikeIdx + halfWidth) * dataTypeNBytes < spikeFile.bytes % check that it's not out of bounds
                tmp = smoothdata(double(memMapData(1:nChannels,thisSpikeIdx-halfWidth:thisSpikeIdx+halfWidth)),2,'gaussian',5);
                tmp = (tmp - mean(tmp(:,1:10),2))';
                tmp(:,end+1:nChannels) = nan(size(tmp,1),nChannels-size(tmp,2));
                % Subtract first 10 samples to level spikes
                spikeMap(:,:,iSpike) = tmp;
            end
        end
    end
    % Extract unit parameters -
    nwavs = sum(sum(~isnan(nanmean(spikeMap,2)),1) == spikeWidth); % Actual number of waves
    % Cross-validate: first versus second half of session
    for cv = 1:2
        if cv==1
            wavidx = floor(1:nwavs/2);
        else
            wavidx = floor(nwavs/2+1:nwavs);
        end
        % Find maximum channels:
        [~,MaxChannel(uid,cv)] = nanmax(nanmax(abs(nanmean(spikeMap(:,:,wavidx),3)),[],1));

        % Extract channel positions that are relevant and extract mean location
        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,cv),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);

        % Mean location:
        mu = sum(repmat(nanmax(abs(nanmean(spikeMap(:,ChanIdx,wavidx),3)),[],1),size(Locs,2),1).*Locs',2)./sum(repmat(nanmax(abs(nanmean(spikeMap(:,ChanIdx,wavidx),3)),[],1),size(Locs,2),1),2);
        ProjectedLocation(:,uid,cv)=mu;

        %     % Mean waveform - first extract the 'weight' for each channel, based on
        %     % how close they are to the projected location (closer = better)
        Distance2MaxChan = sqrt(nansum(abs(Locs-channelpos(MaxChannel(uid,cv),:)).^2,2));
        % Difference in amplitude from maximum amplitude
        spdctmp = nanmax(abs(nanmean(spikeMap(:,MaxChannel(uid,cv),wavidx),3)),[],1)-nanmax(abs(nanmean(spikeMap(:,ChanIdx,wavidx),3)),[],1);
        % Spatial decay (average oer micron)
        spatialdecay(uid,cv) = nanmean(spdctmp./Distance2MaxChan');

        % Use this waveform - weighted average across channels:
        Distance2MaxProj = sqrt(nansum(abs(Locs-ProjectedLocation(:,uid,cv)').^2,2));
        weight = (TakeChannelRadius-Distance2MaxProj)./TakeChannelRadius;
        ProjectedWaveform(:,uid,cv) = nansum(nanmean(spikeMap(:,ChanIdx,wavidx),3).*repmat(weight,1,size(spikeMap,1))',2)./sum(weight);

        % How long is the waveform
        wvdurtmp = find(abs(ProjectedWaveform(:,uid,cv))>abs(nanmean(ProjectedWaveform(1:20,uid,cv)))+2.5*nanstd(ProjectedWaveform(1:20,uid,cv))); % More than 2. std from baseline
        wvdurtmp(find(logical([0; diff(wvdurtmp)>5]),1,'last'):end) = []; %Remove samples that are far away from other samples; noise
        waveformduration(uid,cv) = length(wvdurtmp);

        % Peak Time
        [~,PeakTime(uid,cv)] = nanmax(abs(ProjectedWaveform(wvdurtmp(1):wvdurtmp(end),uid,cv)));
        PeakTime(uid,cv) = PeakTime(uid,cv)+wvdurtmp(1)-1;
        % Mean Location per individual time point:
        ProjectedLocationPerTP(:,uid,wvdurtmp,cv) = cell2mat(arrayfun(@(tp) sum(repmat(abs(nanmean(spikeMap(tp,ChanIdx,wavidx),3)),size(Locs,2),1).*Locs',2)./sum(repmat(abs(nanmean(spikeMap(tp,ChanIdx,wavidx),3)),size(Locs,2),1),2),wvdurtmp','Uni',0));

        % Save spikes for these channels
        %         MultiDimMatrix(wvdurtmp,1:length(ChanIdx),uid,cv) = nanmean(spikeMap(wvdurtmp,ChanIdx,wavidx),3);

    end
    if ~exist(fullfile(rawdatapath(1).folder,'RawWaveforms'))
        mkdir(fullfile(rawdatapath(1).folder,'RawWaveforms'))
    end
    save(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(AllClusterIDs(Good_Idx(uid))) '_RawSpikes.mat']),'spikeMap')
end
fprintf('\n')
disp(['Extracting raw waveforms and parameters took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])

%% Metrics
% PeakTime = nan(nclus,2); % Peak time first versus second half
% MaxChannel = nan(nclus,2); % Max channel first versus second half
% waveformduration = nan(nclus,2); % Waveformduration first versus second half
% spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second half
disp('Computing Metric similarity between pairs of units...')
timercounter = tic;
PeakTimeDiff = arrayfun(@(uid) cell2mat(arrayfun(@(uid2) abs(PeakTime(uid,1)-PeakTime(uid2,2)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
PeakTimeDiff = cat(1,PeakTimeDiff{:});
PeakTimeDiff = 1-((PeakTimeDiff-nanmin(PeakTimeDiff(:)))./(nanmax(PeakTimeDiff(:))-nanmin(PeakTimeDiff(:))));

waveformdurationDiff = arrayfun(@(uid) cell2mat(arrayfun(@(uid2) abs(waveformduration(uid,1)-waveformduration(uid2,2)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
waveformdurationDiff = cat(1,waveformdurationDiff{:});
waveformdurationDiff = 1-((waveformdurationDiff-nanmin(waveformdurationDiff(:)))./(nanmax(waveformdurationDiff(:))-nanmin(waveformdurationDiff(:))));

spatialdecayDiff = arrayfun(@(uid) cell2mat(arrayfun(@(uid2) abs(spatialdecay(uid,1)-spatialdecay(uid2,2)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
spatialdecayDiff = cat(1,spatialdecayDiff{:});
spatialdecayDiff = 1-((spatialdecayDiff-nanmin(spatialdecayDiff(:)))./(nanmax(spatialdecayDiff(:))-nanmin(spatialdecayDiff(:))));

disp(['Calculating other metrics took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])


%% Waveform similarity
disp('Computing waveform similarity between pairs of units...')
timercounter = tic;
RawWVMSE = arrayfun(@(uid) cell2mat(arrayfun(@(uid2) nanmean((ProjectedWaveform(:,uid,1)-ProjectedWaveform(:,uid2,2)).^2)./nanmean(nanmean(cat(2,ProjectedWaveform(:,uid,1),ProjectedWaveform(:,uid2,2)),2).^2),1:nclus,'Uni',0)),1:nclus,'Uni',0);
RawWVMSE = cat(1,RawWVMSE{:});
WVCorr = arrayfun(@(uid) cell2mat(arrayfun(@(uid2) corr(ProjectedWaveform(:,uid,1),ProjectedWaveform(:,uid2,2)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
WVCorr = cat(1,WVCorr{:});

% Make WVCorr a normal distribution
WVCorr = erfinv(WVCorr);
WVCorr = (WVCorr-nanmin(WVCorr(:)))./(nanmax(WVCorr(:))-nanmin(WVCorr(:)));

% Normalize distribution
RawWVMSElog = log10(RawWVMSE);
WavformSimilarity = (RawWVMSElog-nanmin(RawWVMSElog(:)))./(nanmax(RawWVMSElog(:))-nanmin(RawWVMSElog(:)));
WavformSimilarity = 1-WavformSimilarity;

disp(['Calculating waveform similarity took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])

figure('name','Waveform similarity measures')
subplot(1,3,1)
h=imagesc(WVCorr);
title('Waveform correlations')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
colorbar
makepretty

subplot(1,3,2)
h=imagesc(WavformSimilarity);
title('Waveform mean squared errors')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
colorbar
makepretty

subplot(1,3,3)
h=imagesc((WVCorr+WavformSimilarity)/2);
title('Average Waveform scores')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
colorbar
makepretty
%% Location differences between pairs of units:
flag=0;
while flag<2
    figure('name','Projection locations all units')
    scatter(channelpos(:,1),channelpos(:,2),10,[0 0 0],'filled')
    hold on
    scatter(nanmean(ProjectedLocation(1,:,:),3),nanmean(ProjectedLocation(2,:,:),3),10,GoodRecSesID)
    colormap redblue
    makepretty
    xlabel('XPos (um)')
    ylabel('YPos (um)')
    xlim([min(channelpos(:,1))-50 max(channelpos(:,1))+50])
    drawnow

    disp('Computing location distances between pairs of units...')
    timercounter = tic;
    LocDist = arrayfun(@(X) cell2mat(arrayfun(@(Y) pdist(cat(2,ProjectedLocation(:,X,1),ProjectedLocation(:,Y,2))'),1:nclus,'Uni',0)),1:nclus,'Uni',0);
    LocDist = cat(1,LocDist{:}); % Normal difference
    % Difference in distance at different time points
    LocDistSign = arrayfun(@(uid) arrayfun(@(uid2)  cell2mat(arrayfun(@(X) pdist(cat(2,squeeze(ProjectedLocationPerTP(:,uid,X,1)),squeeze(ProjectedLocationPerTP(:,uid2,X,2)))'),1:spikeWidth,'Uni',0)),1:nclus,'Uni',0),1:nclus,'Uni',0);
    LocDistSign = cat(1,LocDistSign{:});
    % Average error
    MSELoc = cell2mat(cellfun(@(X) nanmean(abs(X)),LocDistSign,'Uni',0));
    % Minimal error
    LocalMinError = cell2mat(cellfun(@(X) nanmin(X),LocDistSign,'Uni',0));

    % Normalize each of them from 0 to 1, 1 being the 'best'
    LocDist = 1-((LocDist-nanmin(LocDist(:)))./(nanmax(LocDist(:))-nanmin(LocDist(:))));
    MSELoc = 1-((MSELoc-nanmin(MSELoc(:)))./(nanmax(MSELoc(:))-nanmin(MSELoc(:))));
    LocalMinError = 1-((LocalMinError-nanmin(LocalMinError(:)))./(nanmax(LocalMinError(:))-nanmin(LocalMinError(:))));

    figure('name','Distance Measures')
    subplot(2,3,1)
    h=imagesc(LocDist);
    title('LocationDistance')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(2,3,2)
    h=imagesc(MSELoc);
    title('average trajectory error')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(2,3,3)
    h=imagesc(LocalMinError);
    title('minimal trajectory error')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(2,3,4)
    h=imagesc((LocDist+MSELoc)/2);
    title('LocationDistance and mean error')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(2,3,6)
    h=imagesc((MSELoc+LocalMinError)/2);
    title('Mean and local error')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(2,3,6)
    h=imagesc((LocDist+LocalMinError)/2);
    title('LocationDistance and local error')
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    LocationCombined = (LocDist+LocalMinError+MSELoc)/3;
    disp(['Extracting projected location took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])

    %% Part 1 --> Sum Total scores components - Have a high threshold to have the 'real matches'
    figure('name','Total Score components');
    for sid = 1:length(Scores2Include)
        eval(['tmp = ' Scores2Include{sid} ';'])
        subplot(round(sqrt(length(Scores2Include))),ceil(sqrt(length(Scores2Include))),sid)
        h=imagesc(tmp,[0 1]);
        title(Scores2Include{sid})
        xlabel('Unit Y')
        ylabel('Unit Z')
        hold on
        line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
        line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
        colormap(flipud(gray))
        colorbar
        makepretty
    end

    %% Calculate total score
    figure;
    if length(Scores2Include)>1
        for scid=1:length(Scores2Include)


            ScoresTmp = Scores2Include;
            ScoresTmp(scid)=[];

            TotalScore = zeros(nclus,nclus);
            for scid2=1:length(ScoresTmp)
                eval(['TotalScore=TotalScore+' ScoresTmp{scid2} ';'])
            end
            base = length(ScoresTmp)-1;

            TotalScoreAcrossDays = TotalScore;
            TotalScoreAcrossDays(X==Y)=nan;


            subplot(2,length(Scores2Include)+1,scid)
            h=imagesc(triu(TotalScore,1),[0 base+1]);
            title(['without ' Scores2Include{scid}])
            xlabel('Unit Y')
            ylabel('Unit Z')
            hold on
            line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
            line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
            colormap(flipud(gray))
            colorbar
            makepretty

            % Thresholds
            ThrsOpt = quantile(diag(TotalScore),0.1); %Can be lenient, will select bst ones later
            subplot(2,length(Scores2Include)+1,scid+(length(Scores2Include)+1))

            imagesc(triu(TotalScore>ThrsOpt,1))
            hold on
            title(['Thresholding at ' num2str(ThrsOpt)])
            xlabel('Unit Y')
            ylabel('Unit Z')
            hold on
            line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
            line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
            colormap(flipud(gray))
            colorbar
            %     axis square
            makepretty
        end
    end
    TotalScore = zeros(nclus,nclus);
    Predictors = zeros(nclus,nclus,0);
    for scid2=1:length(Scores2Include)
        eval(['TotalScore=TotalScore+' Scores2Include{scid2} ';'])
        Predictors = cat(3,Predictors,eval(Scores2Include{scid2}));
    end
    subplot(2,length(Scores2Include)+1,length(Scores2Include)+1)
    h=imagesc(triu(TotalScore,1),[0 length(Scores2Include)]);
    title(['Total Score'])
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    makepretty

    % Make initial threshold --> to be optimized
    ThrsOpt = quantile(diag(TotalScore),0.05); %Can be lenient, will select bst ones later
    subplot(2,length(Scores2Include)+1,2*(length(Scores2Include)+1))
    imagesc(triu(TotalScore>ThrsOpt,1))
    hold on
    title(['Thresholding at ' num2str(ThrsOpt)])
    xlabel('Unit Y')
    ylabel('Unit Z')
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    colormap(flipud(gray))
    colorbar
    % axis square
    makepretty
    % Find all pairs
    % first factor authentication: score above threshold
    ThrsScore = ThrsOpt;
    % Take into account:
    label = triu(TotalScore>ThrsOpt,1) & triu(~SameSesMat,1);
    [uid,uid2] = find(label);
    Pairs = cat(2,uid,uid2);
    Pairs = sortrows(Pairs);
    Pairs=unique(Pairs,'rows');
    Pairs(Pairs(:,1)==Pairs(:,2),:)=[];

    %For chronic
    PairsPyKS = [];
    % for uid = 1:nclus
    %     pairstmp = find(AllClusterIDs(Good_Idx)==AllClusterIDs(Good_Idx(uid)))';
    %     if length(pairstmp)>1
    %         PairsPyKS = cat(1,PairsPyKS,pairstmp);
    %     end
    % end
    % PairsPyKS=unique(PairsPyKS,'rows');
    %
    % [Int,A,B] = intersect(Pairs,PairsPyKS,'rows');
    % PercDetected = size(Int,1)./size(PairsPyKS,1).*100;
    % disp(['Detected ' num2str(PercDetected) '% of PyKS matched units'])
    %
    % PercOver = (size(Pairs,1)-size(Int,1))./size(PairsPyKS,1)*100;
    % disp(['Detected ' num2str(PercOver) '% more units than just PyKS matched units'])

    % % interesting: Not detected
    % NotB = 1:size(PairsPyKS,1);
    % NotB(B) = [];
    % OnlyDetectedByPyKS = PairsPyKS(NotB,:);
    %
    % % Too much detected
    % NotA = 1:size(Pairs,1);
    % NotA(A) = [];
    % NotdetectedByPyKS = Pairs(NotA,:);

    %% Functional score for optimization: compute Fingerprint for the matched units - based on Célian Bimbard's noise-correlation finger print method but applied to across session correlations
    % Not every recording day will have the same units. Therefore we will
    % correlate each unit's activity with average activity across different
    % depths    
    timevec = floor(min(sp.st)):binsz:ceil(max(sp.st));
    edges = floor(min(sp.st))-binsz/2:binsz:ceil(max(sp.st))+binsz/2;
    disp('Calculate activity correlations')

    % Use a bunch of units with high total scores as reference population
    Pairs = sortrows(Pairs')';
    % Now sort based on how well they match
    [PairScore,sortid] = sort(cell2mat(arrayfun(@(X) TotalScore(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
    Pairs = Pairs(sortid,:);
    % Only use every 'unit' once --> take the highest scoring matches
    [val,id1,id2]=unique(Pairs(:,1),'stable');
    Pairs = Pairs(id1,:);
    [val,id1,id2]=unique(Pairs(:,2),'stable');
    Pairs = Pairs(id1,:);

    % Correlation on first day
    srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,1),'UniformOutput',0);
    srMatches = cat(1,srMatches{:});
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),1:SessionSwitch-1,'UniformOutput',0);
    srAll = cat(1,srAll{:});
    SessionCorrelation_Pair1 = corr(srMatches(:,1:size(srMatches,2)./2)',srAll(:,1:size(srMatches,2)./2)');
    for pid = 1:size(Pairs,1)
        SessionCorrelation_Pair1(pid,Pairs(pid,1)) = nan;
    end


    figure('name','Cross-correlation Fingerprints')
    subplot(1,3,1)
    imagesc(SessionCorrelation_Pair1')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title('Day 1')
    makepretty

    % Correlation on first day
    srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,2),'UniformOutput',0);
    srMatches = cat(1,srMatches{:});
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),SessionSwitch:nclus,'UniformOutput',0);
    srAll = cat(1,srAll{:});
    SessionCorrelation_Pair2 = corr(srMatches(:,size(srMatches,2)./2+1:end)',srAll(:,size(srMatches,2)./2+1:end)');
    for pid = 1:size(Pairs,1)
        SessionCorrelation_Pair2(pid,Pairs(pid,2)-SessionSwitch+1) = nan;
    end

    subplot(1,3,2)
    imagesc(SessionCorrelation_Pair2')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title('Day 2')
    makepretty

    % Add both together
    SessionCorrelations = cat(2,SessionCorrelation_Pair1,SessionCorrelation_Pair2)';

    % Correlate 'fingerprints'
    FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
    FingerprintR = cat(1,FingerprintR{:});


    subplot(1,3,3)
    imagesc(FingerprintR)
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    clim([0.5 1])
    colormap(flipud(gray))
    xlabel('All units across both days')
    ylabel('All units across both days')
    title('Correlation Fingerprint')
    makepretty

    %
    SigMask = zeros(nclus,nclus);
    RankScoreAll = nan(size(SigMask));
    for pid=1:nclus
        for pid2 = 1:nclus
            if pid2<SessionSwitch
                tmp1 = FingerprintR(pid,1:SessionSwitch-1);
                addthis=0;
            else
                tmp1 = FingerprintR(pid,SessionSwitch:end);
                addthis = SessionSwitch-1;
            end
            [val,ranktmp] = sort(tmp1,'descend');

            tmp1(pid2-addthis)=[];

            if FingerprintR(pid,pid2)>nanmean(tmp1)+2*nanstd(tmp1)
                SigMask(pid,pid2)=1;
            end
            RankScoreAll(pid,pid2) = find(ranktmp==pid2-addthis);

            if 0%(any(ismember(Pairs(:,1),pid)) & any(ismember(Pairs(:,2),pid2))) || RankScoreAll(pid,pid2)==1
                tmpfig = figure;
                subplot(2,1,1); plot(SessionCorrelations(pid,:)); hold on; plot(SessionCorrelations(pid2,:))
                xlabel('Unit')
                ylabel('Cross-correlation')
                legend('Day1','Day2')
                makepretty
                title(['r=' num2str(FingerprintR(pid,pid2)) ', rank = ' num2str(RankScoreAll(pid,pid2))])


                subplot(2,1,2)
                tmp1 = FingerprintR(1:end,pid);
                tmp1(pid)=[];
                tmp2 = FingerprintR(pid,1:end)';
                tmp2(pid)=[];

                histogram([tmp1;tmp2]);
                hold on
                line([FingerprintR(pid,pid2) FingerprintR(pid,pid2)],get(gca,'ylim'),'color',[1 0 0])
                makepretty

                pause(1)
                delete(tmpfig)
            end

        end
    end

    figure;
    scatter(TotalScore(:),FingerprintR(:),14,RankScoreAll(:),'filled')
    colormap(cat(1,[0 0 0],winter))
    xlabel('TotalScore')
    ylabel('Cross-correlation fingerprint')
    makepretty

    %% three ways to define candidate scores
    % Total score larger than threshold
    CandidatePairs = TotalScore>ThrsOpt & RankScoreAll==1;
    CandidatePairs(tril(true(size(CandidatePairs))))=0;
    figure('name','Potential Matches')
    imagesc(CandidatePairs)
    colormap(flipud(gray))
    xlim([SessionSwitch nclus])
    ylim([1 SessionSwitch-1])
    xlabel('Units day 1')
    ylabel('Units day 2')
    title('Potential Matches')
    makepretty

    %% Calculate median drift on this population (between days)
    drift = nanmedian(cell2mat(arrayfun(@(uid) (nanmean(ProjectedLocation(:,Pairs(uid,1),:),3)-nanmean(ProjectedLocation(:,Pairs(uid,2),:),3)),1:size(Pairs,1),'Uni',0)),2);
    disp(['Median drift calculated: X=' num2str(drift(1)) ', Y=' num2str(drift(2))])
    if flag
        break
    end
    ProjectedLocation(1,GoodRecSesID==2,:)=ProjectedLocation(1,GoodRecSesID==2,:)+drift(1);
    ProjectedLocation(2,GoodRecSesID==2,:)=ProjectedLocation(2,GoodRecSesID==2,:)+drift(2);
    ProjectedLocationPerTP(1,GoodRecSesID==2,:,:) = ProjectedLocationPerTP(1,GoodRecSesID==2,:,:) + drift(1);
    ProjectedLocationPerTP(2,GoodRecSesID==2,:,:) = ProjectedLocationPerTP(2,GoodRecSesID==2,:,:) + drift(2);

    flag = flag+1;
    close all

end
%% Prepare naive bayes - inspect probability distributions
figure('name','Parameter Scores');
stepsize = 0.01;
Edges = [0:stepsize:1];
for scid=1:length(Scores2Include)
    eval(['ScoresTmp = ' Scores2Include{scid} ';'])
    ScoresTmp(tril(true(size(ScoresTmp))))=nan;
    subplot(length(Scores2Include),2,(scid-1)*2+1)
    histogram(ScoresTmp(~CandidatePairs),Edges)
    if scid==1
        title('Candidate non-Matches')
    end
    ylabel(Scores2Include{scid})
    makepretty
    
    subplot(length(Scores2Include),2,scid*2)
    histogram(ScoresTmp(CandidatePairs),Edges)
    if scid==1
        title('Candidate Matches')
    end
    makepretty
end

%
figure('name','Projected Location Distance to [0 0]')
Dist2Tip = sqrt(nansum(ProjectedLocation.^2,1));
% Dist2TipMatrix = nan(size(CandidatePairs));

Dist2TipMatrix = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,Dist2Tip(X),Dist2Tip(Y)),1:nclus,'Uni',0)),1:nclus,'Uni',0);
Dist2TipMatrix = cat(3,Dist2TipMatrix{:});
Dist2TipMatrix = reshape(Dist2TipMatrix,2,[]);
subplot(1,2,1)
[N,C] = hist3(Dist2TipMatrix(:,~CandidatePairs(:))');
imagesc(C{1},C{2},N)
colormap(flipud(gray))
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Non-matches')
makepretty

subplot(1,2,2)
[N,C] = hist3(Dist2TipMatrix(:,CandidatePairs(:))');
imagesc(C{1},C{2},N)
colormap(flipud(gray))
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Matches')
makepretty


% Waveform duration
figure('name','WaveDur')
waveformdurationMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,waveformduration(X),waveformduration(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
waveformdurationMat = cat(3,waveformdurationMat{:});
subplot(1,2,1)
[N,C] = hist3(waveformdurationMat(:,~CandidatePairs(:))');
imagesc(C{1},C{2},N)
colormap(flipud(gray))
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Non-matches')
makepretty

subplot(1,2,2)
[N,C] = hist3(waveformdurationMat(:,CandidatePairs(:))');
imagesc(C{1},C{2},N)
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Matches')
makepretty


% SpatialDecaySlope
figure('name','Spatial Decay Slope')
SpatDecMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,spatialdecay(X),spatialdecay(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
SpatDecMat = cat(3,SpatDecMat{:});
subplot(1,2,1)
[N,C] = hist3(SpatDecMat(:,~CandidatePairs(:))');
imagesc(C{1},C{2},N)
colormap(flipud(gray))
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Non-matches')
makepretty

subplot(1,2,2)
[N,C] = hist3(SpatDecMat(:,CandidatePairs(:))');
imagesc(C{1},C{2},N)
colormap(flipud(gray))

xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Candidate Matches')
makepretty
%% Naive bayes classifier
% Usually this means there's no variance in the match distribution
% (which in a way is great). Create some small variance
flag = 0;
npairs = 0;
MinLoss=1;
npairslatest = 0;
maxrun = 1; % Probably we don't want to keep optimizing, as this can be a bit circular (?)
runid=0;
BestMdl = [];
while flag<2 && runid<maxrun
    flag = 0;
    runid=runid+1
    if ApplyExistingBayesModel && exist(fullfile(SaveDir,MiceOpt{midx},'UnitMatchModel.mat'))
        load(fullfile(SaveDir,MiceOpt{midx},'UnitMatchModel.mat'),'BestMdl')
        % Apply naive bays classifier
        Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters

        if isfield(BestMdl,'Parameterkernels')
            [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels);
        else
            [label, posterior, cost] = predict(BestMdl,Tbl);
        end
    else
        Tbl = array2table(reshape(Predictors(Pairs(:,1),Pairs(:,2),:),[],size(Predictors(Pairs(:,1),Pairs(:,2),:),3)),'VariableNames',Scores2Include); %All parameters
        % Use Rank as 'correct' label
        label = reshape(CandidatePairs(Pairs(:,1),Pairs(:,2)),[],1);
        if MakeOwnNaiveBayes
            % Work in progress
            [Parameterkernels,Performance] = CreateNaiveBayes(Tbl,label);
            if Performance>0.9
                flag = flag+1;
            end
            % Apply naive bays classifier
            Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
            [label, posterior] = ApplyNaiveBayes(Tbl,Parameterkernels);
            BestMdl.Parameterkernels = Parameterkernels;

        else % This uses matlab package. Warning: normal distributions assumed?
            try
                Mdl = fitcnb(Tbl,label);
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
            Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
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


        end
    end

    label = reshape(label,size(Predictors,1),size(Predictors,2));
    [r, c] = find(triu(label,1)==1 & triu(~SameSesMat,1)); %Find matches

    Pairs = cat(2,r,c);
    Pairs = sortrows(Pairs);
    Pairs=unique(Pairs,'rows');
    Pairs(Pairs(:,1)==Pairs(:,2),:)=[];
    MatchProbability = reshape(posterior(:,2),size(Predictors,1),size(Predictors,2));
    figure; imagesc(label)

    % Functional score for optimization: compute Fingerprint for the matched units - based on Célian Bimbard's noise-correlation finger print method but applied to across session correlations
    % Not every recording day will have the same units. Therefore we will
    % correlate each unit's activity with average activity across different
    % depths
    disp('Recalculate activity correlations')

    % Use a bunch of units with high total scores as reference population
    [PairScore,sortid] = sort(cell2mat(arrayfun(@(X) MatchProbability(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
    Pairs = Pairs(sortid,:);
    % Only use every 'unit' once
    [val,id1,id2]=unique(Pairs(:,1),'stable');
    Pairs = Pairs(id1,:);
    [val,id1,id2]=unique(Pairs(:,2),'stable');
    Pairs = Pairs(id1,:);
    npairs = size(Pairs,1);
    if npairs==npairslatest
        flag = flag+1;
    end
    npairslatest=npairs;

    disp(['Npairs = ' num2str(npairs)])

    if flag==2
        disp('This will not get any better... quite while ahead')
        break
    end

    % Correlation on first day
    srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,1),'UniformOutput',0);
    srMatches = cat(1,srMatches{:});
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),1:SessionSwitch-1,'UniformOutput',0);
    srAll = cat(1,srAll{:});
    SessionCorrelation_Pair1 = corr(srMatches(:,1:size(srMatches,2)./2)',srAll(:,1:size(srMatches,2)./2)');
    for pid = 1:size(Pairs,1)
        SessionCorrelation_Pair1(pid,Pairs(pid,1)) = nan;
    end


    figure('name','Cross-correlation Fingerprints')
    subplot(1,3,1)
    imagesc(SessionCorrelation_Pair1')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title('Day 1')
    makepretty

    % Correlation on second day
    srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,2),'UniformOutput',0);
    srMatches = cat(1,srMatches{:});
    srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),SessionSwitch:nclus,'UniformOutput',0);
    srAll = cat(1,srAll{:});
    SessionCorrelation_Pair2 = corr(srMatches(:,size(srMatches,2)./2+1:end)',srAll(:,size(srMatches,2)./2+1:end)');
    for pid = 1:size(Pairs,1)
        SessionCorrelation_Pair2(pid,Pairs(pid,2)-SessionSwitch+1) = nan;
    end

    subplot(1,3,2)
    imagesc(SessionCorrelation_Pair2')
    colormap(flipud(gray))
    xlabel('Candidate Units to be matched')
    ylabel('All units')
    title('Day 2')
    makepretty

    % Add both together
    SessionCorrelations = cat(2,SessionCorrelation_Pair1,SessionCorrelation_Pair2)';

    % Correlate 'fingerprints'
    FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
    FingerprintR = cat(1,FingerprintR{:});

    subplot(1,3,3)
    imagesc(FingerprintR)
    hold on
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    clim([0.5 1])
    colormap(flipud(gray))
    xlabel('All units across both days')
    ylabel('All units across both days')
    title('Correlation Fingerprint')
    makepretty

    %
    SigMask = zeros(nclus,nclus);
    RankScoreAll = nan(size(SigMask));
    for pid=1:nclus
        for pid2 = 1:nclus
            if pid2<SessionSwitch
                tmp1 = FingerprintR(pid,1:SessionSwitch-1);
                addthis=0;
            else
                tmp1 = FingerprintR(pid,SessionSwitch:end);
                addthis = SessionSwitch-1;
            end
            [val,ranktmp] = sort(tmp1,'descend');

            tmp1(pid2-addthis)=[];

            if FingerprintR(pid,pid2)>nanmean(tmp1)+2*nanstd(tmp1)
                SigMask(pid,pid2)=1;
            end
            RankScoreAll(pid,pid2) = find(ranktmp==pid2-addthis);
        end
    end
    %

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

    % Total score larger than threshold
    CandidatePairs = label==1 & RankScoreAll==1;
    CandidatePairs(tril(true(size(CandidatePairs))))=0;

end

%% Extract final pairs:
disp('Extracting final pairs of units...')
Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
if isfield(BestMdl,'Parameterkernels')
    [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels);
else
    [label, posterior, cost] = predict(BestMdl,Tbl);
end
MatchProbability = reshape(posterior(:,2),size(Predictors,1),size(Predictors,2));
label = reshape(label,nclus,nclus);
% label = MatchProbability>0.95; % We want to be very confident
[r, c] = find(triu(label,1)==1 & triu(~SameSesMat,1)); %Find matches
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs=unique(Pairs,'rows');
Pairs(Pairs(:,1)==Pairs(:,2),:)=[];
figure; imagesc(label)
colormap(flipud(gray))
xlabel('Unit_i')
ylabel('Unit_j')
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
title('Identified matches')
makepretty
%%
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
    line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
    line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
    title(['p>' num2str(takethisprob(pid))])
end
%% inspect probability distributions
figure('name','Parameter Scores');
Edges = [0:0.01:1];
for scid=1:length(Scores2Include)
    eval(['ScoresTmp = ' Scores2Include{scid} ';'])
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

%
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
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(Dist2TipMatrix(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Matches')

% Waveform duration
figure('name','WaveDur')
waveformdurationMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,waveformduration(X),waveformduration(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
waveformdurationMat = cat(3,waveformdurationMat{:});
subplot(1,2,1)
[N,C] = hist3(waveformdurationMat(:,~label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(waveformdurationMat(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Matches')

% SpatialDecaySlope
figure('name','Spatial Decay Slope')
SpatDecMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,spatialdecay(X),spatialdecay(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
SpatDecMat = cat(3,SpatDecMat{:});
subplot(1,2,1)
[N,C] = hist3(SpatDecMat(:,~label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Non-matches')

subplot(1,2,2)
[N,C] = hist3(SpatDecMat(:,label(:))');
imagesc(N)
colormap(flipud(gray))
makepretty
xlabel('Unit 1')
ylabel('Unit 2')
zlabel('Counts')
title('Identified Matches')
%%
disp('Recalculate activity correlations')

% Use a bunch of units with high total scores as reference population
[PairScore,sortid] = sort(cell2mat(arrayfun(@(X) MatchProbability(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
Pairs = Pairs(sortid,:);
% Only use every 'unit' once
[val,id1,id2]=unique(Pairs(:,1),'stable');
Pairs = Pairs(id1,:);
[val,id1,id2]=unique(Pairs(:,2),'stable');
Pairs = Pairs(id1,:);

% Correlation on first day
srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,1),'UniformOutput',0);
srMatches = cat(1,srMatches{:});
srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),1:SessionSwitch-1,'UniformOutput',0);
srAll = cat(1,srAll{:});
SessionCorrelation_Pair1 = corr(srMatches(:,1:size(srMatches,2)./2)',srAll(:,1:size(srMatches,2)./2)');
for pid = 1:size(Pairs,1)
    SessionCorrelation_Pair1(pid,Pairs(pid,1)) = nan;
end


figure('name','Cross-correlation Fingerprints')
subplot(1,3,1)
imagesc(SessionCorrelation_Pair1')
colormap(flipud(gray))
xlabel('Candidate Units to be matched')
ylabel('All units')
title('Day 1')
makepretty

% Correlation on first day
srMatches = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),Pairs(:,2),'UniformOutput',0);
srMatches = cat(1,srMatches{:});
srAll = arrayfun(@(X) histcounts(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(X)) & sp.RecSes == GoodRecSesID(X)),edges),SessionSwitch:nclus,'UniformOutput',0);
srAll = cat(1,srAll{:});
SessionCorrelation_Pair2 = corr(srMatches(:,size(srMatches,2)./2+1:end)',srAll(:,size(srMatches,2)./2+1:end)');
for pid = 1:size(Pairs,1)
    SessionCorrelation_Pair2(pid,Pairs(pid,2)-SessionSwitch+1) = nan;
end

subplot(1,3,2)
imagesc(SessionCorrelation_Pair2')
colormap(flipud(gray))
xlabel('Candidate Units to be matched')
ylabel('All units')
title('Day 2')
makepretty

% Add both together
SessionCorrelations = cat(2,SessionCorrelation_Pair1,SessionCorrelation_Pair2)';

% Correlate 'fingerprints'
FingerprintR = arrayfun(@(X) cell2mat(arrayfun(@(Y) corr(SessionCorrelations(X,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))',SessionCorrelations(Y,~isnan(SessionCorrelations(X,:))&~isnan(SessionCorrelations(Y,:)))'),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
FingerprintR = cat(1,FingerprintR{:});

subplot(1,3,3)
imagesc(FingerprintR)
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
clim([0.5 1])
colormap(flipud(gray))
xlabel('All units across both days')
ylabel('All units across both days')
title('Correlation Fingerprint')
makepretty

%
SigMask = zeros(nclus,nclus);
RankScoreAll = nan(size(SigMask));
for pid=1:nclus
    for pid2 = 1:nclus
        if pid2<SessionSwitch
            tmp1 = FingerprintR(pid,1:SessionSwitch-1);
            addthis=0;
        else
            tmp1 = FingerprintR(pid,SessionSwitch:end);
            addthis = SessionSwitch-1;
        end
        [val,ranktmp] = sort(tmp1,'descend');

        tmp1(pid2-addthis)=[];

        if FingerprintR(pid,pid2)>nanmean(tmp1)+2*nanstd(tmp1)
            SigMask(pid,pid2)=1;
        end
        RankScoreAll(pid,pid2) = find(ranktmp==pid2-addthis);
    end
end

%%
figure;
subplot(1,3,1)
imagesc(RankScoreAll==1)
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
title('Rankscore <= 2')
makepretty

subplot(1,3,2)
imagesc(label==1)
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
title('Match label == 1')
makepretty

subplot(1,3,3)
imagesc(label==1 & RankScoreAll<=2)
hold on
line([SessionSwitch SessionSwitch],get(gca,'ylim'),'color',[1 0 0])
line(get(gca,'xlim'),[SessionSwitch SessionSwitch],'color',[1 0 0])
colormap(flipud(gray))
title('Rank<=2 and Match label == 1')
makepretty


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

%% Extract final pairs:
[r, c] = find(triu(label,1)==1); %Find matches
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs=unique(Pairs,'rows');
Pairs(Pairs(:,1)==Pairs(:,2),:)=[];

%% TotalScore Pair versus no pair
SelfScore = MatchProbability(logical(eye(size(MatchProbability))));
scorematches = nan(size(Pairs,1),1); %First being TotalScore, second being TemplateMatch
scoreNoMatches = MatchProbability;
scoreNoMatches(logical(eye(size(MatchProbability))))=nan;

for id = 1:size(Pairs,1)
    scorematches(id,1) = MatchProbability(Pairs(id,1),Pairs(id,2));
    scoreNoMatches(Pairs(id,1),Pairs(id,2),:)=nan;
    scoreNoMatches(Pairs(id,2),Pairs(id,1),:)=nan;
end
ThrsScore = min(MatchProbability(label==1));
figure;
subplot(1,2,1)
histogram(scoreNoMatches(:),[0:0.01:1]); hold on
title('Non Matches')
xlabel('Match Probability')
ylabel('Nr Pairs')
makepretty
subplot(1,2,2)
histogram(SelfScore(:),[0:0.01:1]); hold on
histogram(scorematches(:),[0:0.01:1]);
line([ThrsScore ThrsScore],get(gca,'ylim'),'color',[1 0 0])

% histogram(scorematches(:,1),[0:0.02:6])
xlabel('Matching Probability')
ylabel('Nr Pairs')
legend('Self Score','Matches','Threshold','Location','best')
makepretty

save(fullfile(SaveDir,MiceOpt{midx},'MatchingScores.mat'),'BestMdl','SessionSwitch','GoodRecSesID','AllClusterIDs','Good_Idx','WavformSimilarity','WVCorr','LocationCombined','waveformdurationDiff','PeakTimeDiff','spatialdecayDiff','TotalScore','label','MatchProbability')
save(fullfile(SaveDir,MiceOpt{midx},'UnitMatchModel.mat'),'BestMdl')
%% ISI violations (for over splits matching)
ISIViolationsScore = nan(1,size(Pairs,1));
fprintf(1,'Computing functional properties similarity. Progress: %3d%%',0)
for pairid= 1:size(Pairs,1)
    if GoodRecSesID(Pairs(pairid,1)) == GoodRecSesID(Pairs(pairid,2))
        idx1 = sp.spikeTemplates == AllClusterIDs(Good_Idx(Pairs(pairid,1)))&sp.RecSes == GoodRecSesID(Pairs(pairid,1));
        idx2 = sp.spikeTemplates == AllClusterIDs(Good_Idx(Pairs(pairid,2)))&sp.RecSes == GoodRecSesID(Pairs(pairid,2));
        DifScore = diff(sort([sp.st(idx1); sp.st(idx2)]));
        ISIViolationsScore(pairid) = sum(DifScore.*1000<1.5)./length(DifScore);
        fprintf(1,'\b\b\b\b%3.0f%%',pairid/size(Pairs,1)*100)

    end
end
fprintf('\n')
disp(['Removing ' num2str(sum(ISIViolationsScore>0.05)) ' matched oversplits, as merging them will violate ISI >5% of the time'])
Pairs(ISIViolationsScore>0.05,:)=[];

%% Average in 3rd dimension (halfs of a session)
ProjectedWaveform = nanmean(ProjectedWaveform,3); %Average over first and second half of session
ProjectedLocation = nanmean(ProjectedLocation,3);
ProjectedLocationPerTP = nanmean(ProjectedLocationPerTP,4);
%% Figures
% Pairs = Pairs(any(ismember(Pairs,[8,68,47,106]),2),:);
AllClusterIDs(Good_Idx(Pairs))
for pairid=1:size(Pairs,1)
    uid = Pairs(pairid,1);
    uid2 = Pairs(pairid,2);

    pathparts = strsplit(AllRawPaths{GoodRecSesID(uid)},'\');
    rawdatapath = dir(fullfile('\\',pathparts{1:end-1}));
    if isempty(rawdatapath)
        rawdatapath = dir(fullfile(pathparts{1:end-1}));
    end

    % Load raw data
    SM1=load(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(num2str(AllClusterIDs(Good_Idx(uid)))) '_RawSpikes.mat']));
    SM1 = SM1.spikeMap; %Average across these channels

    pathparts = strsplit(AllRawPaths{GoodRecSesID(uid2)},'\');
    rawdatapath = dir(fullfile('\\',pathparts{1:end-1}));
    if isempty(rawdatapath)
        rawdatapath = dir(fullfile(pathparts{1:end-1}));
    end

    SM2=load(fullfile(rawdatapath(1).folder,'RawWaveforms',['Unit' num2str(num2str(AllClusterIDs(Good_Idx(uid2)))) '_RawSpikes.mat']));
    SM2 = SM2.spikeMap; %Average across these channels

    tmpfig = figure;
    subplot(3,3,[1,4])
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,cv),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);
    for id = 1:length(Locs)
        plot(Locs(id,1)*5+[1:size(SM1,1)],Locs(id,2)*10+nanmean(SM1(:,ChanIdx(id),:),3),'b-','LineWidth',1)
        hold on
    end
    plot(ProjectedLocation(1,uid)*5+[1:size(SM1,1)],ProjectedLocation(2,uid)*10+ProjectedWaveform(:,uid),'b--','LineWidth',2)

    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid2,cv),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);
    for id = 1:length(Locs)
        plot(Locs(id,1)*5+[1:size(SM2,1)],Locs(id,2)*10+nanmean(SM2(:,ChanIdx(id),:),3),'r-','LineWidth',1)
        hold on
    end
    plot(ProjectedLocation(1,uid2)*5+[1:size(SM1,1)],ProjectedLocation(2,uid2)*10+ProjectedWaveform(:,uid2),'r--','LineWidth',2)

    makepretty
    set(gca,'xticklabel',arrayfun(@(X) num2str(X./5),cellfun(@(X) str2num(X),get(gca,'xticklabel')),'UniformOutput',0))
    set(gca,'yticklabel',arrayfun(@(X) num2str(X./10),cellfun(@(X) str2num(X),get(gca,'yticklabel')),'UniformOutput',0))
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    title(['unit' num2str(AllClusterIDs(Good_Idx(uid))) ' versus unit' num2str(AllClusterIDs(Good_Idx(uid2))) ', ' 'RecordingDay ' num2str(GoodRecSesID(uid)) ' versus ' num2str(GoodRecSesID(uid2)) ', Probability=' num2str(round(MatchProbability(uid,uid2).*100)) '%'])

    subplot(3,3,[2])

    hold on
    h(1) = plot(squeeze(ProjectedLocationPerTP(1,uid,:)),squeeze(ProjectedLocationPerTP(2,uid,:)),'b-');
    scatter(squeeze(ProjectedLocationPerTP(1,uid,:)),squeeze(ProjectedLocationPerTP(2,uid,:)),30,1:spikeWidth,'filled')
    colormap(flipud(gray))

    h(2) = plot(squeeze(ProjectedLocationPerTP(1,uid2,:)),squeeze(ProjectedLocationPerTP(2,uid2,:)),'r-');
    scatter(squeeze(ProjectedLocationPerTP(1,uid2,:)),squeeze(ProjectedLocationPerTP(2,uid2,:)),30,1:spikeWidth,'filled')
    colormap(flipud(gray))
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    ydif = diff(get(gca,'ylim'));
    xdif = diff(get(gca,'xlim'));
    stretch = (ydif-xdif)./2;
    set(gca,'xlim',[min(get(gca,'xlim')) - stretch, max(get(gca,'xlim')) + stretch])
    %     legend([h(1),h(2)],{['Unit ' num2str(uid)],['Unit ' num2str(uid2)]})
    hc= colorbar;
    hc.Label.String = 'timesample';
    makepretty
    title(['Location score: ' num2str(round(LocationCombined(uid,uid2)*100)./100)])

    subplot(3,3,5)
    plot(channelpos(:,1),channelpos(:,2),'k.')
    hold on
    h(1)=plot(channelpos(MaxChannel(uid),1),channelpos(MaxChannel(uid),2),'b.','MarkerSize',15);
    h(2) = plot(channelpos(MaxChannel(uid2),1),channelpos(MaxChannel(uid2),2),'r.','MarkerSize',15);
    xlabel('X position')
    ylabel('um from tip')
    makepretty
    title(['Chan ' num2str(MaxChannel(uid)) ' versus ' num2str(MaxChannel(uid2))])

    subplot(3,3,3)
    hold on
    SM1 = squeeze(nanmean(SM1(:,MaxChannel(uid),:),2));
    SM2 = squeeze(nanmean(SM2(:,MaxChannel(uid2),:),2));
    h(1)=plot(nanmean(SM1(:,1:2:end),2),'b-');
    h(2)=plot(nanmean(SM1(:,2:2:end),2),'b--');
    h(3)=plot(nanmean(SM2(:,1:2:end),2),'r-');
    h(4)=plot(nanmean(SM2(:,2:2:end),2),'r--');
    makepretty
    title(['Waveform Similarity=' num2str(round(WavformSimilarity(uid,uid2)*100)./100) ', dur=' num2str(round(waveformdurationDiff(uid,uid2)*100)./100) ', peak=' ...
        num2str(round(PeakTimeDiff(uid,uid2)*100)./100) ', decay='  num2str(round(spatialdecayDiff(uid,uid2)*100)./100)])


    % Scatter spikes of each unit
    subplot(3,3,6)
    idx1=find(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid)) & sp.RecSes == GoodRecSesID(uid));
    scatter(sp.st(idx1)./60,sp.spikeAmps(idx1),4,[0 0 1],'filled')
    hold on
    idx2=find(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid2)) &  sp.RecSes == GoodRecSesID(uid2));
    scatter(sp.st(idx2)./60,-sp.spikeAmps(idx2),4,[1 0 0],'filled')
    xlabel('Time (min)')
    ylabel('Abs(Amplitude)')
    title(['Amplitude distribution'])
    xlims = get(gca,'xlim');
    ylims = max(abs(get(gca,'ylim')));
    % Other axis
    [h1,edges,binsz]=histcounts(sp.spikeAmps(idx1));
    %Normalize between 0 and 1
    h1 = ((h1-nanmin(h1))./(nanmax(h1)-nanmin(h1)))*10+xlims(2)+10;
    plot(h1,edges(1:end-1),'b-');
    [h2,edges,binsz]=histcounts(sp.spikeAmps(idx2));
    %Normalize between 0 and 1
    h2 = ((h2-nanmin(h2))./(nanmax(h2)-nanmin(h2)))*10+xlims(2)+10;
    plot(h2,-edges(1:end-1),'r-');
    ylabel('Amplitude')
    ylim([-ylims ylims])

    makepretty


    % compute ACG
    [ccg, ~] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
        ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function
    ACG = ccg(:, 1, 1);
    [ccg, ~] = CCGBz([double(sp.st(idx2)); double(sp.st(idx2))], [ones(size(sp.st(idx2), 1), 1); ...
        ones(size(sp.st(idx2), 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function
    ACG2 = ccg(:, 1, 1);
    [ccg, ~] = CCGBz([double(sp.st([idx1;idx2])); double(sp.st([idx1;idx2]))], [ones(size(sp.st([idx1;idx2]), 1), 1); ...
        ones(size(sp.st([idx1;idx2]), 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function

    subplot(3,3,7); plot(ACG,'b');
    hold on
    plot(ACG2,'r')
    title(['AutoCorrelogram'])
    makepretty
    subplot(3,3,8)

    if exist('NatImgCorr','var')
        if GoodRecSesID(uid)==1 % Recording day 1
            tmp1 = squeeze(D0(OriginalClusID(Good_Idx(uid))+1,:,:));
        else % Recordingday 2
            tmp1 = squeeze(D1(OriginalClusID(Good_Idx(uid))+1,:,:));
        end
        if GoodRecSesID(uid2)==1 % Recording day 1
            tmp2 = squeeze(D0(OriginalClusID(Good_Idx(uid2))+1,:,:));
        else % Recordingday 2
            tmp2 = squeeze(D1(OriginalClusID(Good_Idx(uid2))+1,:,:));
        end

        plot(nanmean(tmp1,1),'b-');
        hold on
        plot(nanmean(tmp2,1),'r-');
        xlabel('Stimulus')
        ylabel('NrSpks')
        makepretty


        if AllClusterIDs(Good_Idx(uid))  == AllClusterIDs(Good_Idx(uid2))
            if ismember(Good_Idx(uid),Good_ClusUnTracked)
                title(['Visual: Untracked, r=' num2str(round(NatImgCorr(pairid,pairid)*100)/100)])
            elseif ismember(Good_Idx(uid),Good_ClusTracked)
                title(['Visual: Tracked, r=' num2str(round(NatImgCorr(pairid,pairid)*100)/100)])
            else
                title(['Visual: Unknown, r=' num2str(round(NatImgCorr(pairid,pairid)*100)/100)])
            end
        else
            title(['Visual: Unknown, r=' num2str(round(NatImgCorr(pairid,pairid)*100)/100)])
        end
    else
        isitot = diff(sort([sp.st(idx1); sp.st(idx2)]));
        histogram(isitot,'FaceColor',[0 0 0])
        hold on
        line([1.5/1000 1.5/1000],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
        title([num2str(round(sum(isitot*1000<1.5)./length(isitot)*1000)/10) '% ISI violations']); %The higher the worse (subtract this percentage from the Total score)
        xlabel('ISI (ms)')
        ylabel('Nr. Spikes')
        makepretty
    end

    subplot(3,3,9)

    plot(SessionCorrelations(uid,:),'b-'); hold on; plot(SessionCorrelations(uid2,:),'r-')
    xlabel('Unit')
    ylabel('Cross-correlation')
    title(['Fingerprint r=' num2str(round(FingerprintR(uid,uid2)*100)/100) ', rank=' num2str(RankScoreAll(uid,uid2))])
    makepretty

    disp(['UniqueID ' num2str(AllClusterIDs(Good_Idx(uid))) ' vs ' num2str(AllClusterIDs(Good_Idx(uid2)))])
    disp(['Peakchan ' num2str(MaxChannel(uid)) ' versus ' num2str(MaxChannel(uid2))])
    disp(['RecordingDay ' num2str(GoodRecSesID(uid)) ' versus ' num2str(GoodRecSesID(uid2))])

    drawnow
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    saveas(gcf,fullfile(SaveDir,MiceOpt{midx},[num2str(round(MatchProbability(uid,uid2).*100)) 'ClusID' num2str(AllClusterIDs(Good_Idx(uid))) 'vs' num2str(AllClusterIDs(Good_Idx(uid2))) '.fig']))
    saveas(gcf,fullfile(SaveDir,MiceOpt{midx},[num2str(round(MatchProbability(uid,uid2).*100)) 'ClusID' num2str(AllClusterIDs(Good_Idx(uid))) 'vs' num2str(AllClusterIDs(Good_Idx(uid2))) '.bmp']))

    close(tmpfig)
end
%% Unused bits and pieces
if 0
    % look for natural images data
    % AL data from Kush:
    D0 = readNPY('H:\Anna_TMP\image_analysis\responses\day_0\template_responses.npy');
    D1 = readNPY('H:\Anna_TMP\image_analysis\responses\day_1\template_responses.npy');
    % matrix: number of spikes from 0 to 0.7seconds after stimulus onset: n_templates X n_reps X n_Images
    % Get rid of the units not currently looked at (We have only shank 0 here)
    D0(1:end-length(unique(AllClusterIDs)),:,:)=[];
    D1(1:end-length(unique(AllClusterIDs)),:,:)=[];
    NatImgCorr = nan(nclus,nclus);
    nrep = size(D0,2);
    for uid=1:nclus
        uid
        if GoodRecSesID(uid)==1 % Recording day 1
            tmp1 = squeeze(D0(OriginalClusID(Good_Idx(uid))+1,:,:));
        else % Recordingday 2
            tmp1 = squeeze(D1(OriginalClusID(Good_Idx(uid))+1,:,:));
        end
        parfor uid2 = uid:nclus
            if GoodRecSesID(uid2)==1 % Recording day 1
                tmp2 = squeeze(D0(OriginalClusID(Good_Idx(uid2))+1,:,:));
            else % Recordingday 2
                tmp2 = squeeze(D1(OriginalClusID(Good_Idx(uid2))+1,:,:));
            end
            %
            %         figure; subplot(2,2,1); imagesc(tmp1); title(['Day ' num2str(GoodRecSesID(uid)) ', Unit ' num2str(OriginalClusID(Good_Idx(uid)))])
            %         colormap gray
            %         colorbar; xlabel('Condition'); ylabel('Repeat')
            %         hold on; subplot(2,2,2); imagesc(tmp2); title(['Day ' num2str(GoodRecSesID(uid2)) ', Unit ' num2str(OriginalClusID(Good_Idx(uid2)))])
            %         colormap gray
            %         colorbar; xlabel('Condition'); ylabel('Repeat')
            %         subplot(2,2,3); hold on

            % Is the unit's response predictable?
            tmpcor = nan(1,nrep);
            for cv = 1:nrep
                % define training and test
                trainidx = circshift(1:nrep,-(cv-1));
                testidx = trainidx(1);
                trainidx(1)=[];

                % Define response:
                train = nanmean(tmp1(trainidx,:),1);
                test = tmp2(testidx,:);

                % Between error
                tmpcor(1,cv) = corr(test',train');
                %             scatter(train,test,'filled')

                %             plot(train);

            end
            NatImgCorr(uid,uid2) = nanmean(nanmean(tmpcor,2));

            %         xlabel('train')
            %         ylabel('test')
            %         title(['Average Correlation ' num2str(round(NatImgCorr(uid,uid2)*100)/100)])
            %        lims = max(cat(1,get(gca,'xlim'),get(gca,'ylim')),[],1);
            %         set(gca,'xlim',lims,'ylim',lims)
            %         makepretty

        end
    end
    % Mirror these
    for uid2 = 1:nclus
        for uid=uid2+1:nclus
            NatImgCorr(uid,uid2)=NatImgCorr(uid2,uid);
        end
    end
    NatImgCorr = arrayfun(@(Y) arrayfun(@(X) NatImgCorr(Pairs(X,1),Pairs(Y,2)),1:size(Pairs,1),'UniformOutput',0),1:size(Pairs,1),'UniformOutput',0)
    NatImgCorr = cell2mat(cat(1,NatImgCorr{:}));

    % Kush's verdict:
    Good_ClusTracked = readNPY('H:\Anna_TMP\image_analysis\cluster_ids\good_clusters_tracked.npy'); % this is an index, 0 indexed so plus 1
    Good_ClusUnTracked = readNPY('H:\Anna_TMP\image_analysis\cluster_ids\good_clusters_untracked.npy') % this is an index, 0 indexed so plus 1

    Good_ClusTracked(Good_ClusTracked>max(AllClusterIDs))=[]; %
    Good_ClusUnTracked(Good_ClusUnTracked>max(AllClusterIDs)) = [];

    NotIncluded = [];
    TSGoodGr = nan(1,length(Good_ClusTracked));
    for uid = 1:length(Good_ClusTracked)
        idx = find(AllClusterIDs(Good_Idx) == Good_ClusTracked(uid));
        if length(idx)==2
            TSGoodGr(uid) = TotalScore(idx(1),idx(2));
        else
            NotIncluded = [NotIncluded  Good_ClusTracked(uid)];
        end
    end
    TSBadGr = nan(1,length(Good_ClusUnTracked));
    for uid = 1:length(Good_ClusUnTracked)
        idx = find(AllClusterIDs(Good_Idx) == Good_ClusUnTracked(uid));
        if length(idx)==2
            TSBadGr(uid) = TotalScore(idx(1),idx(2));
        else
            NotIncluded = [NotIncluded  Good_ClusUnTracked(uid)];
        end
    end
    figure; histogram(TSBadGr,[5:0.05:6]); hold on; histogram(TSGoodGr,[5:0.05:6])
    xlabel('Total Score')
    ylabel('Nr. Matches')
    legend({'Not tracked','Tracked'})
    makepretty
    %Tracked?
    Tracked = zeros(nclus,nclus);
    for uid=1:nclus
        parfor uid2 = 1:nclus
            if uid==uid2
                Tracked(uid,uid2)=0;
            elseif AllClusterIDs(Good_Idx(uid))  == AllClusterIDs(Good_Idx(uid2))
                if ismember(Good_Idx(uid),Good_ClusUnTracked) ||  ismember(Good_Idx(uid2),Good_ClusUnTracked)
                    Tracked(uid,uid2)=-1;
                elseif ismember(Good_Idx(uid),Good_ClusTracked)||  ismember(Good_Idx(uid2),Good_ClusTracked)
                    Tracked(uid,uid2)=1;
                else
                    Tracked(uid,uid2) = 0.5;
                end
            end
        end
    end
end

%% Cross-correlation
if 0
MaxCorrelation = nan(nclus,nclus);
EstTimeshift = nan(nclus,nclus);
EstDrift = nan(nclus,nclus,2);

% Create grid-space in time and space domain
[SpaceMat,TimeMat] = meshgrid(cat(2,-[size(MultiDimMatrix,2)-1:-1:0],1:size(MultiDimMatrix,2)-1),cat(2,-[spikeWidth-1:-1:0],[1:spikeWidth-1]));
for uid = 1:nclus
    parfor uid2 = 1:nclus
        tmp1 = MultiDimMatrix(:,:,uid,1);
        tmp2 = MultiDimMatrix(:,:,uid2,2);
        %Convert nan to 0
        tmp1(isnan(tmp1))=0;
        tmp2(isnan(tmp2))=0;

        %         2D cross correlation
        c = xcorr2(tmp1,tmp2);

        % Find maximum correlation
        [MaxCorrelation(uid,uid2),indx] = max(c(:)); 

        % Index in time domain
        EstTimeshift(uid,uid2)=TimeMat(indx); % Time shift index
        EstDrift(uid,uid2)=SpaceMat(indx);

        % Shift tmp2 by those pixels/time points
        tmp2 = circshift(tmp2,EstTimeshift(uid,uid2),1);
        tmp2 = circshift(tmp2,EstDrift(uid,uid2),2);

        MaxCorrelation(uid,uid2) = corr(tmp1(tmp1(:)~=0&tmp2(:)~=0),tmp2(tmp1(:)~=0&tmp2(:)~=0));

    end
end
end
