function  [UniqueID, MatchTable, WaveformInfo, AllSessionCorrelations] = UnitMatch(clusinfo,param)
%% Match units on neurophysiological evidence
% Input:
% - clusinfo (this is phy output, see also prepareinfo/spikes toolbox)
% - param: parameters for ephys extraction

% Output:
% - UniqueID (Units that are found to be a match are given the same
% UniqueID. This UniqueID can be used for further analysis
% - MatchTable: Probability, rank score and cross-correlation correlation
% (Fingerprint correlation) of all possible unit pairs
% - WaveformInfo: most important waveform information that is useful to
% keep
% - AllSessionCorrelations: Cross correlations between sessions, useful for
% functional confirmation of matches

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
global stepsize
stepsize = 0.01; % Of probability distribution
MakePlotsOfPairs = 1; % Plots all pairs for inspection
Scores2Include = {'AmplitudeSim','WavformSim','LocTrajectorySim','spatialdecaySim','LocDistSim'}; %

% Scores2Include = {'AmplitudeSim','WavformMSE','WVCorr','LocTrajectorySim','spatialdecaySim','LocDistSim'}; %
TakeChannelRadius = 75; %in micron around max channel
maxdist = 200; % Maximum distance at which units are considered as potential matches
binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows
RemoveRawWavForms = 0; %Remove averaged waveforms again to save space --> Currently only two averages saved so shouldn't be a problem to keep it, normally speaking
% Scores2Include = {'WavformSimilarity','LocationCombined','spatialdecayDiff','AmplitudeDiff'};%}
MakeOwnNaiveBayes = 1; % if 0, use standard matlab version, which assumes normal distributions --> not recommended
ApplyExistingBayesModel = 0; %If 1, use probability distributions made available by us
maxrun = 1; % This is whether you want to use Bayes' output to create a new potential candidate set to optimize the probability distributions. Probably we don't want to keep optimizing?, as this can be a bit circular (?)
drawmax = inf; % Maximum number of drawed matches (otherwise it takes forever!)
VisibleSetting = 'off'; %Do we want to see the figures being plot online?

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
%% Extract all cluster info
AllClusterIDs = clusinfo.cluster_id;
% nses = length(AllDecompPaths);
% OriginalClusID = AllClusterIDs; % Original cluster ID assigned by KS
UniqueID = 1:length(AllClusterIDs); % Initial assumption: All clusters are unique
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
% Initialize
ProjectedLocation = nan(2,nclus,2);
ProjectedLocationPerTP = nan(2,nclus,spikeWidth,2);
ProjectedWaveform = nan(spikeWidth,nclus,2); % Just take waveform on maximal channel
PeakTime = nan(nclus,2); % Peak time first versus second half
MaxChannel = nan(nclus,2); % Max channel first versus second half
waveformduration = nan(nclus,2); % Waveformduration first versus second half
Amplitude = nan(nclus,2); % Maximum (weighted) amplitude, first versus second half
spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second half
WaveIdx = false(nclus,spikeWidth,2);
%Calculate how many channels are likely to be included
% fakechannel = [channelpos{1}(1,1) nanmean(channelpos{1}(:,2))];
% ChanIdx = find(cell2mat(arrayfun(@(Y) norm(fakechannel-channelpos{1}(Y,:)),1:size(channelpos{1},1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
% % MultiDimMatrix = nan(spikeWidth,length(ChanIdx),nclus,2); % number time points (max), number of channels (max?), per cluster and cross-validated

% Take geographically close channels (within 50 microns!), not just index!
timercounter = tic;
fprintf(1,'Extracting waveform information. Progress: %3d%%',0)
for uid = 1:nclus
    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    % load data
    spikeMap = readNPY(Path4UnitNPY{uid});

    % Detrending
    spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
    spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
    spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order

    channelpos = Allchannelpos{recsesGood(uid)};

    % Extract channel positions that are relevant and extract mean location
    [~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);

    % Extract unit parameters -
    % Cross-validate: first versus second half of session
    for cv = 1:2
        % Find maximum channels:
        [~,MaxChannel(uid,cv)] = nanmax(nanmax(abs(spikeMap(35:70,ChanIdx,cv)),[],1)); %Only over relevant channels, in case there's other spikes happening elsewhere simultaneously
        MaxChannel(uid,cv) = ChanIdx(MaxChannel(uid,cv));

        % Mean location:
        mu = sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1).*Locs',2)./sum(repmat(nanmax(abs(nanmean(spikeMap(:,ChanIdx,cv),3)),[],1),size(Locs,2),1),2);
        ProjectedLocation(:,uid,cv) = mu;
        % Use this waveform - weighted average across channels:
        Distance2MaxProj = sqrt(nansum(abs(Locs-ProjectedLocation(:,uid,cv)').^2,2));
        weight = (TakeChannelRadius-Distance2MaxProj)./TakeChannelRadius;
        ProjectedWaveform(:,uid,cv) = nansum(spikeMap(:,ChanIdx,cv).*repmat(weight,1,size(spikeMap,1))',2)./sum(weight);
        % Find significant timepoints
        wvdurtmp = find(abs(ProjectedWaveform(:,uid,cv))>abs(nanmean(ProjectedWaveform(1:20,uid,cv)))+2.5*nanstd(ProjectedWaveform(1:20,uid,cv))); % More than 2. std from baseline
        if isempty(wvdurtmp)
            wvdurtmp = waveidx;
        end
        wvdurtmp(~ismember(wvdurtmp,waveidx)) = []; %okay over achiever, gonna cut you off there
        if isempty(wvdurtmp)
            % May again be empty
            wvdurtmp = waveidx;
        end

        % Peak Time - to be safe take a bit of smoothing
        [~,PeakTime(uid,cv)] = nanmax(abs(smooth(ProjectedWaveform(wvdurtmp(1):wvdurtmp(end),uid,cv),2)));
        PeakTime(uid,cv) = PeakTime(uid,cv)+wvdurtmp(1)-1;
    end
    % Give each unit the best opportunity to correlate the waveform; cross
    % correlate to find the difference in peak
    [tmpcor, lags] = xcorr(ProjectedWaveform(:,uid,1),ProjectedWaveform(:,uid,2));
    [~,maxid] = max(tmpcor);
 
    % Shift accordingly
    ProjectedWaveform(:,uid,2) = circshift(ProjectedWaveform(:,uid,2),lags(maxid));
    spikeMap(:,:,2) = circshift(spikeMap(:,:,2),lags(maxid),1);
    if lags(maxid)>0
        ProjectedWaveform(1:lags(maxid),uid,2) = nan;
        spikeMap(1:lags(maxid),:,2) = nan;
    elseif lags(maxid)<0
        ProjectedWaveform(spikeWidth+lags(maxid):spikeWidth,uid,2) = nan;
        spikeMap(spikeWidth+lags(maxid):spikeWidth,:,2) = nan;
    end
    for cv = 1:2
        % Shift data so that peak is at timepoint x
        if PeakTime(uid,1)~=NewPeakLoc % Yes, take the 1st CV on purpose!
            ProjectedWaveform(:,uid,cv) = circshift(ProjectedWaveform(:,uid,cv),-(PeakTime(uid,1)-NewPeakLoc));
            spikeMap(:,:,cv) = circshift(spikeMap(:,:,cv),-(PeakTime(uid,1)-NewPeakLoc),1);
            if PeakTime(uid,1)-NewPeakLoc<0
                ProjectedWaveform(1:-(PeakTime(uid,1)-NewPeakLoc),uid,cv) = nan;
                spikeMap(1:-(PeakTime(uid,1)-NewPeakLoc),:,cv) = nan;
            else
                ProjectedWaveform(spikeWidth-(PeakTime(uid,cv)-NewPeakLoc):spikeWidth,uid,cv) = nan;
                spikeMap(spikeWidth-(PeakTime(uid,1)-NewPeakLoc):spikeWidth,:,cv) = nan;
            end
        end

        %     % Mean waveform - first extract the 'weight' for each channel, based on
        %     % how close they are to the projected location (closer = better)
        Distance2MaxChan = sqrt(nansum(abs(Locs-channelpos(MaxChannel(uid,cv),:)).^2,2));
        % Difference in amplitude from maximum amplitude
        spdctmp = (nanmax(abs(spikeMap(:,MaxChannel(uid,cv),cv)),[],1)-nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1))./nanmax(abs(spikeMap(:,MaxChannel(uid,cv),cv)),[],1);
        % Spatial decay (average oer micron)
        spatialdecay(uid,cv) = nanmean(spdctmp./Distance2MaxChan');
        Peakval = ProjectedWaveform(PeakTime(uid,cv),uid,cv);
        Amplitude(uid,cv) = Peakval;

        % Full width half maximum
        wvdurtmp = find(abs(sign(Peakval)*ProjectedWaveform(waveidx,uid,cv))>0.25*sign(Peakval)*Peakval);
        if ~isempty(wvdurtmp)
            wvdurtmp = [wvdurtmp(1):wvdurtmp(end)]+waveidx(1)-1;
            waveformduration(uid,cv) = length(wvdurtmp);
        else
            waveformduration(uid,cv) = nan;
        end

        % Mean Location per individual time point:
        ProjectedLocationPerTP(:,uid,wvdurtmp,cv) = cell2mat(arrayfun(@(tp) sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1).*Locs',2)./sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1),2),wvdurtmp,'Uni',0));
        WaveIdx(uid,wvdurtmp,cv) = 1;
        % Save spikes for these channels
        %         MultiDimMatrix(wvdurtmp,1:length(ChanIdx),uid,cv) = nanmean(spikeMap(wvdurtmp,ChanIdx,wavidx),3);

    end
    if  norm(channelpos(MaxChannel(uid,1),:)-channelpos(MaxChannel(uid,2),:))>TakeChannelRadius
        keyboard
    end
end
fprintf('\n')
disp(['Extracting raw waveforms and parameters took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
%% Metrics
% PeakTime = nan(nclus,2); % Peak time first versus second half
% MaxChannel = nan(nclus,2); % Max channel first versus second half
% waveformduration = nan(nclus,2); % Waveformduration first versus second half
% spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second half
disp('Computing Metric similarity between pairs of units...')
timercounter = tic;
x1 = repmat(PeakTime(:,1),[1 numel(PeakTime(:,1))]);
x2 = repmat(PeakTime(:,2),[1 numel(PeakTime(:,2))]);
PeakTimeSim = abs(x1 - x2');
%Normalize between 0 and 1 (values that make sense after testing, without having outliers influence this)
PeakTimeSim =1-PeakTimeSim./quantile(PeakTimeSim(:),0.99);
PeakTimeSim(PeakTimeSim<0)=0;

% can't find much better for this one
waveformTimePointSim = nan(nclus,nclus);
for uid = 1:nclus
    for uid2 = 1:nclus
        waveformTimePointSim(uid,uid2) = sum(ismember(find(WaveIdx(uid,:,1)),find(WaveIdx(uid2,:,2))))./sum((WaveIdx(uid,:,1)));
    end
end

x1 = repmat(spatialdecay(:,1),[1 numel(spatialdecay(:,1))]);
x2 = repmat(spatialdecay(:,2),[1 numel(spatialdecay(:,2))]);
spatialdecaySim = abs(x1 - x2');
% Make (more) normal
spatialdecaySim = sqrt(spatialdecaySim);
spatialdecaySim = 1-((spatialdecaySim-nanmin(spatialdecaySim(:)))./(quantile(spatialdecaySim(:),0.99)-nanmin(spatialdecaySim(:))));
spatialdecaySim(spatialdecaySim<0)=0;

% Ampitude difference
x1 = repmat(Amplitude(:,1),[1 numel(Amplitude(:,1))]);
x2 = repmat(Amplitude(:,2),[1 numel(Amplitude(:,2))]);
AmplitudeSim = abs(x1 - x2');
% Make (more) normal
AmplitudeSim = sqrt(AmplitudeSim);
AmplitudeSim = 1-((AmplitudeSim-nanmin(AmplitudeSim(:)))./(quantile(AmplitudeSim(:),.99)-nanmin(AmplitudeSim(:))));
AmplitudeSim(AmplitudeSim<0)=0;
disp(['Calculating other metrics took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])

%% Waveform similarity
disp('Computing waveform similarity between pairs of units...')
timercounter = tic;
% Normalize between 0 and 1
% 1st cv
x1 = ProjectedWaveform(waveidx,:,1);
% x1(WaveIdx(:,:,1)'==0)=nan; 
% 2nd cv
x2 = ProjectedWaveform(waveidx,:,2);
% x2(WaveIdx(:,:,2)'==0)=nan; 
% Correlate
WVCorr = corr(x1,x2,'rows','pairwise');
% Make WVCorr a normal distribution
WVCorr = atanh(WVCorr);
WVCorr = (WVCorr-nanmin(WVCorr(:)))./(nanmax(WVCorr(~isinf(WVCorr(:))))-nanmin(WVCorr(:)));

ProjectedWaveformNorm = cat(3,x1,x2);
ProjectedWaveformNorm = (ProjectedWaveformNorm-nanmin(ProjectedWaveformNorm,[],1))./(nanmax(ProjectedWaveformNorm,[],1)-nanmin(ProjectedWaveformNorm,[],1));
x1 = repmat(ProjectedWaveformNorm(:,:,1),[1 1 size(ProjectedWaveformNorm,2)]);
x2 = permute(repmat(ProjectedWaveformNorm(:,:,2),[1 1 size(ProjectedWaveformNorm,2)]),[1 3 2]);
RawWVMSE = squeeze(nanmean((x1 - x2).^2));

% sort of Normalize distribution
RawWVMSENorm = sqrt(RawWVMSE);
WavformMSE = (RawWVMSENorm-nanmin(RawWVMSENorm(:)))./(quantile(RawWVMSENorm(:),0.99)-nanmin(RawWVMSENorm(:)));
WavformMSE = 1-WavformMSE;
WavformMSE(WavformMSE<0) = 0;

disp(['Calculating waveform similarity took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])
figure('name','Waveform similarity measures')
subplot(1,3,1)
imagesc(WVCorr);
title('Waveform correlations')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
colorbar
makepretty

subplot(1,3,2)
imagesc(WavformMSE);
title('Waveform mean squared errors')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
colorbar
makepretty

WavformSim = (WVCorr+WavformMSE)/2;
subplot(1,3,3)
imagesc(WavformSim);
title('Average Waveform scores')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
colorbar
makepretty

%% Get traces for fingerprint correlations
disp('Computing neural traces...')
timercounter = tic;
Unit2Take = AllClusterIDs(Good_Idx);
srAllDays = cell(1,ndays);
for did = 1:ndays
    Unit2TakeIdxAll = find(recsesAll(Good_Idx) == did);

    % Load sp for correct day
    tmppath = strsplit(Path4UnitNPY{Unit2TakeIdxAll(1)},'RawWaveforms');
    tmp = matfile(fullfile(tmppath{1},'PreparedData.mat'));
    sp = tmp.sp;
    
    % Define edges for this dataset
    edges = floor(min(sp.st))-binsz/2:binsz:ceil(max(sp.st))+binsz/2;

    % bin data to create PSTH
    srAllDays{did} = nan(numel(Unit2TakeIdxAll),numel(edges)-1);
    for uid = 1:numel(Unit2TakeIdxAll)
        srAllDays{did}(uid,:) =  histcounts(sp.st(sp.spikeTemplates == Unit2Take(Unit2TakeIdxAll(uid))),edges);
    end
end
disp(['Calculating neural traces took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])

% clear variables to save space
clear tmp
clear sp
%% Location differences between pairs of units: - This is done twice to account for large drift between sessions
channelpos_AllCat = unique(cat(1,Allchannelpos{:}),'rows');

flag=0;
while flag<2
    figure('name','Projection locations all units')
    scatter(channelpos_AllCat(:,1),channelpos_AllCat(:,2),10,[0 0 0],'filled')
    hold on
    scatter(nanmean(ProjectedLocation(1,:,:),3),nanmean(ProjectedLocation(2,:,:),3),10,GoodRecSesID)
    colormap jet
    makepretty
    xlabel('XPos (um)')
    ylabel('YPos (um)')
    xlim([min(channelpos_AllCat(:,1))-50 max(channelpos_AllCat(:,1))+50])
    drawnow
    saveas(gcf,fullfile(SaveDir,'ProjectedLocation.fig'))
    saveas(gcf,fullfile(SaveDir,'ProjectedLocation.bmp'))

    disp('Computing location distances between pairs of units...')
    timercounter = tic;
    LocDist = sqrt((ProjectedLocation(1,:,1)'-ProjectedLocation(1,:,2)).^2 + ...
        (ProjectedLocation(2,:,1)'-ProjectedLocation(2,:,2)).^2);

    disp('Computing location distances between pairs of units, per individual time point of the waveform...')
    % Difference in distance at different time points
    x1 = repmat(squeeze(ProjectedLocationPerTP(:,:,waveidx,1)),[1 1 1 size(ProjectedLocationPerTP,2)]);
    x2 = permute(repmat(squeeze(ProjectedLocationPerTP(:,:,waveidx,2)),[1 1 1 size(ProjectedLocationPerTP,2)]),[1 4 3 2]);
    LocDistSign = squeeze(sqrt(nansum((x1-x2).^2,1)));
    w = squeeze(isnan(abs(x1(1,:,:,:)-x2(1,:,:,:))));
%     LocDistSign(w) = nan;
    % Average location + variance in location (captures the trajectory of a waveform in space)
    LocDistSim = squeeze(nanmean(LocDistSign,2)+nanvar(LocDistSign,[],2)); 
 

    % Variance in error, corrected by average error. This captures whether
    % the trajectory is consistenly separate
    MSELoc = squeeze(nanvar(LocDistSign,[],2)./nanmean(LocDistSign,2)+nanmean(LocDistSign,2));
    LocDistSign = squeeze(nanvar(LocDistSign,[],2))./sqrt(squeeze(sum(w,2))); % Variance in distance between the two traces

    disp('Computing location angle (direction) differences between pairs of units, per individual time point of the waveform...')
    % Difference in angle between two time points
    x1 = ProjectedLocationPerTP(:,:,waveidx(2):waveidx(end),:);
    x2 = ProjectedLocationPerTP(:,:,waveidx(1):waveidx(end-1),:);
    LocAngle = squeeze(atan(abs(x1(1,:,:,:)-x2(1,:,:,:))./abs(x1(2,:,:,:)-x2(2,:,:,:))));

    % Tried this out: Count the 'significant' waveform more; we're more certain about this trajectory than about parts of the waveform that are in the noise
    %     y1 = repmat(WaveIdx(:,waveidx(2:end),1),[1 1 nclus])+1;
    %     y2 = permute(repmat(WaveIdx(:,waveidx(2:end),2),[1 1 nclus]),[3 2 1])+1;
    %     x1 = repmat(LocAngle(:,:,1),[1 1 nclus]).*y1;
    %     x2 = permute(repmat(LocAngle(:,:,2),[1 1 nclus]),[3 2 1]).*y2; %
    %     w = (y1+y2);
    %     LocAngleSim = sqrt(squeeze(nansum(abs(x1-x2),2)));

    % Actually just taking the weighted sum of angles is better?
    x1 = repmat(LocAngle(:,:,1),[1 1 nclus]);
    x2 = permute(repmat(LocAngle(:,:,2),[1 1 nclus]),[3 2 1]); %
    w = ~isnan(abs(x1-x2));
    LocAngleSim = sqrt(squeeze(nansum(abs(x1-x2),2)./nansum(w,2)));
 
    % Normalize each of them from 0 to 1, 1 being the 'best'
    % If distance > maxdist micron it will never be the same unit:
    LocDistSim = 1-((LocDistSim-nanmin(LocDistSim(:)))./(maxdist-nanmin(LocDistSim(:)))); %Average difference
    LocDistSim(LocDistSim<0)=0;
    LocDist = 1-((LocDist-nanmin(LocDist(:)))./(maxdist-nanmin(LocDist(:))));
    LocDist(LocDist<0)=0;
    MSELoc = 1-((MSELoc-nanmin(MSELoc(:)))./(nanmax(MSELoc(:))-nanmin(MSELoc(:))));
    LocAngleSim = 1-((LocAngleSim-nanmin(LocAngleSim(:)))./(nanmax(LocAngleSim(:))-nanmin(LocAngleSim(:))));
    LocDistSign = 1-((LocDistSign-nanmin(LocDistSign(:)))./(quantile(LocDistSign(~isinf(LocDistSign)),0.99)-nanmin(LocDistSign(:))));
    LocDistSign(LocDistSign<0)=0;
    LocTrajectorySim = (LocAngleSim+LocDistSign)./2; % Trajectory Similarity is sum of distance + sum of angles
    LocTrajectorySim = (LocTrajectorySim-nanmin(LocTrajectorySim(:))./(nanmax(LocTrajectorySim(:))-nanmin(LocTrajectorySim(:))));
    %
    figure('name','Distance Measures')
    subplot(4,2,1)
    imagesc(LocDist);
    title('LocationDistance')
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    makepretty

    subplot(4,2,2)
    histogram(LocDist(:));
    xlabel('Score')
    makepretty

    subplot(4,2,3)
    imagesc(LocDistSim);
    title('LocationDistanceAveraged')
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    makepretty
    subplot(4,2,4)
    histogram(LocDistSim(:));
    xlabel('Score')
    makepretty

    subplot(4,2,5)
    imagesc(MSELoc);
    title('average trajectory error')
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    makepretty
    subplot(4,2,6)
    h=histogram(MSELoc(:));
    xlabel('Score')
    makepretty

    subplot(4,2,7)
    imagesc(LocTrajectorySim);
    title('Mean trajectory difference')
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    makepretty
    subplot(4,2,8)
    h=histogram(LocTrajectorySim(:));
    xlabel('Score')
    makepretty
    LocDistSim(isnan(LocDistSim))=0;
    LocationCombined = nanmean(cat(3,LocDistSim,LocTrajectorySim),3);
    disp(['Extracting projected location took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

    %% These are the parameters to include:
    figure('name','Total Score components');
    for sid = 1:length(Scores2Include)
        eval(['tmp = ' Scores2Include{sid} ';'])
        subplot(round(sqrt(length(Scores2Include))),ceil(sqrt(length(Scores2Include))),sid)
        try
            imagesc(tmp,[quantile(tmp(:),0.1) 1]);
        catch
            imagesc(tmp,[0 1]);
        end
        title(Scores2Include{sid})
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty
    end
    saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.fig'))
    saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.bmp'))

    %% Calculate total score
    disp('Computing total score...')
    timercounter = tic;

    priorMatch = 1-(nclus*ndays)./(nclus*nclus);
    leaveoutmatches = false(nclus,nclus,length(Scores2Include)); %Used later
    figure;
    if length(Scores2Include)>1
        for scid=1:length(Scores2Include)
            ScoresTmp = Scores2Include(scid);
            %             ScoresTmp(scid)=[];

            TotalScore = zeros(nclus,nclus);
            for scid2=1:length(ScoresTmp)
                eval(['TotalScore=TotalScore+' ScoresTmp{scid2} ';'])
            end
            base = length(ScoresTmp)-1;

            TotalScoreAcrossDays = TotalScore;
            TotalScoreAcrossDays(X==Y)=nan;

            subplot(2,length(Scores2Include),scid)
            h=imagesc(triu(TotalScore,1),[0 base+1]);
            title([Scores2Include{scid}])
            xlabel('Unit_i')
            ylabel('Unit_j')
            hold on
            arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
            arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
            colormap(flipud(gray))
            colorbar
            makepretty

            % Thresholds
            ThrsOpt = quantile(TotalScore(:),priorMatch); %Select best ones only later
            if ThrsOpt == max(TotalScore(:))
                ThrsOpt = ThrsOpt-0.1;
            end
            subplot(2,length(Scores2Include),scid+(length(Scores2Include)))
            leaveoutmatches(:,:,scid)=TotalScore>ThrsOpt;
            imagesc(triu(TotalScore>ThrsOpt,1))
            hold on
            title(['Thresholding at ' num2str(ThrsOpt)])
            xlabel('Unit_i')
            ylabel('Unit_j')
            hold on
            arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
            arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
            colormap(flipud(gray))
            colorbar
            %     axis square
            makepretty
        end
    end
    TotalScore = zeros(nclus,nclus);
    Predictors = zeros(nclus,nclus,0);
    for scid2=1:length(Scores2Include)
        Predictors = cat(3,Predictors,eval(Scores2Include{scid2}));
        eval(['TotalScore=TotalScore+' Scores2Include{scid2} ';'])
    end
    figure('name','TotalScore')
    subplot(2,1,1)
    imagesc(TotalScore,[0 length(Scores2Include)]);
    title('Total Score')
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    makepretty

    % Make initial threshold --> to be optimized
    ThrsOpt = quantile(TotalScore(:),priorMatch); %Select best ones only later
    subplot(2,1,2)
    imagesc(TotalScore>ThrsOpt)
    hold on
    title(['Thresholding at ' num2str(ThrsOpt)])
    xlabel('Unit_i')
    ylabel('Unit_j')
    hold on
    arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
    colormap(flipud(gray))
    colorbar
    % axis square
    makepretty
    saveas(gcf,fullfile(SaveDir,'TotalScore.fig'))
    saveas(gcf,fullfile(SaveDir,'TotalScore.bmp'))
    % Find all pairs
    % first factor authentication: score above threshold
    ThrsScore = ThrsOpt;
    % Take into account:
    label = TotalScore>ThrsOpt;
    [uid,uid2] = find(label);
    Pairs = cat(2,uid,uid2);
    Pairs = sortrows(Pairs);
    Pairs = unique(Pairs,'rows');
    Pairs(Pairs(:,1) == Pairs(:,2),:)=[];

    disp(['Computing total score took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

    %% Functional score for optimization: compute Fingerprint for the matched units - based on Célian Bimbard's noise-correlation finger print method but applied to across session correlations
    % Not every recording day will have the same units. Therefore we will
    % correlate each unit's activity with average activity across different
    % depths
    % Use a bunch of units with high total scores as reference population
    timercounter = tic;
    disp('Computing fingerprints correlations...')
    [~,sortid] = sort(cell2mat(arrayfun(@(X) TotalScore(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
    Pairs = Pairs(sortid,:);
    [FingerprintR,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(srAllDays,Pairs,Unit2Take,recsesGood);
%     CrossCorrelationFingerPrint_BU

    figure;
    h=scatter(TotalScore(:),FingerprintR(:),14,RankScoreAll(:),'filled','AlphaData',0.1);
    colormap(cat(1,[0 0 0],winter))
    xlabel('TotalScore')
    ylabel('Cross-correlation fingerprint')
    makepretty
    disp(['Computing fingerprints correlations took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

    %% three ways to define candidate scores
    % Total score larger than threshold
    CandidatePairs = TotalScore>ThrsOpt & RankScoreAll == 1 & SigMask == 1; %
    %     CandidatePairs(tril(true(size(CandidatePairs))))=0;
   

    %% Calculate median drift on this population (between days)
    if ndays>1
        for did = 1:ndays-1
            [uid,uid2] = find(CandidatePairs);
            BestPairs = cat(2,uid,uid2);
            idx = find(BestPairs(:,1)>=SessionSwitch(did)&BestPairs(:,1)<SessionSwitch(did+1) & BestPairs(:,2)>=SessionSwitch(did+1)&BestPairs(:,2)<SessionSwitch(did+2));
            if isempty(idx)
                disp('No pairs found to do any drift correction...')
                flag = 1;
                break
            end
            drift = nanmedian(nanmean(ProjectedLocation(:,BestPairs(idx,1),:),3)-nanmean(ProjectedLocation(:,BestPairs(idx,2),:),3),2);
            disp(['Median drift recording ' num2str(did) ' calculated: X=' num2str(drift(1)) ', Y=' num2str(drift(2))])
            if flag
                break
            end
            ProjectedLocation(1,GoodRecSesID==did+1,:)=ProjectedLocation(1,GoodRecSesID==did+1,:)+drift(1);
            ProjectedLocation(2,GoodRecSesID==did+1,:)=ProjectedLocation(2,GoodRecSesID==did+1,:)+drift(2);
            ProjectedLocationPerTP(1,GoodRecSesID==did+1,:,:) = ProjectedLocationPerTP(1,GoodRecSesID==did+1,:,:) + drift(1);
            ProjectedLocationPerTP(2,GoodRecSesID==did+1,:,:) = ProjectedLocationPerTP(2,GoodRecSesID==did+1,:,:) + drift(2);
            close all

        end
    else
        break
    end
    flag = flag+1;

end

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
xlabel('Units day 1')
ylabel('Units day 2')
title('Potential Matches')
makepretty
[uid,uid2] = find(CandidatePairs);
Pairs = cat(2,uid,uid2);
Pairs = sortrows(Pairs);
Pairs = unique(Pairs,'rows');

%% Naive bayes classifier
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
while flag<2 && runid<maxrun

    timercounter = tic;
    disp('Getting the Naive Bayes model...')

    flag = 0;
    runid = runid+1;
    filedir = which('UnitMatch');
    filedir = dir(filedir);
    if ApplyExistingBayesModel && exist(fullfile(filedir.folder,'UnitMatchModel.mat'))
        load(fullfile(SaveDir,'UnitMatchModel.mat'),'BestMdl')
        % Apply naive bays classifier
        Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters

        if isfield(BestMdl,'Parameterkernels')
            [label, posterior] = ApplyNaiveBayes(Tbl,BestMdl.Parameterkernels,[0 1],Priors);
        else
            [label, posterior, cost] = predict(BestMdl,Tbl);
        end
    else
        tmp= reshape(Predictors(Pairs(:,1),Pairs(:,2),:),[],length(Scores2Include));
        Tbl = array2table(reshape(tmp,[],size(tmp,2)),'VariableNames',Scores2Include); %All parameters
        % Use Rank as 'correct' label
        label = reshape(CandidatePairs(Pairs(:,1),Pairs(:,2)),1,[])';
       
        if MakeOwnNaiveBayes
            % Work in progress
            fprintf('Creating the Naive Bayes model...\n')
            timercounter2 = tic;
            [Parameterkernels,Performance] = CreateNaiveBayes(Tbl,label,Priors);
            fprintf('Done in %ds.\n', round(toc(timercounter2)))
            if any(Performance'<MaxPerf)
                flag = flag+1;
            else
                BestMdl.Parameterkernels = Parameterkernels;
            end
            % Apply naive bays classifier
            Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
            [label, posterior] = ApplyNaiveBayes(Tbl,Parameterkernels,[0 1],Priors);
            saveas(gcf,fullfile(SaveDir,'ProbabilityDistribution.fig'))
            saveas(gcf,fullfile(SaveDir,'ProbabilityDistribution.bmp'))
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
            saveas(gcf,fullfile(SaveDir,'ProbabilityDistribution.fig'))
            saveas(gcf,fullfile(SaveDir,'ProbabilityDistribution.bmp'))

        end
    end
    drawnow

    if runid<maxrun % Otherwise a waste of time!
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
        [FingerprintR,RankScoreAll,SigMask,AllSessionCorrelations] = CrossCorrelationFingerPrint(srAllDays,Pairs,Unit2Take,recsesGood);
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

%% If this was stitched pykilosort, we know what pykilosort thought about the matches
PyKSLabel = [];
PairsPyKS = [];
if RunPyKSChronicStitched
    for uid = 1:nclus
        pairstmp = find(AllClusterIDs(Good_Idx)==AllClusterIDs(Good_Idx(uid)))';
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
        eval(['ScoresTmp = ' Scores2Include{scid} ';'])
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
end

%% Extract final pairs:
disp('Extracting final pairs of units...')
timercounter = tic;
   
Tbl = array2table(reshape(Predictors,[],size(Predictors,3)),'VariableNames',Scores2Include); %All parameters
if isfield(BestMdl,'Parameterkernels')
    BestMdl.VariableNames = Scores2Include;
    Edges = [0:stepsize:1];
    ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
    BestMdl.ScoreVector = ScoreVector;
    if RunPyKSChronicStitched
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
Pairs=unique(Pairs,'rows');
% Pairs(Pairs(:,1)==Pairs(:,2),:)=[];
figure; imagesc(label)
colormap(flipud(gray))
xlabel('Unit_i')
ylabel('Unit_j')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
title('Identified matches')
makepretty
saveas(gcf,fullfile(SaveDir,'IdentifiedMatches.fig'))
saveas(gcf,fullfile(SaveDir,'IdentifiedMatches.bmp'))

disp(['Extracting final pair of units took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

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
%% Fingerprint correlations
% disp('Recalculate activity correlations')

% Use a bunch of units with high total scores as reference population
% [PairScore,sortid] = sort(cell2mat(arrayfun(@(X) MatchProbability(Pairs(X,1),Pairs(X,2)),1:size(Pairs,1),'Uni',0)),'descend');
% Pairs = Pairs(sortid,:);
% CrossCorrelationFingerPrint - do we really need this again?

%%
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

%% Extract final pairs:
label = MatchProbability>=param.ProbabilityThreshold;% | (MatchProbability>0.05 & RankScoreAll==1 & SigMask==1);
% label = MatchProbability>=0.99 | (MatchProbability>=0.05 & RankScoreAll==1 & SigMask==1);
[r, c] = find(triu(label,1)); %Find matches
Pairs = cat(2,r,c);
Pairs = sortrows(Pairs);
Pairs=unique(Pairs,'rows');
Pairs(Pairs(:,1)==Pairs(:,2),:)=[];

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
    MatchLkh = cell2mat(arrayfun(@(Y) Parameterkernels(minidx(Y,:),Y,2),1:size(minidx,1),'Uni',0));

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
saveas(gcf,fullfile(SaveDir,'ScoresSelfvsMatch.fig'))
saveas(gcf,fullfile(SaveDir,'ScoresSelfvsMatch.bmp'))

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
waveformdurationMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,waveformduration(X),waveformduration(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
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
SpatDecMat = arrayfun(@(Y) cell2mat(arrayfun(@(X) cat(1,spatialdecay(X),spatialdecay(Y)),1:nclus,'UniformOutput',0)),1:nclus,'UniformOutput',0);
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
save(fullfile(SaveDir,'MatchingScores.mat'),'BestMdl','SessionSwitch','GoodRecSesID','AllClusterIDs','Good_Idx','WavformMSE','WVCorr','LocationCombined','waveformTimePointSim','PeakTimeSim','spatialdecaySim','TotalScore','label','MatchProbability')
save(fullfile(SaveDir,'UnitMatchModel.mat'),'BestMdl')

%% Change these parameters to probabilities of being a match
for pidx = 1:length(Scores2Include)
    [~,IndividualScoreprobability, ~]=ApplyNaiveBayes(Tbl(:,pidx),Parameterkernels(:,pidx,:),Priors);
    eval([Scores2Include{pidx} ' = reshape(IndividualScoreprobability(:,2),nclus,nclus).*100;'])
end

%% Assign same Unique ID
OriUniqueID = UniqueID; %need for plotting
[PairID1,PairID2]=meshgrid(AllClusterIDs(Good_Idx));
[recses1,recses2] = meshgrid(recsesAll(Good_Idx));
[PairID3,PairID4]=meshgrid(OriUniqueID(Good_Idx));

MatchTable = table(PairID1(:),PairID2(:),recses1(:),recses2(:),PairID3(:),PairID4(:),MatchProbability(:),RankScoreAll(:),FingerprintR(:),TotalScore(:),'VariableNames',{'ID1','ID2','RecSes1','RecSes2','UID1','UID2','MatchProb','RankScore','FingerprintCor','TotalScore'});
UniqueID = AssignUniqueID(MatchTable,clusinfo,Path4UnitNPY,param);
% Add Scores2Include to MatchTable
for sid = 1:length(Scores2Include)
    eval(['MatchTable.' Scores2Include{sid} ' = ' Scores2Include{sid} '(:);'])
end

%% Other useful parameters to keep:
WaveformInfo.MaxChannel = MaxChannel;
WaveformInfo.ProjectedLocation = ProjectedLocation;
WaveformInfo.ProjectedWaveform = ProjectedWaveform;
WaveformInfo.ProjectedLocationPerTP = ProjectedLocationPerTP;

%% Figures
if MakePlotsOfPairs
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

    PlotTheseUnits_UM(Pairs(DrawPairs),MatchTable,WaveformInfo,AllSessionCorrelations,param,VisibleSetting)
end
return

