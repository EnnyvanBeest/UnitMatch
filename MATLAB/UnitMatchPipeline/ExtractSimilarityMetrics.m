function param = ExtractSimilarityMetrics(Scores2Include,AllWVBParameters,clusinfo,param,drawthis)


%% Extract fields and parameters
ExtractFields({AllWVBParameters})
waveidx = param.waveidx;
Allchannelpos = param.AllChannelPos;
AllchannelCoord = param.Coordinates;
SaveDir = param.SaveDir;
maxdist = param.maxdist;

if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')

end
GoodRecSesID = clusinfo.RecSesID(Good_Idx);
OriginalClusterIDs = clusinfo.cluster_id;

recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);

nclus = length(Good_Idx);
ndays = length(unique(GoodRecSesID));
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),unique(GoodRecSesID),'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch); nclus+1];
drift = nan;
% Do in batches
batchsz = 1000;
nbatch = ceil(nclus./batchsz);
SaveDrifts = nan(ndays,3,0);
if nargin<5
    if nclus > batchsz
        drawthis = 0;
    else
        drawthis = 1;
    end
end

% Save AUCs to help find best parameters
disp('Computing location distances between pairs of units...')
Loc2D = cat(3,ProjectedLocation(1,:,1)'-ProjectedLocation(1,:,2),ProjectedLocation(2,:,1)'-ProjectedLocation(2,:,2),ProjectedLocation(3,:,1)'-ProjectedLocation(3,:,2));
LocDist = squeeze(vecnorm(Loc2D,2,3)); % Distance


SameIdx = logical(eye(nclus));
WithinIdx = false(nclus,nclus);
for did = 1:ndays
    WithinIdx(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1) = true;
end
WithinIdx(LocDist>param.NeighbourDist) = false;
WithinIdx(logical(eye(size(WithinIdx)))) = false;
labels = [ones(1,sum(SameIdx(:))), zeros(1,sum(WithinIdx(:)))];
paramNames = {'waveformTimePointSim','spatialdecaySim','spatialdecayfitSim','AmplitudeSim','WVCorr','WavformMSE','WavformSim','CentroidDist','CentroidVar','CentroidDistRecentered','CentroidOverlord','TrajAngleSim','TrajDistSim','LocTrajectorySim'};
AUC = nan(1,length(paramNames));

%% Compute Metrics
disp('Computing Metric similarity between pairs of units...')
timercounter = tic;
% x1 = repmat(PeakTime(:,1),[1 numel(PeakTime(:,1))]);
% x2 = repmat(PeakTime(:,2),[1 numel(PeakTime(:,2))]);
% PeakTimeSim = abs(x1 - x2');
% %Normalize between 0 and 1 (values that make sense after testing, without having outliers influence this)
% PeakTimeSim =1-PeakTimeSim./quantile(PeakTimeSim(:),0.99);
% PeakTimeSim(PeakTimeSim<0)=0;

% % can't find much better for this one
% waveformTimePointSim = nan(nclus,nclus);
% for uid = 1:nclus
%     for uid2 = 1:nclus
%         waveformTimePointSim(uid,uid2) = sum(ismember(find(WaveIdx(uid,:,1)),find(WaveIdx(uid2,:,2))))./sum((WaveIdx(uid,:,1)));
%     end
% end
%
% scores = [waveformTimePointSim(SameIdx(:))', waveformTimePointSim(WithinIdx(:))'];
% paramid = find(ismember(paramNames,'waveformTimePointSim'));
% [~,~,~,AUC(paramid)] = perfcurve(labels,scores,1);

% Spatial decay
x1 = repmat(spatialdecay(:,1),[1 numel(spatialdecay(:,1))]);
x2 = repmat(spatialdecay(:,2),[1 numel(spatialdecay(:,2))]);
spatialdecaySim = abs(x1 - x2')./nanmean(cat(3,x1,x2'),3);
clear x1 x2
% Make (more) normal
spatialdecaySim = sqrt(spatialdecaySim);
spatialdecaySim = 1-((spatialdecaySim-nanmin(spatialdecaySim(:)))./(quantile(spatialdecaySim(:),0.99)-nanmin(spatialdecaySim(:))));
spatialdecaySim(spatialdecaySim<0)=0;

scores = [spatialdecaySim(SameIdx(:))', spatialdecaySim(WithinIdx(:))'];
paramid = find(ismember(paramNames,'spatialdecaySim'));
[~,~,~,AUC(paramid)] = perfcurve(labels,scores,1);

% spatial decay fit
x1 = repmat(spatialdecayfit(:,1),[1 numel(spatialdecayfit(:,1))]);
x2 = repmat(spatialdecayfit(:,2),[1 numel(spatialdecayfit(:,2))]);
spatialdecayfitSim = abs(x1 - x2');%./nanmean(cat(3,x1,x2'),3);
clear x1 x2
% Make (more) normal
spatialdecayfitSim = sqrt(spatialdecayfitSim);
spatialdecayfitSim = 1-((spatialdecayfitSim-nanmin(spatialdecayfitSim(:)))./(quantile(spatialdecayfitSim(:),0.99)-nanmin(spatialdecayfitSim(:))));
spatialdecayfitSim(spatialdecayfitSim<0 | isnan(spatialdecayfitSim))=0;


scores = [spatialdecayfitSim(SameIdx(:))', spatialdecayfitSim(WithinIdx(:))'];
paramid = find(ismember(paramNames,'spatialdecayfitSim'));
[~,~,~,AUC(paramid)] = perfcurve(labels,scores,1);

% Ampitude difference
x1 = repmat(Amplitude(:,1),[1 numel(Amplitude(:,1))]);
x2 = repmat(Amplitude(:,2),[1 numel(Amplitude(:,2))]);
AmplitudeSim = abs(x1 - x2')./nanmean(abs(cat(3,x1,x2')),3);
clear Amplitude
clear x1 x2
% Remove extreme outliers
% AmplitudeSim((AmplitudeSim)>quantile(AmplitudeSim(:),0.9999)) = quantile(AmplitudeSim(:),0.9999);

% Make (more) normal
AmplitudeSim = sqrt(AmplitudeSim);
AmplitudeSim = 1-((AmplitudeSim-nanmin(AmplitudeSim(:)))./(quantile(AmplitudeSim(:),.99)-nanmin(AmplitudeSim(:))));
AmplitudeSim(AmplitudeSim<0)=0;

scores = [AmplitudeSim(SameIdx(:))', AmplitudeSim(WithinIdx(:))'];
paramid = find(ismember(paramNames,'AmplitudeSim'));
[~,~,~,AUC(paramid)] = perfcurve(labels,scores,1);

%% Waveform similarity
timercounter = tic;
disp('Computing waveform similarity between pairs of units...')
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
WVCorr = (WVCorr-quantile(WVCorr(:),0.005))./(quantile(WVCorr(:),0.995)-quantile(WVCorr(:),0.005)); %Give WVCorr a better chance
WVCorr(WVCorr<0)=0;
WVCorr(WVCorr>1)=1;
scores = [WVCorr(SameIdx(:))', WVCorr(WithinIdx(:))'];
paramid = find(ismember(paramNames,'WVCorr'));
[x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);


ProjectedWaveformNorm = cat(3,x1,x2);
clear x1 x2
ProjectedWaveformNorm = (ProjectedWaveformNorm-nanmin(ProjectedWaveformNorm,[],1))./(nanmax(ProjectedWaveformNorm,[],1)-nanmin(ProjectedWaveformNorm,[],1));
x1 = repmat(ProjectedWaveformNorm(:,:,1),[1 1 size(ProjectedWaveformNorm,2)]);
x2 = permute(repmat(ProjectedWaveformNorm(:,:,2),[1 1 size(ProjectedWaveformNorm,2)]),[1 3 2]);
RawWVMSE = squeeze(nanmean((x1 - x2).^2));
clear x1 x2

% sort of Normalize distribution
RawWVMSENorm = sqrt(RawWVMSE);
WavformMSE = (RawWVMSENorm-nanmin(RawWVMSENorm(:)))./(quantile(RawWVMSENorm(:),0.99)-nanmin(RawWVMSENorm(:)));
WavformMSE = 1-WavformMSE;
WavformMSE(WavformMSE<0) = 0;
scores = [WavformMSE(SameIdx(:))', WavformMSE(WithinIdx(:))'];
paramid = find(ismember(paramNames,'WavformMSE'));
[x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

disp(['Calculating waveform similarity took ' num2str(round(toc(timercounter))) ' seconds for ' num2str(nclus) ' units'])
figure('name','Waveform similarity measures')
subplot(1,3,1)
imagesc(WVCorr,[0.5 0.9]);
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
imagesc(WavformMSE,[0.5 0.9]);
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
scores = [WavformSim(SameIdx(:))', WavformSim(WithinIdx(:))'];
paramid = find(ismember(paramNames,'WavformSim'));
[x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

subplot(1,3,3)
imagesc(WavformSim,[0.5 0.9]);
title('Average Waveform scores')
xlabel('Unit Y')
ylabel('Unit Z')
hold on
arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
colormap(flipud(gray))
colorbar
makepretty

%% Location differences between pairs of units: - This is done twice to account for large drift between sessions
channelpos_AllCat = unique(cat(1,AllchannelCoord{:}),'rows');

%% Flip trajectory if necessary?
% for which dimensions do we allow flipping?
AllowFlipping = false(size(channelpos_AllCat,2),nclus); % Dimension x channel
for uid = 1:nclus
    if param.RunPyKSChronicStitched
        channelpos = AllchannelCoord{1};
    else
        channelpos = AllchannelCoord{recsesGood(uid)};
    end
    if isnan(MaxChannel(uid,1))
        continue
    end
    %Load channels
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = round(channelpos(ChanIdx,:)./50).*50;
    AllowFlipping(cell2mat(arrayfun(@(X) length(unique(Locs(:,X))),1:size(Locs,2),'Uni',0))<=2 & cell2mat(arrayfun(@(X) length(unique(channelpos(ChanIdx,X))),1:size(Locs,2),'Uni',0))>1,:) = true;
end
FlipDim = find(any(AllowFlipping,2));
% Housekeeping
clear AllowFlipping
clear AllchannelCoord
% Flip trajectory for flip dimensions
ProjectedLocationPerTPAllFlips = nan(size(ProjectedLocationPerTP,1),size(ProjectedLocationPerTP,2),size(ProjectedLocationPerTP,3),size(ProjectedLocationPerTP,4),length(FlipDim));
for flipid = 1:length(FlipDim)
    tmpdat = squeeze(ProjectedLocationPerTP(FlipDim(flipid),:,:,:));
    range = cat(2,nanmin(tmpdat,[],2), nanmax(tmpdat,[],2));

    % change values
    newvals = nanmin(tmpdat,[],2) + (nanmax(tmpdat,[],2) - tmpdat);
    ProjectedLocationPerTPAllFlips(:,:,:,:,flipid) = ProjectedLocationPerTP;
    ProjectedLocationPerTPAllFlips(FlipDim(flipid),:,:,:,flipid) = newvals;
end
% Housekeeping
clear newvals
clear tmpdat
clear range
ProjectedLocationPerTPAllFlips = cat(5,ProjectedLocationPerTP,ProjectedLocationPerTPAllFlips); % add them all together

flag=0;
while flag<2
    timercounter = tic;
    if flag && drawthis
        figure('name','Projection locations all units')
        scatter3(channelpos_AllCat(:,1),channelpos_AllCat(:,2),channelpos_AllCat(:,3),10,[0 0 0],'filled')
        hold on
        scatter3(nanmean(ProjectedLocation(1,:,:),3),nanmean(ProjectedLocation(2,:,:),3),nanmean(ProjectedLocation(3,:,:),3),10,GoodRecSesID)
        colormap jet
        %     set(gca,'CameraPosition',[-1.9009e+04 -3.7880e+04 1.0753e+04])
        makepretty
        ylabel('YPos (um)')
        zlabel('Depth (um)')
        xlabel('XPos (um)')

        %     xlim([min(channelpos_AllCat(:,1))-50 max(channelpos_AllCat(:,1))+50])
        drawnow
        saveas(gcf,fullfile(SaveDir,'ProjectedLocation.fig'))
        saveas(gcf,fullfile(SaveDir,'ProjectedLocation.bmp'))
    end
    disp('Computing location distances between pairs of units, per individual time point of the waveform...')
    % Difference in distance between centroids of two halves of the recording

    % Clean up variables; this gets intense
    clear x
    clear X
    clear Y
    clear y

    % memory efficient?
    EuclDist = nan(nclus,length(waveidx),length(FlipDim)+1,nclus);
    for batchid1 = 1:nbatch
        idx = (batchid1-1)*batchsz+1:batchsz*batchid1;
        idx(idx>nclus) = [];
        for batchid2 = 1:nbatch

            idx2 = (batchid2-1)*batchsz+1:batchsz*batchid2;
            idx2(idx2>nclus) = [];

            x1 = repmat(squeeze(ProjectedLocationPerTPAllFlips(:,idx,waveidx,1,:)),[1 1 1 1 numel(idx2)]);
            x2 = permute(repmat(squeeze(ProjectedLocationPerTPAllFlips(:,idx2,waveidx,2,:)),[1 1 1 1 numel(idx)]),[1,5,3,4,2]);%Switch the two nclus around

            w = squeeze(isnan(abs(x1(1,:,:,:,:)-x2(1,:,:,:,:))));
            tmpEu = squeeze(vecnorm(x1-x2,2,1)); % Distance
            tmpEu(w) = nan;
            EuclDist(idx,:,:,idx2) = tmpEu; % Euclidean distance
        end
    end


    clear x1 x2
    % Average location
    CentroidDist = squeeze(nanmin(squeeze(EuclDist(:,param.NewPeakLoc-param.waveidx==0,:,:)),[],2));%


    % Normalize each of them from 0 to 1, 1 being the 'best'
    % If distance > maxdist micron it will never be the same unit:
    CentroidDist = 1-((CentroidDist-nanmin(CentroidDist(:)))./(maxdist-nanmin(CentroidDist(:)))); %Average difference
    CentroidDist(CentroidDist<0)=0;
    CentroidDist(isnan(CentroidDist))=0;

    scores = [CentroidDist(SameIdx(:))', CentroidDist(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'CentroidDist'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);


    % Variance in error, corrected by average error. This captures whether
    % the trajectory is consistenly separate
    CentroidVar = squeeze(min(squeeze(nanvar(EuclDist,[],2)),[],2));%./nanmean(EuclDist,2)+nanmean(EuclDist,2));
    CentroidVar = sqrt(CentroidVar);
    CentroidVar = 1-((CentroidVar-nanmin(CentroidVar(:)))./(quantile(CentroidVar(:),0.99)-nanmin(CentroidVar(:)))); %Average difference
    CentroidVar(CentroidVar<0) = 0;
    CentroidVar(isnan(CentroidVar)) = 0;
    scores = [CentroidVar(SameIdx(:))', CentroidVar(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'CentroidVar'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);
    % @Célian suggestion: recenter to 0 first, then calculate the
    % difference in distance (to account for uncorrected drift)
    %     ProjectedLocationPerTPRecentered = permute(permute(ProjectedLocationPerTP,[1,2,4,3]) - ProjectedLocation,[1,2,4,3]);
    %     x1 = repmat(squeeze(ProjectedLocationPerTPRecentered(:,:,waveidx,1)),[1 1 1 size(ProjectedLocationPerTPRecentered,2)]);
    %     x2 = permute(repmat(squeeze(ProjectedLocationPerTPRecentered(:,:,waveidx,2)),[1 1 1 size(ProjectedLocationPerTPRecentered,2)]),[1 4 3 2]);
    %     EuclDist2 = squeeze(sqrt(nansum((x1-x2).^2,1))); % Euclidean distance
    %     w = squeeze(isnan(abs(x1(1,:,:,:)-x2(1,:,:,:))));
    %     EuclDist2(w) = nan;
    disp('Computing location distances between pairs of units, per individual time point of the waveform, Recentered...')
    ProjectedLocationPerTPRecentered = permute(permute(ProjectedLocationPerTPAllFlips,[1,2,4,3,5]) - ProjectedLocation,[1,2,4,3,5]);
    EuclDist2 = nan(nclus,length(waveidx),length(FlipDim)+1,nclus);
    for batchid1 = 1:nbatch
        idx = (batchid1-1)*batchsz+1:batchsz*batchid1;
        idx(idx>nclus) = [];
        for batchid2 = 1:nbatch

            idx2 = (batchid2-1)*batchsz+1:batchsz*batchid2;
            idx2(idx2>nclus) = [];

            x1 = repmat(squeeze(ProjectedLocationPerTPRecentered(:,idx,waveidx,1,:)),[1 1 1 1 numel(idx2)]);
            x2 = permute(repmat(squeeze(ProjectedLocationPerTPRecentered(:,idx2,waveidx,2,:)),[1 1 1 1 numel(idx)]),[1,5,3,4,2]);%Switch the two nclus around

            w = squeeze(isnan(abs(x1(1,:,:,:,:)-x2(1,:,:,:,:))));
            tmpEu = squeeze(vecnorm(x1-x2,2,1)); % Distance
            tmpEu(w) = nan;
            EuclDist2(idx,:,:,idx2) = tmpEu; % Euclidean distance
        end
    end


    clear x1 x2
    % Average location
    CentroidDistRecentered = squeeze(nanmin(nanmean(EuclDist2,2),[],3));% minimum across flips
    CentroidDistRecentered = 1-(CentroidDistRecentered-nanmin(CentroidDistRecentered(:)))./(quantile(CentroidDistRecentered(:),0.99)-nanmin(CentroidDistRecentered(:)));
    CentroidDistRecentered(CentroidDistRecentered<0|isnan(CentroidDistRecentered))=0;
    scores = [CentroidDistRecentered(SameIdx(:))', CentroidDistRecentered(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'CentroidDistRecentered'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);



    CentroidOverlord = (CentroidDistRecentered+CentroidVar)/2;
    scores = [CentroidOverlord(SameIdx(:))', CentroidOverlord(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'CentroidOverlord'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

    disp('Computing location angle (direction) differences between pairs of units, per individual time point of the waveform...')
    x1 = ProjectedLocationPerTPAllFlips(:,:,waveidx(2):waveidx(end),:,:);
    x2 = ProjectedLocationPerTPAllFlips(:,:,waveidx(1):waveidx(end-1),:,:);
    % The distance traveled (Eucledian)
    TrajDist = squeeze(vecnorm(x1-x2,2,1));
    %Use  TrajDist to select angle ehere there is a minimum amount movement
    good_ang = zeros(size(TrajDist));
    good_ang(TrajDist >= param.min_angledist) = 1;
    % Difference in angle between two time points
    LocAngle = nan(size(TrajDist,1),size(TrajDist,2),size(TrajDist,3),size(TrajDist,4),0);
    countid=1;
    for dimid1=1:size(ProjectedLocationPerTPAllFlips,1)
        for dimid2=2:size(ProjectedLocationPerTPAllFlips,1)
            if dimid2<=dimid1
                continue
            end
            LocAngle(:,:,:,:,countid) = squeeze(atan(abs(x1(dimid1,:,:,:,:)-x2(dimid1,:,:,:,:))./abs(x1(dimid2,:,:,:,:)-x2(dimid2,:,:,:,:)))) .* good_ang;
            countid = countid + 1;
        end
    end

    % Sum the angles across dimensions
    LocAngle = nansum(LocAngle,5);

    clear x1 x2

    % Actually just taking the weighted sum of angles is better
    x1 = repmat(squeeze(LocAngle(:,:,1,:)),[1 1 1 nclus]);
    x2 = permute(repmat(squeeze(LocAngle(:,:,2,1)),[1 1 1 nclus]),[4 2 3 1]); %switch nclus around
    AngleSubtraction = abs(x1-x2);
    AngleSubtraction(isnan(abs(x1-x2))) = 2*pi; %punish points with nan
    clear x1 x2
    TrajAngleSim = squeeze(nanmin(nansum(AngleSubtraction,2),[],3)); % sum of angles, minimum across flips
    TrajAngleSim = 1-((TrajAngleSim-nanmin(TrajAngleSim(:)))./(quantile(TrajAngleSim(:),0.99)-nanmin(TrajAngleSim(:))));
    TrajAngleSim(TrajAngleSim<0 | isnan(TrajAngleSim))=0;
    scores = [TrajAngleSim(SameIdx(:))', TrajAngleSim(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'TrajAngleSim'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

    % Continue distance traveled
    x1 = repmat(squeeze(TrajDist(:,:,1,:)),[1 1 1 nclus]);
    x2 = permute(repmat(squeeze(TrajDist(:,:,2,:)),[1 1 1 nclus]),[4 2 3 1]); % switch nclus around
    % Distance similarity (subtract for each pair of units)
    TrajDistCompared = abs(x1-x2);%
    clear x1 x2
    TrajDistSim = squeeze(nanmin(nansum(TrajDistCompared,2),[],3)); %and take minimum across flips
    TrajDistSim = sqrt(TrajDistSim); % Make more normal
    TrajDistSim = 1-((TrajDistSim-nanmin(TrajDistSim(:)))./(quantile(TrajDistSim(:),0.99)-nanmin(TrajDistSim(:))));
    TrajDistSim(TrajDistSim<0 | isnan(TrajDistSim))=0;
    scores = [TrajDistSim(SameIdx(:))', TrajDistSim(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'TrajDistSim'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);


    LocTrajectorySim = (TrajAngleSim+TrajDistSim)./2; % Trajectory Similarity is sum of distance + sum of angles
    LocTrajectorySim = (LocTrajectorySim-nanmin(LocTrajectorySim(:))./(quantile(LocTrajectorySim(:),0.99)-nanmin(LocTrajectorySim(:))));
    scores = [LocTrajectorySim(SameIdx(:))', LocTrajectorySim(WithinIdx(:))'];
    paramid = find(ismember(paramNames,'LocTrajectorySim'));
    [x,y,~,AUC(paramid)] = perfcurve(labels,scores,1);

    %
    if flag  && drawthis
        figure('name','Distance Measures')
        subplot(5,2,1)
        imagesc(CentroidDist);
        title('CentroidDist')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty

        subplot(5,2,2)
        tmp = CentroidDist;
        hc1 = histcounts(diag(tmp),[0:0.01:1])./size(tmp,1);
        hc2 = histcounts(tmp(~logical(eye(size(tmp)))),[0:0.01:1])./(numel(tmp)-size(tmp,1));
        plot([0.005:0.01:0.995],hc2); hold on; plot([0.005:0.01:0.995],hc1)
        xlabel('Score')
        makepretty

        subplot(5,2,3)
        imagesc(TrajAngleSim);
        title('TrajAngleSim')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty
        subplot(5,2,4)
        tmp = TrajAngleSim;
        hc1 = histcounts(diag(tmp),[0:0.01:1])./size(tmp,1);
        hc2 = histcounts(tmp(~logical(eye(size(tmp)))),[0:0.01:1])./(numel(tmp)-size(tmp,1));
        plot([0.005:0.01:0.995],hc2); hold on; plot([0.005:0.01:0.995],hc1)
        xlabel('Score')
        makepretty

        subplot(5,2,5)
        imagesc(TrajDistSim);
        title('TrajDistSim')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty
        subplot(5,2,6)
        tmp = TrajDistSim;
        hc1 = histcounts(diag(tmp),[0:0.01:1])./size(tmp,1);
        hc2 = histcounts(tmp(~logical(eye(size(tmp)))),[0:0.01:1])./(numel(tmp)-size(tmp,1));
        plot([0.005:0.01:0.995],hc2); hold on; plot([0.005:0.01:0.995],hc1)
        xlabel('Score')
        makepretty

        subplot(5,2,7)
        imagesc(CentroidVar);
        title('CentroidVar')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty
        subplot(5,2,8)
        tmp = CentroidVar;
        hc1 = histcounts(diag(tmp),[0:0.01:1])./size(tmp,1);
        hc2 = histcounts(tmp(~logical(eye(size(tmp)))),[0:0.01:1])./(numel(tmp)-size(tmp,1));
        plot([0.005:0.01:0.995],hc2); hold on; plot([0.005:0.01:0.995],hc1)
        xlabel('Score')
        makepretty

        subplot(5,2,9)
        imagesc(CentroidDistRecentered);
        title('CentroidDistRecentered')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty

        subplot(5,2,10)
        tmp = CentroidDistRecentered;
        hc1 = histcounts(diag(tmp),[0:0.01:1])./size(tmp,1);
        hc2 = histcounts(tmp(~logical(eye(size(tmp)))),[0:0.01:1])./(numel(tmp)-size(tmp,1));
        plot([0.005:0.01:0.995],hc2); hold on; plot([0.005:0.01:0.995],hc1)
        xlabel('Score')
        makepretty
    end

    disp(['Extracting projected location took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
    % Average EuclDist
    EuclDist = squeeze(nanmin(squeeze(EuclDist(:,param.NewPeakLoc-param.waveidx==0,:,:)),[],2));%

    % Plotting order (sort units based on distance)
    [~,SortingOrder] = arrayfun(@(X) sort(EuclDist(1,SessionSwitch(X):SessionSwitch(X+1)-1)),1:ndays,'Uni',0);
    SortingOrder = arrayfun(@(X) squeeze(SortingOrder{X}+SessionSwitch(X)-1),1:ndays,'Uni',0);
    if size(SortingOrder{1},1)==1
        SortingOrder = cat(2,SortingOrder{:});
    else
        SortingOrder = cat(1,SortingOrder{:});
    end
    %% These are the parameters to include:
    if drawthis && flag
        figure('name','Total Score components');
        for sid = 1:length(Scores2Include)
            eval(['tmp = ' Scores2Include{sid} ';'])
            tmp(EuclDist>param.NeighbourDist)=nan; %Remove pairs that are never gonna happen

            subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),sid)
            try
                imagesc(tmp(SortingOrder,SortingOrder));
            catch
                imagesc(tmp(SortingOrder,SortingOrder),[0 1]);
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
        try
            saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.fig'))
        catch ME
        end
        try
            saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.bmp'))
        catch ME
        end


        figure('name','Diagonal versus off-diagonal')
        for sid = 1:length(Scores2Include)
            eval(['tmp = ' Scores2Include{sid} ';'])
            % Take centroid dist > maxdist out
            tmp(EuclDist>param.NeighbourDist)=nan;
            % Take between session out
            for did = 1:ndays
                for did2 = 1:ndays
                    if did==did2
                        continue
                    end
                    tmp(SessionSwitch(did):SessionSwitch(did+1),SessionSwitch(did2):SessionSwitch(did2+1))=nan;
                end
            end
            hd = histcounts(diag(tmp),0:0.01:1)./nclus;
            hnd = histcounts(tmp(~eye(size(tmp))),0:0.01:1)./sum(~isnan(tmp(~eye(size(tmp)))));
            subplot(ceil(sqrt(length(Scores2Include))),round(sqrt(length(Scores2Include))),sid)
            plot(0.005:0.01:0.995,hd,'r-'); hold on; plot(0.005:0.01:0.995,hnd,'b-')

            title(Scores2Include{sid})
            makepretty

        end
        legend('Within session diagonal','Within session off-diagonal (<maxdist)','Location', 'best')
        saveas(gcf,fullfile(SaveDir,'WithinSessionDistributions.fig'))
        saveas(gcf,fullfile(SaveDir,'WithinSessionDistributions.bmp'))

    end


    %% Calculate total score
    [X,Y]=meshgrid(recsesAll(Good_Idx));

    IncludeThesePairs = find(EuclDist<maxdist);

    disp('Computing total score...')
    timercounter = tic;

    TotalScore = zeros(nclus,nclus);
    Predictors = zeros(nclus,nclus,0);
    for sid=1:length(Scores2Include)
        eval(['tmp = ' Scores2Include{sid} ';'])
        % Take centroid dist > maxdist out
        %         tmp(EuclDist>param.NeighbourDist)=nan;
        Predictors = cat(3,Predictors,tmp);
        TotalScore=TotalScore+tmp;
    end

    % Normalize TotalScore
    TotalScore = (TotalScore-min(TotalScore(:)))./(max(TotalScore(:))-min(TotalScore(:)));

    if flag  && drawthis
        figure('name','Predictor Matrix of max dist pairs')
        %  Score correlation matrix
        tmp = reshape(Predictors,nclus*nclus,[]);
        [~,ax] = plotmatrix(tmp(IncludeThesePairs,:));
        for sid = 1:size(ax,1)
            set(ax(sid,1),'YTick',0.5,'YTickLabel',Scores2Include{sid})
            set(ax(end,sid),'XTick',0.5,'XTickLabel',Scores2Include{sid},'XTickLabelRotation',45)
        end
        makepretty

        figure('name','TotalScore')
        subplot(2,2,1)
        imagesc(TotalScore(SortingOrder,SortingOrder),[0 1]);
        title('Total Score')
        xlabel('Unit_i')
        ylabel('Unit_j')
        hold on
        arrayfun(@(X) line([SessionSwitch(X) SessionSwitch(X)],get(gca,'ylim'),'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        arrayfun(@(X) line(get(gca,'xlim'),[SessionSwitch(X) SessionSwitch(X)],'color',[1 0 0]),2:length(SessionSwitch),'Uni',0)
        colormap(flipud(gray))
        colorbar
        makepretty
    end
    %         hold on
    %         scatter(find(SortingOrder==Pairs(1)),find(SortingOrder==Pairs(3)),5,cols(3,:),'filled')
    %         scatter(find(SortingOrder==Pairs(1)),find(SortingOrder==Pairs(2)),5,cols(2,:),'filled')


    %% Make initial threshold --> to be optimized
    %     Take initial distributions for same and neighbors within day
    stepsz = 0.01;
    ScoreVector = stepsz./2:stepsz:1-stepsz./2;
    Bins = 0:stepsz:1;

    tmp = TotalScore;
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
    hd = histcounts(diag(tmp),Bins)./nclus;
    hnd = histcounts(tmp(~eye(size(tmp))),Bins)./sum(~isnan(tmp(~eye(size(tmp)))));
    hp = histcounts(tmp(:),Bins)./sum(~isnan(tmp(:)));
    if size(ScoreVector) ~= size(hd)
        ScoreVector = ScoreVector';
    end
    ThrsOpt = ScoreVector(find(smoothdata(hd)>smoothdata(hnd)&ScoreVector>0.6,1,'first'));
    [muw, sw] = normfit(tmp(~isnan(tmp) & tmp<ThrsOpt));



    tmp = TotalScore;
    % Take centroid dist > maxdist out
    tmp(EuclDist>param.NeighbourDist)=nan;
    % Take within session out
    for did = 1:ndays
        tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)=nan;
    end
    ha = histcounts(tmp(:),Bins)./sum(~isnan(tmp(:)));
    [mua, sa] = normfit(tmp(~isnan(tmp)  & tmp<ThrsOpt));

    if ~isnan(mua) && mua<muw && ~flag
        ThrsOpt = ThrsOpt-abs(muw-mua); % Correct for general scores being lower across days (e.g. unresolved drift)
    end
    if flag  && drawthis
        subplot(2,2,3)
        plot(ScoreVector,hd,'-','color',[0 0.7 0]); hold on; plot(ScoreVector,hnd,'b-')
        plot(ScoreVector,ha,'-','color',[1 0 0]);
        %     plot(ScoreVector,hp,'--','color',[1 0 0]);

        line([ThrsOpt ThrsOpt],get(gca,'ylim'),'LineStyle','--','color',[0 0 0])
        xlabel('TotalScore')
        ylabel('number of pairs')
        makepretty
        title('Total score distributions')
        legend('Same Unit','Neighbors','Across','Threshold','Location', 'best')


        subplot(2,2,2)
        imagesc(TotalScore(SortingOrder,SortingOrder)>ThrsOpt)
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
    end
    %     hold on
    %     scatter(find(SortingOrder==Pairs(1)),find(SortingOrder==Pairs(3)),5,cols(3,:),'filled')
    %     scatter(find(SortingOrder==Pairs(1)),find(SortingOrder==Pairs(2)),5,cols(2,:),'filled')


    param.nExpectedMatches = sum(TotalScore(:)>ThrsOpt);
    if param.nExpectedMatches==0
        keyboard
    end
    priorMatch = 1-(param.nExpectedMatches./length(IncludeThesePairs)); %Punish multiple days (unlikely to find as many matches after a few days) % Times 2 for symmetry

    %% Cumulative density function
    if flag  && drawthis
        subplot(2,2,4)
        [h,stats] = cdfplot(TotalScore(IncludeThesePairs));
        h.Color = [0 0 0];
        hold on
        tmp = TotalScore;
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
        [h,stats] = cdfplot(tmp(~isnan(tmp))); %Neighbours
        h.Color = [0 0 0.5];

        tmp = TotalScore;
        % Take centroid dist > maxdist out
        tmp(EuclDist>param.NeighbourDist)=nan;
        % Take within session out
        for did = 1:ndays
            tmp(SessionSwitch(did):SessionSwitch(did+1)-1,SessionSwitch(did):SessionSwitch(did+1)-1)=nan;
        end
        if any(~isnan(tmp(:)))
            [h,stats] = cdfplot(tmp(:)); %pAcross sessions
            h.Color = [1 0 0];
        end
        %     [h,stats] = cdfplot(tmp(tmp(:)<ThrsOpt)); %putative matches
        %     h.Color = [0.5 0.2 0];

        xlabel('TotalScore')
        ylabel('Cumulative density')
        line([ThrsOpt,ThrsOpt],[0 1],'color',[1 0 0],'LineStyle','--')
        makepretty

        legend('All pairs','Same unit','Neighbors','Across','threshold')

        saveas(gcf,fullfile(SaveDir,'TotalScore.fig'))
        saveas(gcf,fullfile(SaveDir,'TotalScore.bmp'))
    end

    %% More plots
    if flag  && drawthis
        leaveoutmatches = false(nclus,nclus,length(Scores2Include)); %Used later
        figure;
        if length(Scores2Include)>1
            for scid=1:length(Scores2Include)
                ScoresTmp = Scores2Include(scid);
                %             ScoresTmp(scid)=[];

                TotalScoreLO = zeros(nclus,nclus);
                for scid2=1:length(ScoresTmp)
                    eval(['TotalScoreLO=TotalScoreLO+' ScoresTmp{scid2} ';'])
                end
                base = length(ScoresTmp)-1;

                TotalScoreAcrossDays = TotalScoreLO;
                TotalScoreAcrossDays(X==Y)=nan;

                subplot(2,length(Scores2Include),scid)
                h=imagesc(TotalScoreLO,[0 base+1]);
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
                ThrsOptLO = quantile(TotalScoreLO(IncludeThesePairs),priorMatch); %Select best ones only later
                if ThrsOptLO == max(TotalScoreLO(IncludeThesePairs))
                    ThrsOptLO = ThrsOptLO-0.1;
                end
                subplot(2,length(Scores2Include),scid+(length(Scores2Include)))
                leaveoutmatches(:,:,scid)=TotalScoreLO>ThrsOptLO;
                imagesc(TotalScoreLO>ThrsOptLO)
                hold on
                title(['Thresholding at ' num2str(ThrsOptLO)])
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
    end
    %% three ways to define candidate scores
    % Total score larger than threshold
    CandidatePairs = TotalScore>ThrsOpt;%
    TheseDrifts = nan(ndays,3);
    %% Calculate median drift on this population (between days)
    if ndays>1
        for did = 1:ndays-1
            [uid,uid2] = find(CandidatePairs);
            BestPairs = cat(2,uid,uid2);
            idx = find(BestPairs(:,1)>=SessionSwitch(did)&BestPairs(:,1)<SessionSwitch(did+1) & BestPairs(:,2)>=SessionSwitch(did+1)&BestPairs(:,2)<SessionSwitch(did+2));
            if isempty(idx) || length(idx)<3
                disp('No pairs found to do any drift correction...')
                continue
            end
            drift = nanmedian(nanmean(ProjectedLocation(:,BestPairs(idx,1),:),3)-nanmean(ProjectedLocation(:,BestPairs(idx,2),:),3),2);
            disp(['Median drift recording ' num2str(did) ' calculated: X=' num2str(drift(1)) ', Y=' num2str(drift(2)) ', Z=' num2str(drift(3))])
            TheseDrifts(did,:) = drift;
            if ~flag

                ProjectedLocation(1,GoodRecSesID==did+1,:)=ProjectedLocation(1,GoodRecSesID==did+1,:)+drift(1);
                ProjectedLocation(2,GoodRecSesID==did+1,:)=ProjectedLocation(2,GoodRecSesID==did+1,:)+drift(2);
                ProjectedLocation(3,GoodRecSesID==did+1,:)=ProjectedLocation(3,GoodRecSesID==did+1,:)+drift(3);

                ProjectedLocationPerTP(1,GoodRecSesID==did+1,:,:) = ProjectedLocationPerTP(1,GoodRecSesID==did+1,:,:) + drift(1);
                ProjectedLocationPerTP(2,GoodRecSesID==did+1,:,:) = ProjectedLocationPerTP(2,GoodRecSesID==did+1,:,:) + drift(2);
                ProjectedLocationPerTP(3,GoodRecSesID==did+1,:,:) = ProjectedLocationPerTP(3,GoodRecSesID==did+1,:,:) + drift(3);

                ProjectedLocationPerTPAllFlips(1,GoodRecSesID==did+1,:,:,:) = ProjectedLocationPerTPAllFlips(1,GoodRecSesID==did+1,:,:,:) + drift(1);
                ProjectedLocationPerTPAllFlips(2,GoodRecSesID==did+1,:,:,:) = ProjectedLocationPerTPAllFlips(2,GoodRecSesID==did+1,:,:,:) + drift(2);
                ProjectedLocationPerTPAllFlips(3,GoodRecSesID==did+1,:,:,:) = ProjectedLocationPerTPAllFlips(3,GoodRecSesID==did+1,:,:,:) + drift(3);

                close all

            end

        end
    else
        break
    end
    flag = flag+1;
    SaveDrifts = cat(3,SaveDrifts,TheseDrifts);
end
%% Store in param
param.drift = SaveDrifts;
%% Assign to workspace
assignin('caller','TotalScore',TotalScore)
assignin('caller','Predictors',Predictors)
assignin('caller','ProjectedLocationPerTP',ProjectedLocationPerTP)
assignin('caller','ProjectedLocation',ProjectedLocation)
assignin('caller','EuclDist',EuclDist)
assignin('caller','SortingOrder',SortingOrder)
AUCStruct.AUC = AUC;
AUCStruct.ParamNames = paramNames;
save(fullfile(SaveDir,'AUC.mat'),'AUCStruct')
return
if 0 % THis can be used to look at some example projections
    % Find all pairs
    % first factor authentication: score above threshold
    % Take into account:
    label = TotalScore>ThrsOpt;
    [uid,uid2] = find(label);
    Pairs = cat(2,uid,uid2);
    Pairs = sortrows(Pairs);
    Pairs = unique(Pairs,'rows');
    Pairs(Pairs(:,1) == Pairs(:,2),:)=[];
    %% Plot
    Pairs = [198, 469, 47] % Example (AL032, take 10)

    %     Pairs = [10,450,11] % Example
    cols =  [0 0 0; 1 0 0; 0 0 0.7];

    figure
    subplot(1,3,3)
    for uidx=1:length(Pairs)
        uid = Pairs(uidx);
        channelpos = Allchannelpos{recsesGood(uid)};
        % Load raw data
        try
            spikeMap = readNPY(fullfile(param.KSDir{recsesGood(uid)},'RawWaveforms',['Unit' num2str(OriginalClusterIDs(Good_Idx(uid))) '_RawSpikes.npy'])); %0-indexed
        catch
            keyboard
        end
        % Detrending
        spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
        spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
        spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order
        %Load channels
        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius.*0.4); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);

        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled','marker','s')
        hold on

        takesamples = param.waveidx;
        takesamples = unique(takesamples(~isnan(squeeze(ProjectedLocationPerTP(2,uid,takesamples,1)))));
        h(1) = plot(squeeze(ProjectedLocationPerTP(2,uid,takesamples,1)),squeeze(ProjectedLocationPerTP(3,uid,takesamples,1)),'-','color',cols(uidx,:));
        scatter(squeeze(ProjectedLocationPerTP(2,uid,takesamples(1),1)),squeeze(ProjectedLocationPerTP(3,uid,takesamples(1),1)),30,cols(uidx,:),'filled')
        plot(ProjectedLocation(2,uid,1),ProjectedLocation(3,uid,1),'.','MarkerSize',25,'color',cols(uidx,:))

    end

    %     MyColMap = hsv(length(takesamples)*2);
    %     colormap(MyColMap(length(takesamples)+1:end,:))

    xlabel('Xpos (\mum)')
    ylabel('Ypos (\mum)')
    ylims = [min(Locs(:,2))-15 max(Locs(:,2))+15];

    xlims = [min(ProjectedLocation(2,Pairs,1))-diff(ylims)/2 min(ProjectedLocation(2,Pairs,1))+diff(ylims)/2];
    ylabel('')
    set(gca,'xlim',xlims,'ylim',ylims)
    axis square
    %     hpos = get(gca,'Position')
    %         legend([h(1),h(2)],{['Unit ' num2str(uid)],['Unit ' num2str(uid2)]})
    %     hc= colorbar;
    %     try
    %         hc.Label.String = 'timesample';
    %     catch ME
    %         disp(ME)
    %         keyboard
    %     end
    title('Centroid trajectory')
    makepretty



    subplot(1,3,2)
    patch([-(1/30)*(41-34),-(1/30)*(41-34),-(1/30)*(41-56),-(1/30)*(41-56)],[-150 50 50 -150],[0.5 0.5 0.5],'edgecolor','none','FaceAlpha',0.4)

    for uidx=1:length(Pairs)
        uid = Pairs(uidx);

        hold on
        plot(-(1/30)*(41-(1:size(ProjectedWaveform,1))),squeeze(ProjectedWaveform(:,uid,1)),'color',cols(uidx,:))
    end
    hold on
    axis square
    xlabel('Time (ms)')
    ylabel('\muV')
    %     ylim([])
    makepretty
    title('Average Waveforms')

    %     subplot(1,3,2)
    %     for uidx=1:length(Pairs)
    %         uid = Pairs(uidx);
    %         channelpos = Allchannelpos{recsesGood(uid)};
    %         % Load raw data
    %         try
    %             spikeMap = readNPY(fullfile(param.KSDir{recsesGood(uid)},'RawWaveforms',['Unit' num2str(OriginalClusterIDs(Good_Idx(uid))) '_RawSpikes.npy'])); %0-indexed
    %         catch
    %             keyboard
    %         end
    %         % Detrending
    %         spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
    %         spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
    %         spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order
    %         %Load channels
    %         ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    %         Locs = channelpos(ChanIdx,:);
    %
    %         scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled')
    %         hold on
    %         scatter(ProjectedLocation(1,uid,1),ProjectedLocation(2,uid,1),20,cols(uidx,:),'filled')
    %         plot(ProjectedLocation(1,uid,1)+0.1*[1:size(ProjectedWaveform,1)],ProjectedLocation(2,uid,1)+0.1*ProjectedWaveform(:,uid,1),'color',cols(uidx,:))
    %     end
    %
    %     xlabel('Xpos (\mum)')
    %     ylabel('')
    %     set(gca,'xlim',xlims,'ylim',ylims,'YTickLabel',[])
    %
    %
    %     axis square
    %     title('Average waveforms')
    %     makepretty

    subplot(1,3,1)
    for uidx=1:length(Pairs)
        uid = Pairs(uidx);
        channelpos = Allchannelpos{recsesGood(uid)};
        % Load raw data
        try
            spikeMap = readNPY(fullfile(param.KSDir{recsesGood(uid)},'RawWaveforms',['Unit' num2str(OriginalClusterIDs(Good_Idx(uid))) '_RawSpikes.npy'])); %0-indexed
        catch
            keyboard
        end
        % Detrending
        spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
        spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
        spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order
        %Load channels
        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius.*2.5); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);

        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled','marker','s')
        hold on
        for chanid = 1:length(ChanIdx)
            plot(Locs(chanid,1)+0.1*[1:size(ProjectedWaveform,1)],Locs(chanid,2)+0.1*spikeMap(:,ChanIdx(chanid),1),'color',cols(uidx,:))
        end
    end

    xlabel('Xpos (\mum)')
    ylabel('Ypos (\mum)')
    ylims = [min(Locs(:,2))-15 max(Locs(:,2))+15];

    xlims = [min(ProjectedLocation(2,Pairs,1))-diff(ylims)/2 min(ProjectedLocation(2,Pairs,1))+diff(ylims)/2];
    set(gca,'xlim',xlims,'ylim',ylims)
    axis square
    title('all waveforms')
    makepretty

    figure('name','Probeview')
    scatter(Allchannelpos{1}(:,1),Allchannelpos{1}(:,2),20,[0.5 0.5 0.5],'filled','marker','s')

    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(Pairs(1),1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);
    hold on
    scatter(Locs(:,1),Locs(:,2),20,[1 0 0],'filled','marker','s')


    %
    %     figure
    %     subplot(2,1,1)
    %     for uid = 2:length(Pairs)
    %         tmp = sqrt(nansum((squeeze([ProjectedLocationPerTP(1,Pairs(1),param.waveidx,1),ProjectedLocationPerTP(2,Pairs(1),param.waveidx,1)])-squeeze([ProjectedLocationPerTP(1,Pairs(uid),param.waveidx,1),ProjectedLocationPerTP(2,Pairs(uid),param.waveidx,1)])).^2,1));
    %         tmp(tmp==0)=nan; % =nan
    %         plot(param.waveidx,tmp,'|','color',cols(uid,:))
    %         hold on
    %     end
    %     ylabel('\Deltad_i_j')
    %     makepretty
    %
    %
    %     subplot(2,1,2)
    %     for uid = 2:length(Pairs)
    %         plot(param.waveidx(2:end),squeeze(AngleSubtraction(Pairs(1),:,Pairs(uid))),'|','color',cols(uid,:))
    %         hold on
    %     end
    %     ylabel('\Deltad_i_j')
    %     makepretty


    %% Similarity scores for the two pairs
    figure('name','Similarity scores')
    subplot(1,2,1)
    bar(cat(1,squeeze(Predictors(Pairs(1),Pairs(3),:)),squeeze(TotalScore(Pairs(1),Pairs(3))./6)),'FaceColor',cols(3,:),'EdgeColor','none')
    set(gca,'XTick',1:size(Predictors,3)+1,'XTickLabel',{'C','W','V','D','A',char(920),'T'},'YAxisLocation','right','YTickLabelRotation',90,'XTickLabelRotation',90)
    ylabel('Normalized Score')
    makepretty

    subplot(1,2,2)
    bar(cat(1,squeeze(Predictors(Pairs(1),Pairs(2),:)),squeeze(TotalScore(Pairs(1),Pairs(2))./6)),'FaceColor',cols(2,:),'EdgeColor','none')
    set(gca,'XTick',1:size(Predictors,3)+1,'XTickLabel',{'C','W','V','D','A',char(920),'T'},'YAxisLocation','right','YTickLabelRotation',90,'XTickLabelRotation',90)
    ylabel('Normalized Score')
    makepretty





end
