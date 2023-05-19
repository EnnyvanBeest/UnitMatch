function ExtractSimilarityMetrics(Scores2Include,AllWVBParameters,clusinfo,param,drawthis)

if nargin<5
    drawthis=1;
end
%% Extract fields and parameters
ExtractFields({AllWVBParameters})
waveidx = param.waveidx;
Allchannelpos = param.channelpos;
SaveDir = param.SaveDir;
maxdist = param.maxdist;

Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
GoodRecSesID = clusinfo.RecSesID(Good_Idx);
OriginalClusterIDs = clusinfo.cluster_id;

recsesAll = clusinfo.RecSesID;
recsesGood = recsesAll(Good_Idx);
DepthOnProbe = clusinfo.depth;
if length(DepthOnProbe) == length(recsesAll)/2
    DepthOnProbe = [DepthOnProbe DepthOnProbe]; %Stitched
end
DepthOnProbe = DepthOnProbe(Good_Idx);
[X,Y]=meshgrid(recsesAll(Good_Idx));
nclus = length(Good_Idx);
ndays = length(unique(recsesAll));
SessionSwitch = arrayfun(@(X) find(GoodRecSesID==X,1,'first'),1:ndays,'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch) nclus+1];

%% Compute Metrics
disp('Computing Metric similarity between pairs of units...')
timercounter = tic;
x1 = repmat(PeakTime(:,1),[1 numel(PeakTime(:,1))]);
x2 = repmat(PeakTime(:,2),[1 numel(PeakTime(:,2))]);
PeakTimeSim = abs(x1 - x2');
%Normalize between 0 and 1 (values that make sense after testing, without having outliers influence this)
PeakTimeSim =1-PeakTimeSim./quantile(PeakTimeSim(:),0.99);
PeakTimeSim(PeakTimeSim<0)=0;

if any(ismember(Scores2Include,'waveformTimePointSim'))
    % can't find much better for this one
    waveformTimePointSim = nan(nclus,nclus);
    for uid = 1:nclus
        for uid2 = 1:nclus
            waveformTimePointSim(uid,uid2) = sum(ismember(find(WaveIdx(uid,:,1)),find(WaveIdx(uid2,:,2))))./sum((WaveIdx(uid,:,1)));
        end
    end
end

if any(ismember(Scores2Include,'spatialdecaySim'))
    x1 = repmat(spatialdecay(:,1),[1 numel(spatialdecay(:,1))]);
    x2 = repmat(spatialdecay(:,2),[1 numel(spatialdecay(:,2))]);
    spatialdecaySim = abs(x1 - x2');
    % Make (more) normal
    spatialdecaySim = sqrt(spatialdecaySim);
    spatialdecaySim = 1-((spatialdecaySim-nanmin(spatialdecaySim(:)))./(quantile(spatialdecaySim(:),0.99)-nanmin(spatialdecaySim(:))));
    spatialdecaySim(spatialdecaySim<0)=0;
end

if any(ismember(Scores2Include,'AmplitudeSim'))
    % Ampitude difference
    x1 = repmat(Amplitude(:,1),[1 numel(Amplitude(:,1))]);
    x2 = repmat(Amplitude(:,2),[1 numel(Amplitude(:,2))]);
    AmplitudeSim = abs(x1 - x2');
    % Make (more) normal
    AmplitudeSim = sqrt(AmplitudeSim);
    AmplitudeSim = 1-((AmplitudeSim-nanmin(AmplitudeSim(:)))./(quantile(AmplitudeSim(:),.99)-nanmin(AmplitudeSim(:))));
    AmplitudeSim(AmplitudeSim<0)=0;
end

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
channelpos_AllCat = unique(cat(1,Allchannelpos{:}),'rows');

flag=0;
while flag<2
    timercounter = tic;
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
    LocDist = sqrt((ProjectedLocation(1,:,1)'-ProjectedLocation(1,:,2)).^2 + ...
        (ProjectedLocation(2,:,1)'-ProjectedLocation(2,:,2)).^2);
    LocDist = 1-((LocDist-nanmin(LocDist(:)))./(maxdist-nanmin(LocDist(:))));
    LocDist(LocDist<0)=0;

    disp('Computing location distances between pairs of units, per individual time point of the waveform...')
    % Difference in distance between centroids of two halfs of the recording
    x1 = repmat(squeeze(ProjectedLocationPerTP(:,:,waveidx,1)),[1 1 1 size(ProjectedLocationPerTP,2)]);
    x2 = permute(repmat(squeeze(ProjectedLocationPerTP(:,:,waveidx,2)),[1 1 1 size(ProjectedLocationPerTP,2)]),[1 4 3 2]);
    EuclDist = squeeze(sqrt(nansum((x1-x2).^2,1))); % Euclidean distance
    w = squeeze(isnan(abs(x1(1,:,:,:)-x2(1,:,:,:))));
    EuclDist(w) = nan;
    % Average location
    CentroidDist = squeeze(nanmean(EuclDist,2));%
    % Normalize each of them from 0 to 1, 1 being the 'best'
    % If distance > maxdist micron it will never be the same unit:
    CentroidDist = 1-((CentroidDist-nanmin(CentroidDist(:)))./(maxdist-nanmin(CentroidDist(:)))); %Average difference
    CentroidDist(CentroidDist<0)=0;
    CentroidDist(isnan(CentroidDist))=0;

    % Variance in error, corrected by average error. This captures whether
    % the trajectory is consistenly separate
    CentroidVar = squeeze(nanvar(EuclDist,[],2));%./nanmean(EuclDist,2)+nanmean(EuclDist,2));
    CentroidVar = sqrt(CentroidVar);
    CentroidVar = 1-((CentroidVar-nanmin(CentroidVar(:)))./(nanmax(CentroidVar(:))-nanmin(CentroidVar(:)))); %Average difference

    % @Célian suggestion: recenter to 0 first, then calculate the
    % difference in distance (to account for uncorrected drift)
    ProjectedLocationPerTPRecentered = permute(permute(ProjectedLocationPerTP,[1,2,4,3]) - ProjectedLocation,[1,2,4,3]);
    x1 = repmat(squeeze(ProjectedLocationPerTPRecentered(:,:,waveidx,1)),[1 1 1 size(ProjectedLocationPerTPRecentered,2)]);
    x2 = permute(repmat(squeeze(ProjectedLocationPerTPRecentered(:,:,waveidx,2)),[1 1 1 size(ProjectedLocationPerTPRecentered,2)]),[1 4 3 2]);
    EuclDist2 = squeeze(sqrt(nansum((x1-x2).^2,1))); % Euclidean distance
    w = squeeze(isnan(abs(x1(1,:,:,:)-x2(1,:,:,:))));
    EuclDist2(w) = nan;
    % Average location
    CentroidDistRecentered = squeeze(nanmean(EuclDist2,2));%
    CentroidDistRecentered = 1-(CentroidDistRecentered-nanmin(CentroidDistRecentered(:)))./(nanmax(CentroidDistRecentered(:))-nanmin(CentroidDistRecentered(:)));
    
    CentroidOverlord = (CentroidDistRecentered+CentroidVar)/2;

    disp('Computing location angle (direction) differences between pairs of units, per individual time point of the waveform...')
    x1 = ProjectedLocationPerTP(:,:,waveidx(2):waveidx(end),:);
    x2 = ProjectedLocationPerTP(:,:,waveidx(1):waveidx(end-1),:);
    % The distance traveled (Eucledian)
    TrajDist = sqrt(squeeze(nansum((x1-x2).^2,1)));
    % Difference in angle between two time points
    LocAngle = squeeze(atan(abs(x1(1,:,:,:)-x2(1,:,:,:))./abs(x1(2,:,:,:)-x2(2,:,:,:))));
    % Actually just taking the weighted sum of angles is better
    x1 = repmat(LocAngle(:,:,1),[1 1 nclus]);
    x2 = permute(repmat(LocAngle(:,:,2),[1 1 nclus]),[3 2 1]); %
    AngleSubtraction = abs(x1-x2);
    AngleSubtraction(isnan(abs(x1-x2))) = 0.5*pi; %punish points with nan
    TrajAngleSim = squeeze(nansum(AngleSubtraction,2)); % sum of angles
    TrajAngleSim = 1-((TrajAngleSim-nanmin(TrajAngleSim(:)))./(nanmax(TrajAngleSim(:))-nanmin(TrajAngleSim(:))));
    TrajAngleSim(TrajAngleSim<0 | isnan(TrajAngleSim))=0;

    % Continue distance traveled
    x1 = repmat(TrajDist(:,:,1),[1 1 nclus]);
    x2 = permute(repmat(TrajDist(:,:,2),[1 1 nclus]),[3 2 1]); %
    % Distance similarity (subtract for each pair of units)
    TrajDistCompared = abs(x1-x2);%
    TrajDistSim = squeeze(nansum(TrajDistCompared,2));
    TrajDistSim = sqrt(TrajDistSim); % Make more normal
    TrajDistSim = 1-((TrajDistSim-nanmin(TrajDistSim(:)))./(nanmax(TrajDistSim(:))-nanmin(TrajDistSim(:))));


    LocTrajectorySim = (TrajAngleSim+TrajDistSim)./2; % Trajectory Similarity is sum of distance + sum of angles
    LocTrajectorySim = (LocTrajectorySim-nanmin(LocTrajectorySim(:))./(nanmax(LocTrajectorySim(:))-nanmin(LocTrajectorySim(:))));
    %
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
    disp(['Extracting projected location took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

    % Average EuclDist
    EuclDist = squeeze(nanmean(EuclDist,2));
    % Plotting order (sort units based on distance)
    [~,SortingOrder] = arrayfun(@(X) sort(DepthOnProbe(SessionSwitch(X):SessionSwitch(X+1)-1)),1:ndays,'Uni',0);
    SortingOrder = arrayfun(@(X) SortingOrder{X}+SessionSwitch(X)-1,1:ndays,'Uni',0);
    SortingOrder = cat(2,SortingOrder{:});

    %% These are the parameters to include:
    if drawthis
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
        saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.fig'))
        saveas(gcf,fullfile(SaveDir,'TotalScoreComponents.bmp'))


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
    IncludeThesePairs = find(EuclDist<maxdist);

    disp('Computing total score...')
    timercounter = tic;
%     priorMatch = 1-((nclus+nclus.*sqrt(ndays-1))./length(IncludeThesePairs)); %Punish multiple days (unlikely to find as many matches after a few days)
    priorMatch = 1-((nclus+nclus.*sqrt(ndays-1))./length(IncludeThesePairs)); %Punish multiple days (unlikely to find as many matches after a few days)

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
            ThrsOpt = quantile(TotalScore(IncludeThesePairs),priorMatch); %Select best ones only later
            if ThrsOpt == max(TotalScore(IncludeThesePairs))
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
    for sid=1:length(Scores2Include)
        eval(['tmp = ' Scores2Include{sid} ';'])
        % Take centroid dist > maxdist out
        %         tmp(EuclDist>param.NeighbourDist)=nan;
        Predictors = cat(3,Predictors,tmp);
        TotalScore=TotalScore+tmp;
    end

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
    imagesc(TotalScore(SortingOrder,SortingOrder),[0 length(Scores2Include)]);
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
    ThrsOpt = quantile(TotalScore(IncludeThesePairs),priorMatch); %Select best ones only later
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


    subplot(2,2,3)
    tmp = TotalScore;
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
    hd = histcounts(diag(tmp),0:0.1:length(Scores2Include))./nclus;
    hnd = histcounts(tmp(~eye(size(tmp))),0:0.1:length(Scores2Include))./sum(~isnan(tmp(~eye(size(tmp)))));
    plot(0.05:0.1:length(Scores2Include)-0.05,hd,'-','color',[0.5 0.5 0.5]); hold on; plot(0.05:0.1:length(Scores2Include)-0.05,hnd,'k-')
    line([ThrsOpt ThrsOpt],get(gca,'ylim'),'LineStyle','--','color',[1 0 0])
    xlabel('TotalScore')
    ylabel('Proportion|Group')
    makepretty
    title('Within session cross-validation')
    legend('Within session diagonal',['Within session off-diagonal (<' num2str(param.NeighbourDist) ')'],'Threshold','Location', 'best')


    subplot(2,2,4)
    tmp = TotalScore;
    % Take centroid dist > maxdist out
    tmp(EuclDist>param.NeighbourDist)=nan;
    % Take within session out
    for did = 1:ndays       
        tmp(SessionSwitch(did):SessionSwitch(did+1),SessionSwitch(did):SessionSwitch(did+1))=nan;
    end
    hb = histcounts(tmp(~isnan(tmp(:))),0:0.1:length(Scores2Include))./sum(tmp(~isnan(tmp(:))));
    hd = histcounts(tmp(tmp(:)>ThrsOpt),0:0.1:length(Scores2Include))./sum(tmp(:)>ThrsOpt);
    hnd = histcounts(tmp(tmp(:)<=ThrsOpt),0:0.1:length(Scores2Include))./sum(tmp(:)<=ThrsOpt);
    hold on
    plot(0.05:0.1:length(Scores2Include)-0.05,hd,'-','color',[0 0.5 0]); hold on; plot(0.05:0.1:length(Scores2Include)-0.05,hnd,'k-'); 
    xlabel('TotalScore')
    ylabel('Proportion|Group')
    makepretty
    title(['Across sessions, EuclDist< ' num2str(param.NeighbourDist) 'um'])
    legend('T>threshold','T<threshold','Location', 'best')


    saveas(gcf,fullfile(SaveDir,'TotalScore.fig'))
    saveas(gcf,fullfile(SaveDir,'TotalScore.bmp'))
    % Find all pairs
    % first factor authentication: score above threshold
    % Take into account:
    label = TotalScore>ThrsOpt;
    [uid,uid2] = find(label);
    Pairs = cat(2,uid,uid2);
    Pairs = sortrows(Pairs);
    Pairs = unique(Pairs,'rows');
    Pairs(Pairs(:,1) == Pairs(:,2),:)=[];

    %% Get correlation matrics for fingerprint correlations
    Unit2Take = OriginalClusterIDs(Good_Idx);
    sessionCorrelationsAll = cell(1,ndays);
    for did = 1:ndays
        % Load sp for correct day
        if length(param.KSDir)>1
            tmp = matfile(fullfile(param.KSDir{did},'PreparedData.mat'));
        else %Stitched
            tmp = matfile(fullfile(param.KSDir{1},'PreparedData.mat'));
        end
        SessionCorrelations = tmp.SessionCorrelations;
        if length(tmp.SessionCorrelations)==ndays
            sessionCorrelationsAll{did} = SessionCorrelations{did};
        elseif length(tmp.SessionCorrelations)==1  %Normal situation
            if iscell(SessionCorrelations)
                sessionCorrelationsAll{did} = SessionCorrelations{1};
            else
                sessionCorrelationsAll{did} = SessionCorrelations;
            end
        else
            disp('This is a weird situation...')
            keyboard
        end
    end
  
    %% three ways to define candidate scores
    % Total score larger than threshold
    CandidatePairs = TotalScore>ThrsOpt;%

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

%% Assign to workspace
assignin('caller','TotalScore',TotalScore)
assignin('caller','Predictors',Predictors)
assignin('caller','drift',drift)
assignin('caller','ProjectedLocationPerTP',ProjectedLocationPerTP)
assignin('caller','ProjectedLocation',ProjectedLocation)
assignin('caller','sessionCorrelationsAll',sessionCorrelationsAll)
assignin('caller','Unit2Take',Unit2Take)
assignin('caller','EuclDist',EuclDist)
assignin('caller','SortingOrder',SortingOrder)

return
if 0 % THis can be used to look at some example projections
    %% Plot
    Pairs = [2,276,4] % Example
    cols =  jet(length(Pairs));

    figure
    subplot(5,1,1:3)
    for uidx=1:length(Pairs)
        uid = Pairs(uidx);
        channelpos = Allchannelpos{recsesGood(uid)};
        % Load raw data
        try
            spikeMap = readNPY(fullfile(param.KSDir{recsesGood(uid)},'RawWaveforms',['Unit' num2str(OriginalClusterIDs(uid)+1) '_RawSpikes.npy'])); %0-indexed to 1-indexed
        catch
            keyboard
        end
        % Detrending
        spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
        spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
        spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order
        %Load channels
        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);

        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled')
        hold on
        scatter(ProjectedLocation(1,uid,1),ProjectedLocation(2,uid,1),20,cols(uidx,:),'filled')

        takesamples = param.waveidx;
        takesamples = unique(takesamples(~isnan(takesamples)));
        h(1) = plot(squeeze(ProjectedLocationPerTP(1,uid,takesamples,1)),squeeze(ProjectedLocationPerTP(2,uid,takesamples,1)),'-','color',cols(uidx,:));
        scatter(squeeze(ProjectedLocationPerTP(1,uid,takesamples,1)),squeeze(ProjectedLocationPerTP(2,uid,takesamples,1)),30,takesamples,'filled')
    end
    colormap(hot)

    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    xlims = [ProjectedLocation(1,uid,1)-30 ProjectedLocation(1,uid,1)+30];
    ylims = [ProjectedLocation(2,uid,1)-30 ProjectedLocation(2,uid,1)+30];
    set(gca,'xlim',xlims,'ylim',ylims)
    axis square
    %     legend([h(1),h(2)],{['Unit ' num2str(uid)],['Unit ' num2str(uid2)]})
    hc= colorbar;
    try
        hc.Label.String = 'timesample';
    catch ME
        disp(ME)
        keyboard
    end
    makepretty

    subplot(5,1,4)
    for uid = 2:length(Pairs)
        tmp = sqrt(nansum((squeeze([ProjectedLocationPerTP(1,Pairs(1),param.waveidx,1),ProjectedLocationPerTP(2,Pairs(1),param.waveidx,1)])-squeeze([ProjectedLocationPerTP(1,Pairs(uid),param.waveidx,1),ProjectedLocationPerTP(2,Pairs(uid),param.waveidx,1)])).^2,1));
        tmp(tmp==0)=nan; % =nan
        plot(param.waveidx,tmp,'|','color',cols(uid,:))
        hold on
    end
    ylabel('\Deltad_i_j')
    makepretty


    subplot(5,1,5)
    for uid = 2:length(Pairs)
        plot(param.waveidx(2:end),squeeze(AngleSubtraction(Pairs(1),:,Pairs(uid))),'|','color',cols(uid,:)) 
        hold on
    end
    ylabel('\Deltad_i_j')
    makepretty

    figure
    for uidx=1:length(Pairs)
        uid = Pairs(uidx);
        
        hold on
        plot(squeeze(ProjectedWaveform(:,uid,1)),'color',cols(uidx,:))
     end


end
