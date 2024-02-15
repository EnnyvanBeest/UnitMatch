function PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo,Conservative,AddDriftBack)
if nargin<6
    AddDriftBack = 1;
end
if nargin<5
    Conservative = 1;
end
if Conservative
    [~,id1,id2] = unique(UniqueIDConversion.UniqueIDConservative(UniqueIDConversion.GoodID==1));
else
    [~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));
end

% Draw findings on probe
neuroncols = jet(length(id1)); % Colours per unit
neuroncols = datasample(neuroncols,length(id1),1,'replace',false); 
% give same units same colour
neuroncols = neuroncols(id2,:);

% Units that only occur once:
neuroncols(logical(arrayfun(@(X) sum(id2==X)==1,id2)),:) = repmat([0.5 0.5 0.5],sum(arrayfun(@(X) sum(id2==X)==1,id2)),1); % Gray if no match is found
drift = [0 0 0];

% Extract good rec ses id
recsesGood = clusinfo.RecSesID(logical(clusinfo.Good_ID));
nrec = unique(recsesGood);
figure('name',['Projection locations all units, conservative = ' num2str(Conservative) ', DriftCorrected = ' num2str(~AddDriftBack)])
for recid = 1:length(nrec)

    if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))       
        drift = drift + UMparam.drift(recid-1,:,1); % Initial drift only, is cumulative across days
    end
    scatter3(UMparam.Coordinates{nrec(recid)}(:,1),UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'filled','Marker','s')
    hold on
    scatter3(nanmean(WaveformInfo.ProjectedLocation(1,recsesGood==nrec(recid),:),3)+drift(1),nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid),:),3)+drift(2),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid),:),3)+drift(3),30,neuroncols(recsesGood==nrec(recid),:),'filled')
end
makepretty
offsetAxes

drift = [0 0 0];
% In separate plots
figure('name',['Projection locations all units, conservative = ' num2str(Conservative) ', DriftCorrected = ' num2str(~AddDriftBack)]')
for recid = 1:length(nrec)
   if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))    
        drift = drift + UMparam.drift(recid-1,:,1); % Initial drift only, is cumulative across days
    end
    subplot(ceil(sqrt(length(nrec))),round(sqrt(length(nrec))),recid)
    scatter3(UMparam.Coordinates{nrec(recid)}(:,1),UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'filled','Marker','s')
    hold on
    scatter3(nanmean(WaveformInfo.ProjectedLocation(1,recsesGood==nrec(recid),:),3)+drift(1),nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid),:),3)+drift(2),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid),:),3)+drift(3),30,neuroncols(recsesGood==nrec(recid),:),'filled')
    makepretty
offsetAxes
title(['Recording ' num2str(recid)])
end
linkaxes


%% Show example neurons that was tracked across most recordings
nRecPerUnit = arrayfun(@(X) sum(id2==X),id1);
TrackedNeuronPop = id1(find(nRecPerUnit>0.3*length(nrec)));

drift = [0 0 0];

figure('name',['Projection locations example units, conservative = ' num2str(Conservative) ', DriftCorrected = ' num2str(~AddDriftBack)]')
for recid = 1:length(nrec)
    if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))
        drift = drift + UMparam.drift(recid-1,:,1); % Initial drift only, is cumulative across days
    end

    subplot(length(nrec),1,recid)
    scatter(UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'filled','Marker','s')
    hold on
    scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & ismember(id2,TrackedNeuronPop),:),3)+drift(2),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& ismember(id2,TrackedNeuronPop),:),3)+drift(3),30,neuroncols(recsesGood==nrec(recid)& ismember(id2,TrackedNeuronPop),:),'filled')
    makepretty
    offsetAxes
    title(['Recording ' num2str(recid)])

end
linkaxes