function PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo)

% Draw findings on probe
[~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));

neuroncols = jet(length(id1)); % Colours per unit
neuroncols = datasample(neuroncols,length(id1),1,'replace',false); 
% give same units same colour
neuroncols = neuroncols(id2,:);

% Units that only occur once:
neuroncols(logical(arrayfun(@(X) sum(id2==X)==1,id2)),:) = repmat([0.5 0.5 0.5],sum(arrayfun(@(X) sum(id2==X)==1,id2)),1); % Gray if no match is found


% Extract good rec ses id
recsesGood = clusinfo.RecSesID(logical(clusinfo.Good_ID));
nrec = unique(recsesGood);
figure('name','Projection locations all units')
for recid = 1:length(nrec)
    scatter3(UMparam.Coordinates{nrec(recid)}(:,1),UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'filled','Marker','s')
    hold on
    scatter3(nanmean(WaveformInfo.ProjectedLocation(1,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid),:),3),30,neuroncols(recsesGood==nrec(recid),:),'filled')
end
makepretty
offsetAxes

% In separate plots
figure('name','Projection locations all units')
for recid = 1:length(nrec)
    subplot(ceil(sqrt(length(nrec))),round(sqrt(length(nrec))),recid)
    scatter3(UMparam.Coordinates{nrec(recid)}(:,1),UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'filled','Marker','s')
    hold on
    scatter3(nanmean(WaveformInfo.ProjectedLocation(1,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid),:),3),30,neuroncols(recsesGood==nrec(recid),:),'filled')
    makepretty
offsetAxes
title(['Recording ' num2str(recid)])
end
linkaxes