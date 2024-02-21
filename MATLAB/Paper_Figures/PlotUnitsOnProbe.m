function PlotUnitsOnProbe(clusinfo,UMparam,UniqueIDConversion,WaveformInfo,AddDriftBack)
if nargin<5
    AddDriftBack = 1;
end

[~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));

% Draw findings on probe
neuroncols = jet(length(id1))+rand(length(id1),3)*2-1; % Colours per unit
neuroncols(neuroncols<0)=0;
neuroncols(neuroncols>1)=1;
neuroncols = datasample(neuroncols,length(id1),1,'replace',false);
% give same units same colour
neuroncols = neuroncols(id2,:);
id2ori = id2;
id1ori = id1;
neuroncolsori = neuroncols;
ExampleFig = figure('name',['Projection locations example units, DriftCorrected = ' num2str(~AddDriftBack)]');

% Extract DeltaDays
try
    days = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    days = cell2mat(days) - days{1};
catch ME
    disp('Can''t read in days')
end

for modethis = 1:2
    neuroncols = neuroncolsori;

    if modethis==2
        SaveSplitUnits = false(length(id2),1);

        Conservative = 1;
        [~,id1,id2] = unique(UniqueIDConversion.UniqueIDConservative(UniqueIDConversion.GoodID==1));
        % 
        splitUnit = find(arrayfun(@(X) length(unique(id2(id2ori==X))),id1ori)>1); % also unique in this one or two different ones?
        for Xid = 1:length(splitUnit)
            tmpselect = id2(id2ori==id1ori(splitUnit(Xid))); % Find new IDs of other pairs
            [tmpselect, iid1, iid2] = unique(tmpselect);
            [maxn,maxid] = max(arrayfun(@(X) sum(iid2==X),iid1)); % Keep the one which happens most the same colour
            tmpselect(maxid)=[];
            for newid = 1:length(tmpselect)
                SaveSplitUnits(id2==tmpselect(newid)) = true;
                neuroncols(id2==tmpselect(newid),:)=repmat(datasample(neuroncols,1,1),sum(id2==tmpselect(newid)),1); % Replace colour for second set of unique neurons
            end
        end
    else
        Conservative = 0;
        [~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));
    end
    % Units that only occur once:
    neuroncols(logical(arrayfun(@(X) sum(id2==X)==1,id2)),:) = repmat([0.5 0.5 0.5],sum(arrayfun(@(X) sum(id2==X)==1,id2)),1); % Gray if no match is found
    drift = [0 0 0];

    % Extract good rec ses id
    recsesGood = clusinfo.RecSesID(logical(clusinfo.Good_ID));
    nrec = unique(recsesGood);
    figure('name',['Projection locations all units, conservative = ' num2str(Conservative) ', DriftCorrected = ' num2str(~AddDriftBack)])
    for recid = 1:length(nrec)

        if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))
            drift = [0 0 0] +  + UMparam.drift(recid-1,:,1); % Initial drift only, is not cumulative across days
        end
        scatter3(UMparam.Coordinates{nrec(recid)}(:,1),UMparam.Coordinates{nrec(recid)}(:,2),UMparam.Coordinates{nrec(recid)}(:,3),30,[0 0 0],'square','filled')
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
            drift = [0 0 0] + UMparam.drift(recid-1,:,1); % Initial drift only, is not cumulative across days
        end
        subplot(ceil(sqrt(length(nrec))),round(sqrt(length(nrec))),recid)
        scatter3(UMparam.Coordinates{nrec(recid)}(:,1)-drift(1),UMparam.Coordinates{nrec(recid)}(:,2)-drift(2),UMparam.Coordinates{nrec(recid)}(:,3)-drift(3),30,[0 0 0],'square','filled')
        hold on
        scatter3(nanmean(WaveformInfo.ProjectedLocation(1,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid),:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid),:),3),30,neuroncols(recsesGood==nrec(recid),:),'filled')
        makepretty
        offsetAxes
        title(['Recording ' num2str(recid)])
    end
    linkaxes


    %% Show example neurons that was tracked across most recordings
    if modethis==1
        nRecPerUnit = arrayfun(@(X) sum(id2==X),id1);
        TrackedNeuronPop = id1(find(nRecPerUnit>0.5*length(nrec)));
        TrackedNeuronPop = ismember(id2,TrackedNeuronPop);
    end

    figure(ExampleFig)
    driftprobe = [0 0 0];

    for recid = 1:length(nrec)
        if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))
            driftprobe = [0 0 0] + UMparam.drift(recid-1,:,1); % Initial drift only, is not cumulative across days
        end

        subplot(length(nrec),2,(recid-1)*2+modethis)
        scatter(UMparam.Coordinates{nrec(recid)}(:,2)-driftprobe(2),UMparam.Coordinates{nrec(recid)}(:,3)-driftprobe(3),30,[0 0 0],'square','filled');
        hold on
        if modethis == 2
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & ~SaveSplitUnits,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnits,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnits,:),'filled');
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & SaveSplitUnits,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnits,:),3),70,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnits,:),'filled');
        else            
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop,:),'filled');
        end


        offsetAxes
        makepretty
        if recid==1
            title(['Conservative = ' num2str(Conservative)])
        end
        if recid < length(nrec)
            set(gca,'XTickLabel',[])
        end
        if modethis == 2
            set(gca,'YTickLabel',[])
        else
            if exist('days','var')
                ylabel(['\DeltaDays ' num2str(days(recid))])
            else
                ylabel(['R=' num2str(recid)])
            end
        end
        % axis equal
        % set(gca,'DataAspectRatio',[1 2000 1]);


    end

    linkaxes

end