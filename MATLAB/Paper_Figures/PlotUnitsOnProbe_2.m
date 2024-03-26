function PlotUnitsOnProbe_2(clusinfo,UMparam,UniqueIDConversion,WaveformInfo,AddDriftBack)
if nargin<5
    AddDriftBack = 0;
end

PlotAll = 0;
PlotMaxRecSes = 7; 

[~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));

% Draw findings on probe
neuroncols = distinguishable_colors(length(id1))-rand(length(id1),3).*0.25; % Colours per unit
neuroncols(neuroncols<0)=0;
neuroncols(neuroncols>1)=1;
neuroncols = datasample(neuroncols,length(id1),1,'replace',false);
% give same units same colour
neuroncols = neuroncols(id2,:);
id2ori = id2;
id1ori = id1;
neuroncolsori = neuroncols;
ExampleFig = figure('name',['Projection locations example units, DriftCorrected = ' num2str(~AddDriftBack)]');

recsesGood = clusinfo.RecSesID(logical(clusinfo.Good_ID));
nrec = unique(recsesGood);
if isempty(PlotMaxRecSes)
    PlotMaxRecSes = nrec;
end
SesIdx = [4, 6, 9, 12, 15, 18, 20]; 
% SesIdx = sort(datasample(nrec,PlotMaxRecSes,'Replace',false));
takethesenot = ~ismember(recsesGood,SesIdx);
WaveformInfo.ProjectedLocation(:,takethesenot,:) = [];
recsesGood(takethesenot) = [];
id2(takethesenot) = [];
nrec = unique(recsesGood);
nRecPerUnit = arrayfun(@(X) sum(id2==X),id2);
UnitstoConsider = id2(recsesGood==median(SesIdx) & nRecPerUnit>4);% All neurons tracked in median session
% TrackedNeuronPop = ismember(id2,UnitstoConsider);



% Extract DeltaDays
try
    if isunix
       % Adjusted regexp pattern for Unix-style paths with forward slashes
        datePattern = '/\d+-\d+-\d+/';
        
        % Extract 'folder' field values into a cell array
        folderPaths = cellfun(@(s) s.folder, UMparam.RawDataPaths, 'UniformOutput', false);
        
        % Use the adjusted regexp to find date strings in each folder path
        dateStrings = cellfun(@(path) regexp(path, datePattern, 'match'), folderPaths, 'UniformOutput', false);
        
        % Flatten the nested cell arrays produced by regexp (if necessary)
        dateStrings = cellfun(@(x) x{1}, dateStrings, 'UniformOutput', false);
        
        % Convert date strings to datenum format
        days = cellfun(@(dateStr) datenum(dateStr), dateStrings, 'UniformOutput', false);
    else
        days = cellfun(@(y) datenum(y), cellfun(@(x) regexp(x.folder,'\\\d*-\d*-\d*\\','match'), UMparam.RawDataPaths, 'uni', 0), 'uni', 0);
    end
    days = cell2mat(days) - days{1};



catch ME
    disp('Can''t read in days')
end

Modes = {'Liberal','Intermediate','Conservative'};
for modethis = 1:3
    neuroncols = neuroncolsori;

    if modethis==3 % Conservative

        [~,id1,id2] = unique(UniqueIDConversion.UniqueIDConservative(UniqueIDConversion.GoodID==1));
        SaveSplitUnitsCons = false(length(id2),1);

        % 
        splitUnit = find(arrayfun(@(X) length(unique(id2(id2ori==X))),id1ori)>1); % also unique in this one or two different ones?
        for Xid = 1:length(splitUnit)
            tmpselect = id2(id2ori==id1ori(splitUnit(Xid))); % Find new IDs of other pairs
            [tmpselect, iid1, iid2] = unique(tmpselect);
            [maxn,maxid] = max(arrayfun(@(X) sum(iid2==X),iid1)); % Keep the one which happens most the same colour
            tmpselect(maxid)=[];
            for newid = 1:length(tmpselect)
                SaveSplitUnitsCons(id2==tmpselect(newid)) = true;
                neuroncols(id2==tmpselect(newid),:)=repmat(datasample(neuroncols,1,1),sum(id2==tmpselect(newid)),1); % Replace colour for second set of unique neurons
            end
        end
    elseif modethis == 2 % Intermediate

        [~,id1,id2] = unique(UniqueIDConversion.UniqueID(UniqueIDConversion.GoodID==1));
        SaveSplitUnitsInterm = false(length(id2),1);

        %
        splitUnit = find(arrayfun(@(X) length(unique(id2(id2ori==X))),id1ori)>1); % also unique in this one or two different ones?
        for Xid = 1:length(splitUnit)
            tmpselect = id2(id2ori==id1ori(splitUnit(Xid))); % Find new IDs of other pairs
            [tmpselect, iid1, iid2] = unique(tmpselect);
            [maxn,maxid] = max(arrayfun(@(X) sum(iid2==X),iid1)); % Keep the one which happens most the same colour
            tmpselect(maxid)=[];
            for newid = 1:length(tmpselect)
                SaveSplitUnitsInterm(id2==tmpselect(newid)) = true;
                neuroncols(id2==tmpselect(newid),:)=repmat(datasample(neuroncols,1,1),sum(id2==tmpselect(newid)),1); % Replace colour for second set of unique neurons
            end
        end

    else

        [~,id1,id2] = unique(UniqueIDConversion.UniqueIDLiberal(UniqueIDConversion.GoodID==1));
    end
    % Units that only occur once:
    neuroncols(logical(arrayfun(@(X) sum(id2==X)==1,id2)),:) = repmat([0.5 0.5 0.5],sum(arrayfun(@(X) sum(id2==X)==1,id2)),1); % Gray if no match is found
    drift = [0 0 0];

    % Extract good rec ses id
    recsesGood = clusinfo.RecSesID(logical(clusinfo.Good_ID));
    nrec = unique(recsesGood);
    if ~isempty(PlotMaxRecSes)
        if modethis == 2
            SaveSplitUnitsInterm(takethesenot) = [];
        elseif modethis == 3
            SaveSplitUnitsCons(takethesenot) = [];
        end

        neuroncols(takethesenot,:) = [];
        id2(takethesenot) = [];
        recsesGood(takethesenot) = [];
        nrec = unique(recsesGood);
        if modethis == 1
            days = days(nrec);
        end
    end
    UnitstoConsider = id2(recsesGood==median(SesIdx) & nRecPerUnit>4);% All neurons tracked in median session

    TrackedNeuronPop = ismember(id2,UnitstoConsider);


    %% Show example neurons that was tracked across most recordings
    figure(ExampleFig)
    driftprobe = [0 0 0];
    for recid = 1:length(nrec)
        if AddDriftBack && recid > 1 && ~ any(isnan(UMparam.drift(recid-1,:,1)))
            driftprobe = [0 0 0] + UMparam.drift(recid-1,:,1); % Initial drift only, is not cumulative across days
        end

        subplot(length(nrec),3,(recid-1)*3+modethis)
        scatter(UMparam.Coordinates{nrec(recid)}(:,2)-driftprobe(2),UMparam.Coordinates{nrec(recid)}(:,3)-driftprobe(3),10,[0.2 0.2 0.2],'square','filled');
        hold on
        if modethis == 3
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & ~SaveSplitUnitsCons,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnitsCons,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnitsCons,:),'filled');
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & SaveSplitUnitsCons,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnitsCons,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnitsCons,:),'filled');
        elseif modethis == 2   
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & ~SaveSplitUnitsInterm,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnitsInterm,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & ~SaveSplitUnitsInterm,:),'filled');
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop & SaveSplitUnitsInterm,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnitsInterm,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop & SaveSplitUnitsInterm,:),'filled');
         else
            scatter(nanmean(WaveformInfo.ProjectedLocation(2,recsesGood==nrec(recid) & TrackedNeuronPop,:),3),nanmean(WaveformInfo.ProjectedLocation(3,recsesGood==nrec(recid)& TrackedNeuronPop,:),3),30,neuroncols(recsesGood==nrec(recid)& TrackedNeuronPop,:),'filled');
        end


        offsetAxes
        makepretty
        if recid==1
            title([Modes{modethis}])
        end
        if recid < length(nrec)
            set(gca,'XTickLabel',[])
        end
        if modethis >= 2
            set(gca,'YTickLabel',[])
        else
            if exist('days','var')
                ylabel(['Day ' num2str(days(recid))])
            else
                ylabel(['R=' num2str(recid)])
            end
        end
        % axis equal
        % set(gca,'DataAspectRatio',[1 2000 1]);


    end

    linkaxes

end


