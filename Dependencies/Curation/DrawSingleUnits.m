%% Create figures for manual curation of cells
KSPath = 'H:\PyKSOutput\EB001\2021-02-24\Probe0';


% Find params.py
spikeStruct = loadParamsPy(fullfile(KSPath,'params.py'));
rawD = spikeStruct.dat_path;
rawD = rawD(strfind(rawD,'"')+1:end);
rawD = rawD(1:strfind(rawD,'"')-1);
tmpdr = rawD;
rawD = dir(rawD);
if isempty(rawD)
    rawD = dir(strrep(tmpdr,'bin','cbin'));
end


%% Load Spike Data
sp = loadKSdir(fullfile(KSPath)); % Load Spikes with PCs
[sp.spikeAmps, sp.spikeDepths, sp.templateDepths, sp.templateXpos, sp.tempAmps, sp.tempsUnW, sp.templateDuration, sp.waveforms] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.xcoords, sp.spikeTemplates, sp.tempScalingAmps); %from the spikes toolbox

%% Remove noise; spikes across all channels'
sp = RemoveNoiseAmplitudeBased(sp);

%% Channel data
myClusFile = dir(fullfile(KSPath,'channel_map.npy'));
channelmaptmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

myClusFile = dir(fullfile(KSPath,'channel_positions.npy'));
channelpostmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
if length(channelmaptmp)<length(channelpostmp)
    channelmaptmp(end+1:length(channelpostmp))=length(channelmaptmp):length(channelpostmp)-1;
end
channelpostmpconv = ChannelIMROConversion(rawD(1).folder,0); % For conversion when not automatically done

%% Load Cluster Info
myClusFile = dir(fullfile(KSPath,'cluster_info.tsv')); % If you did phy (manual curation) we will find this one... We can trust you, right?
if isempty(myClusFile)
    disp('This data is not curated with phy! Hopefully you''re using automated quality metrics to find good units!')
    curratedflag=0;
    myClusFile = dir(fullfile(KSPath,'cluster_group.tsv'));
    clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
    % Convert sp data to correct cluster according to phy (clu and
    % template are not necessarily the same after phy)
    [clusinfo,sp] = ConvertTemplatesAfterPhy(clusinfo,sp);

    cluster_id = clusinfo.cluster_id;
    Label = char(1,length(clusinfo.cluster_id));
    KSLabelfile = tdfread(fullfile(KSPath,'cluster_KSLabel.tsv'));
    Label(ismember(clusinfo.cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,clusinfo.cluster_id));
    totSpkNum = histc(sp.clu,sp.cids);
    Good_IDtmp = ismember(Label,'g')' & totSpkNum'>300; %Identify good clusters

    % Find depth and channel
    depthtmp = nan(length(clusinfo.cluster_id),1);
    xtmp = nan(length(clusinfo.cluster_id),1);
    channeltmp = nan(length(clusinfo.cluster_id),1);
    for clusid=1:length(depthtmp)
        % Depth according to PyKS2 output
        depthtmp(clusid)=round(sp.templateDepths(clusid));%round(nanmean(sp.spikeDepths(find(sp.clu==clusid-1))));
        xtmp(clusid)=sp.templateXpos(clusid);
        [~,minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1,channelpostmp(X,:),[xtmp(clusid),depthtmp(clusid)]),'euclidean'),1:size(channelpostmp,1),'UniformOutput',0)));
        try
            channeltmp(clusid) = channelmaptmp(minidx);
            depthtmp(clusid) = channelpostmpconv(minidx,2);
            xtmp(clusid) = channelpostmpconv(minidx,1);
        catch
            channeltmp(clusid)=channelmaptmp(minidx-1);
            depthtmp(clusid) = channelpostmpconv(minidx-1,2);
            xtmp(clusid) = channelpostmpconv(minidx-1,1);
        end
        sp.spikeDepths(ismember(sp.clu,cluster_id(clusid))) = depthtmp(clusid);

    end
    depth = depthtmp;
    channel = channeltmp;

else
    disp('You did manual curation. You champion. If you have not enough time, maybe consider some automated algorithm...')
    CurationDone = 1;
    save(fullfile(KSPath,'CuratedResults.mat'),'CurationDone')
    clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
    % Convert sp data to correct cluster according to phy (clu and
    % template are not necessarily the same after phy)
    [clusinfo,sp] = ConvertTemplatesAfterPhy(clusinfo,sp);

    curratedflag=1;
    if isfield(clusinfo,'id')
        clusidtmp = clusinfo.id;
        cluster_id = [cluster_id,clusinfo.id];
    elseif isfield(clusinfo,'cluster_id')
        clusidtmp = clusinfo.cluster_id;

        cluster_id=[cluster_id,clusinfo.cluster_id];
    else
        keyboard
        disp('Someone thought it was nice to change the name again...')
    end
    KSLabel = clusinfo.KSLabel;
    Label = [Label,clusinfo.group]; % You want the group, not the KSLABEL!
    % Find depth and channel
    depthtmp = nan(length(clusidtmp),1);
    xtmp = nan(length(clusidtmp),1);
    channeltmp = nan(length(clusidtmp),1);
    for clusid=1:length(depthtmp)
        % Depth according to PyKS2 output
        depthtmp(clusid)=round(sp.templateDepths(clusid));%round(nanmean(sp.spikeDepths(find(sp.clu==clusid-1))));
        xtmp(clusid)=sp.templateXpos(clusid);
        [~,minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1,channelpostmp(X,:),[xtmp(clusid),depthtmp(clusid)]),'euclidean'),1:size(channelpostmp,1),'UniformOutput',0)));
        try
            channeltmp(clusid) = channelmaptmp(minidx);
            depthtmp(clusid) = channelpostmpconv(minidx,2);
            xtmp(clusid) = channelpostmpconv(minidx,1);
        catch
            channeltmp(clusid)=channelmaptmp(minidx-1);
            depthtmp(clusid) = channelpostmpconv(minidx-1,2);
            xtmp(clusid) = channelpostmpconv(minidx-1,1);
        end
        sp.spikeDepths(ismember(sp.clu,clusidtmp(clusid))) = depthtmp(clusid);

    end

    depth = depthtmp;

    channel = clusinfo.ch;
    Good_IDtmp = ismember(cellstr(clusinfo.group),'good');
    channeltmp = clusinfo.ch;
end
channelpos = channelpostmpconv;
%% Draw them
if ~exist(fullfile(KSPath,'IndividualUnitFigures'))
    mkdir(fullfile(KSPath,'IndividualUnitFigures'))
end

cols = jet(2);
Good_IDx = find(Good_IDtmp);
for gid = 1:length(Good_IDx)
    tmpfig = figure('name',['Cluster ' num2str(cluster_id(Good_IDx(gid)))],'visible','off')

    % Show average waveform
    spikeMap = readNPY(fullfile(KSPath,'RawWaveforms',['Unit' num2str(cluster_id(Good_IDx(gid))) '_RawSpikes.npy']));

    % Detrending
    spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
    spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
    spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order

    % Extract channel positions that are relevant and extract mean location
    [~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
    ChanIdx = find(cell2mat(arrayfun(@(Y) vecnorm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<50); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);

    % waveform per channel
    subplot(2,3,[1])
    hold on
    scatter(Locs(:,1)*10,Locs(:,2)*20,40,[0.5 0.5 0.5],'filled') % Indicate sites
    for cv = 1:2
        for id = 1:size(Locs,1)
            plot(Locs(id,1)*10+[1:size(spikeMap,1)],Locs(id,2)*20+spikeMap(:,ChanIdx(id),cv),'-','color',cols(cv,:),'LineWidth',1)
        end
    end
    xlabel('x-pos')
    ylabel('Depth on probe')
    title('Average waveform on different sites')
    makepretty
    
    % Weighted Average waveform
    for cv = 1:2
        mu = sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1).*Locs',2)./sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1),2);
        % Use this waveform - weighted average across channels:
        Distance2MaxProj = vecnorm(Locs-mu',2,2);
        weight = (50-Distance2MaxProj)./50;
        ProjectedWaveform = nansum(spikeMap(:,ChanIdx,cv).*repmat(weight,1,size(spikeMap,1))',2)./sum(weight);

        subplot(2,3,2)
        hold on
        plot(ProjectedWaveform,'color',cols(cv,:))

        % Normalized waveform
        ProjectedWaveform = (ProjectedWaveform-nanmin(ProjectedWaveform(:)))./(nanmax(ProjectedWaveform(:))-nanmin(ProjectedWaveform(:)));
        subplot(2,3,3)
        hold on
        plot(ProjectedWaveform,'color',cols(cv,:))
    end
    subplot(2,3,2)
    ylabel('\muV')
    xlabel('time')
    title('Average waveform')
    makepretty

    subplot(2,3,3)
    ylabel('\muV')
    xlabel('time')
    title('Normalized waveform')
    makepretty

    % Scatter spikes
    idx1=find(sp.spikeTemplates == cluster_id(Good_IDx(gid)));
    if ~isempty(idx1)
        for cv = 1:2
            if cv==1
                idx1 = idx1(1:floor(length(idx1)/2));
            else
                idx1 = idx1(ceil(length(idx1)/2):end);
            end

            subplot(2,3,[4,5])
            hold on

            scatter(sp.st(idx1)./60,sp.spikeAmps(idx1),4,cols(cv,:),'filled')

            xlims = get(gca,'xlim');
            % Other axis
            [h1,edges,binsz]=histcounts(sp.spikeAmps(idx1));
            %Normalize between 0 and 1
            h1 = ((h1-nanmin(h1))./(nanmax(h1)-nanmin(h1)))*10+max(sp.st./60);
            plot(h1,edges(1:end-1),'-','color',cols(cv,:));

            % compute ACG
            [ccg, t] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
                ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', 0.001, 'duration', 0.5, 'norm', 'rate'); %function
            ACG = ccg(:, 1, 1);
            xlabel('Time (s)')
            ylabel('Amplitude')
            title('Spikes & amplitude distribution')
            makepretty

            subplot(2,3,[6])
            hold on
            plot(t,ACG,'color',cols(cv,:));
            title(['AutoCorrelogram'])
            xlim([-0.1 0.1])
        end
    end
    makepretty

    set(tmpfig,'units','normalized','outerposition',[0 0 1 1])


    fname = ['Cluster ' num2str(cluster_id(Good_IDx(gid)))];
    saveas(tmpfig,fullfile(KSPath,'IndividualUnitFigures',[fname '.fig']))
    saveas(tmpfig,fullfile(KSPath,'IndividualUnitFigures',[fname '.bmp']))
    delete(tmpfig)
    






 



    

end