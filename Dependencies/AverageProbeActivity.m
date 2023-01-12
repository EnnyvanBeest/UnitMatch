%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
for idx = 1:length(DataDir)
    DateOpt{idx} = cellfun(@(X) dir(fullfile(DataDir{idx},X,'*-*')),MiceOpt(DataDir2Use==idx),'UniformOutput',0);
end
DateOpt = cat(2,DateOpt{:});
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
RedoAfterClustering=1;
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering
PlotLFP = 0;
plotUnitActivity = 0;
GoodUnits = cell(1,length(MiceOpt));
for midx = length(MiceOpt)
    myKsDir = fullfile(LocalDir,MiceOpt{midx});
    if strcmp(RecordingType,'Chronic')  % These are my chronic mice, one dataset per mouse
        %% Loading data from kilosort/phy easily
        myKsDir = fullfile(LocalDir,MiceOpt{midx});
        subksdirs = dir(fullfile(myKsDir,'Chronic','Probe*')); %This changed because now I suddenly had 2 probes per recording
        if length(subksdirs)<1
            clear subksdirs
            subksdirs.folder = myKsDir; %Should be a struct array
            subksdirs.name = 'Probe0';
        end

        ProbeOpt = (unique({subksdirs(:).name}));
        for probeid = 1:length(ProbeOpt)
            %Saving directory
            thisprobe =  ProbeOpt{probeid}
            myKsDir = fullfile(LocalDir,MiceOpt{midx},'*',thisprobe);
            % Check for multiple subfolders?
            subsesopt = dir(myKsDir);
            subsesopt(~[subsesopt.isdir])=[];
            subsesopt(ismember({subsesopt(:).name},{'..','.phy'}))=[];

            if strcmp(RecordingType{midx},'Chronic')
                if ~RunPyKSChronic %MatchUnitsAcrossDays
                    disp('Unit matching in Matlab')
                    subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),{subsesopt(:).folder},'UniformOutput',0))) = []; %Use separate days and match units via matlab script
                else
                    disp('Using chronic pyks option')
                    subsesopt = subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),{subsesopt(:).folder},'UniformOutput',0))); %Use chronic output from pyks
                end
            end

            % Copy file and then run pykilosort on it
            channelposition = cell(1,length(subsesopt));
            for id = 1:length(subsesopt)
                tmpfile = dir(fullfile(subsesopt(id).folder,subsesopt(id).name,'channel_positions.npy'));
                if isempty(tmpfile)
                    continue
                end
                channelposition{id} = readNPY(fullfile(tmpfile.folder,tmpfile.name));
            end
            subsesopt(cellfun(@isempty,channelposition))=[];
            channelposition(cellfun(@isempty,channelposition))=[];
            subsesoptAll = subsesopt;
            ImroCount = 1;
            IMROTableOpt = {};
            ChannelPosOpt = {};
            subsesoptGroups = {};
            for id = 1:length(channelposition)
                id1 = find(cell2mat(cellfun(@(X) all(channelposition{id}(:)==X(:)),ChannelPosOpt,'Uni',0)));
                if ~isempty(id1)
                    subsesoptGroups{id1} = [subsesoptGroups{id1} id];
                else
                    IMROTableOpt = {IMROTableOpt{:} ['IMRO' num2str(ImroCount)]};
                    subsesoptGroups{ImroCount} = id;
                    ChannelPosOpt{ImroCount} = channelposition{id};
                    ImroCount = ImroCount+1;
                end
            end

            for IMROID = 1:length(subsesoptGroups)
                %% Create saving directory
                clear params
                params.loadPCs=true;
                params.thisdate=[];
                thisIMRO = IMROTableOpt{IMROID};
                thisdate = params.thisdate;
                if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,thisIMRO))
                    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,thisIMRO))
                end
                subsesopt = subsesoptAll(subsesoptGroups{IMROID});

                %% Get cluster information                
                PrepareClusInfo
                %% Save out
                GoodUnits{midx}{probeid}{IMROID} = clusinfo;

                %% Show Activity on probe
                depthbins = min(depth):100:max(depth);
                depthavg = min(depth)+50:100:max(depth)-50;
                spShDepth  = nan(length(depthbins)-1,length(ShankOpt));
                for shid = 1:length(ShankOpt)
                    % Clusters on this shank
                    parfor deptid=1:length(depthbins)-1
                        tmpclu = spikeCluster(ismember(spikeShank,ShankOpt(shid))& spikeDepths>=depthbins(deptid)&spikeDepths<=depthbins(deptid+1));

                        tmpst = spikeTimes(ismember(spikeShank,ShankOpt(shid))& spikeDepths>=depthbins(deptid)&spikeDepths<=depthbins(deptid+1));

                        spShDepth(deptid,shid) = nanmean(histcounts(tmpst,'BinWidth',1))./length(unique(tmpclu));

                    end
                end

                if plotUnitActivity
                    figure('name',[MiceOpt{midx} ' Activity Per Shank'])
                    imagesc(ShankOpt,depthavg,spShDepth,[quantile(spShDepth(:),0.05) quantile(spShDepth(:),0.95)])
                    set(gca,'YDir','normal')
                    colormap hot
                    h=colorbar;
                    h.Label.String = 'Spikes/Sec';
                    xlabel('Shank')
                    ylabel('Depth (micron)')
                end
                %% Get multiunit correlation - Copied from Petersen github
                goodonly=0
                if goodonly
                    spikeID = ismember(spikeCluster,Good_ID);
                else
                    spikeID = true(length(spikeCluster),1);
                    try
                        % Some form of quality control
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes<=quantile(clusinfo.n_spikes,0.01))))=0; %too little spikes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes>=quantile(clusinfo.n_spikes,0.99))))=0; %too many spikes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.amp>=quantile(clusinfo.amp,0.99))))=0; %ridiculous amplitudes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.fr>=quantile(clusinfo.fr,0.99))))=0; %ridiculous firing rates
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.ContamPct>=quantile(clusinfo.ContamPct,0.99))))=0; %ridiculous Contamination percentage
                        spikeID(ismember(spikeCluster,cluster_id(ismember(cellstr(Label),'noise'))))=0; %noise should not count, only MUA and good unit
                    catch
                        disp('This is non curated data, using all units from kilosort output')
                        spikeID = true(length(spikeCluster),1);
                    end
                end
                n_corr_groups = 80;
                startpoint =  min(depth);
                endpoint = max(depth);
                depth_group_edges = linspace(startpoint,endpoint,n_corr_groups+1);
                depth_group = discretize(spikeDepths,depth_group_edges);
                depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
                unique_depths = 1:length(depth_group_edges)-1;

                spike_binning = 1; % seconds
                corr_edges = nanmin(spikeTimes(spikeID)):spike_binning:nanmax(spikeTimes(spikeID));
                corr_centers = corr_edges(1:end-1) + diff(corr_edges);

                nshanks = length(ShankOpt);
                mua_corr = cell(1,nshanks);
                for shid = 1:nshanks
                    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
                    parfor curr_depth = 1:length(unique_depths)
                        binned_spikes_depth(curr_depth,:) = histcounts(spikeTimes(spikeID& depth_group == unique_depths(curr_depth) & spikeShank==shid), corr_edges);
                    end
                    %     % Z-score
                    %     binned_spikes_depth = (binned_spikes_depth - nanmean(binned_spikes_depth(:)))./nanstd(binned_spikes_depth(:));

                    binned_spikes_depth(:,nansum(binned_spikes_depth,1)>quantile(nansum(binned_spikes_depth,1),0.95))=0;
                    mua_corr{shid} = smooth2a(corrcoef(binned_spikes_depth'),3);
                end
                mua_corr = cat(3,mua_corr{:});
                mua_corr = reshape(mua_corr,size(mua_corr,1),[]);
                limup = [quantile(mua_corr(:),0.1) quantile(mua_corr(:),0.95)];
                % Plot multiunit correlation
                figure
                multiunit_ax = subplot(3,9,[1:5,10:14,19:23]);
                h=imagesc(1:length(depth_group_centers)*nshanks,depth_group_centers,mua_corr,limup);
                caxis([0,max(mua_corr(mua_corr~=1))]); colormap(hot);
                set(h,'Alphadata',~isnan(mua_corr))
                set(gca,'Color',[0.5 0.5 0.5])
                hold on
                for shid = 1:nshanks
                    line([size(mua_corr,1)*shid size(mua_corr,1)*shid],[startpoint endpoint],'color',[0.4 0.6 0],'LineStyle','--','LineWidth',3)
                end
                ylim([startpoint,endpoint]);
                xlim([1,length(depth_group_centers)*nshanks]);
                set(gca,'XTickLabel','')
                DChannels = endpoint-startpoint;


                set(multiunit_ax,'YDir','normal');
                title('MUA correlation');
                xlabel(multiunit_ax,'Multiunit depth X Probe');
                ylabel(multiunit_ax,'Multiunit depth');

                %% LFP
                if PlotLFP
                    myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},'*','ephys');
                    lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                    if isempty(lfpD)
                        disp('No LFP data found')
                    elseif length(lfpD)>length(subksdirs)
                        disp('Just take data from the last recording')
                        lfpD = lfpD(end);
                    elseif length(lfpD)<length(subksdirs)
                        disp('Should be a different amount of probes?')
                        keyboard
                    else
                        lfpD = lfpD(probeid);
                    end
                    freqBands = {[1.5 4], [4 12], [12 20], [20 30], [25 100],[150 200]};
                    FreqNames = {'Delta','Theta','Alpha','Beta','Gamma','Ripples'};
                    [lfpByChannel, allPowerEst, F, allPowerVar] = lfpBandPowerNP2(fullfile(lfpD.folder,lfpD.name),freqBands);
                    LFP_On=1;
                    %normalize LFP per frequency
                    lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,1))./nanstd(lfpByChannel,[],1);
                    if size(channelpos,1)>size(lfpByChannel,1)
                        channelpos = channelpostmp; %Take channel opt from last session
                    end

                    %Channel to depth/xpos:
                    xposopt = unique(floor(channelpos(:,1)./100).*100);
                    yposopt = unique(round(channelpos(:,2)./100).*100);

                    PwPerSh=arrayfun(@(Y) arrayfun(@(X) squeeze(nanmean(lfpByChannel((round(channelpos(:,1)./100).*100==X&round(channelpos(:,2)./100).*100==Y),:),1)),xposopt,'UniformOutput',0),yposopt,'UniformOutput',0);
                    PwPerSh = cat(2,PwPerSh{:});
                    PwPerSh = cat(2,PwPerSh{:});
                    PwPerSh = reshape(PwPerSh,length(freqBands),length(xposopt),length(yposopt));


                    figure('name',['LFP ' MiceOpt{midx} ' ' thisdate])
                    for shid = 1:length(ShankOpt)
                        subplot(1,length(ShankOpt),shid)
                        imagesc([],yposopt,squeeze(PwPerSh(:,shid,:))');
                        set(gca,'ydir','normal')
                        ylabel('Micron from tip')
                        xlabel('frequency')
                        colormap hot
                        title(['Shank ' num2str(ShankOpt(shid))])
                        set(gca,'XTick',1:length(FreqNames),'XTickLabel',FreqNames','XTickLabelRotation',45)
                        %                     title(sprintf('%d to %d Hz',round(freqBand{shid}(1)),round(freqBand{shid}(2))))

                    end
                end
            end
        end
    else
        % For every date a different dataset
        Dates4Mouse = DateOpt{midx};
        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'

            thisdate = Dates4Mouse{didx};
            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(LocalDir,MiceOpt{midx},thisdate);
            subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
            if length(subksdirs)<1
                clear subksdirs
                subksdirs.folder = myKsDir; %Should be a struct array
                subksdirs.name = 'Probe0';
            end
            DataThisProbe = cell(1,length(subksdirs));
            for probeid = 1:length(subksdirs)
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name)
                if ~isdir(myKsDir)
                    continue
                end

                %Saving directory
                thisprobe = subksdirs(probeid).name



                %% Get cluster information
                clear params
                params.loadPCs=true;
                params.thisdate = thisdate;
                PrepareClusInfo

                GoodUnits{midx}{didx}{probeid} = clusinfo;
                if ~any(Good_ID)
                    disp('No Good units found.. continue')
                    continue
                end

                % Show probe recording per location
                figure('name',[MiceOpt{midx} ' ' thisdate ' ' thisprobe ' Recording Per Shank'])
                scatter(Shank,depth,8,recses);

                myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                lfpD = dir(fullfile(myLFDir,'*','*','*.lf.*bin')); % ap file from spikeGLX specifically
                if isempty(lfpD)
                    lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
                end

                if isempty(lfpD)
                    disp('No LFP data found')
                elseif length(lfpD)>length(subksdirs)
                    disp('Just take data from the last recording')
                    lfpD = lfpD(end);
                elseif length(lfpD)<length(subksdirs)
                    disp('Should be a different amount of probes?')
                    keyboard
                else
                    lfpD = lfpD(probeid);
                end

                %% Show Activity on probe
                depthbins = min(clusinfo.depth):100:max(clusinfo.depth);
                depthavg = min(clusinfo.depth)+50:100:max(clusinfo.depth)-50;
                spShDepth  = nan(length(depthbins)-1,length(ShankOpt));
                for shid = 1:length(ShankOpt)

                    parfor deptid=1:length(depthbins)-1
                        tmpclu = spikeCluster(ismember(spikeShank,ShankOpt(shid))& spikeDepths>=depthbins(deptid)&spikeDepths<=depthbins(deptid+1));
                        tmpst = spikeTimes(ismember(spikeShank,ShankOpt(shid))& spikeDepths>=depthbins(deptid)&spikeDepths<=depthbins(deptid+1));

                        spShDepth(deptid,shid) = nanmean(histcounts(tmpst,'BinWidth',1))./length(unique(tmpclu));

                    end
                end

                if plotUnitActivity

                    figure('name',[MiceOpt{midx} ' ' thisdate ' ' thisprobe ' Activity Per Shank'])
                    imagesc(ShankOpt,depthavg,spShDepth,[quantile(spShDepth(:),0.05) quantile(spShDepth(:),0.95)])
                    set(gca,'YDir','normal')
                    colormap hot
                    h=colorbar;
                    h.Label.String = 'Spikes/Sec';
                    xlabel('Shank')
                    ylabel('Depth (micron)')
                end

                %% Get multiunit correlation - Copied from Petersen github
                goodonly=0
                if goodonly
                    spikeID = ismember(spikeCluster,Good_ID);
                else
                    spikeID = true(length(spikeCluster),1);
                    try
                        % Some form of quality control
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes<=quantile(clusinfo.n_spikes,0.01))))=0; %too little spikes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes>=quantile(clusinfo.n_spikes,0.99))))=0; %too many spikes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.amp>=quantile(clusinfo.amp,0.99))))=0; %ridiculous amplitudes
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.fr>=quantile(clusinfo.fr,0.99))))=0; %ridiculous firing rates
                        spikeID(ismember(spikeCluster,cluster_id(clusinfo.ContamPct>=quantile(clusinfo.ContamPct,0.99))))=0; %ridiculous Contamination percentage
                        spikeID(ismember(spikeCluster,cluster_id(ismember(cellstr(Label),'noise'))))=0; %noise should not count, only MUA and good unit
                    catch
                        disp('This is non curated data, using all units from kilosort output')
                        spikeID = true(length(spikeCluster),1);
                    end
                end
                n_corr_groups = 80;
                startpoint =  min(depth);
                endpoint = max(depth);
                depth_group_edges = linspace(startpoint,endpoint,n_corr_groups+1);
                depth_group = discretize(spikeDepths,depth_group_edges);
                depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
                unique_depths = 1:length(depth_group_edges)-1;

                spike_binning = 1; % seconds
                corr_edges = nanmin(spikeTimes(spikeID)):spike_binning:nanmax(spikeTimes(spikeID));
                corr_centers = corr_edges(1:end-1) + diff(corr_edges);

                nshanks = length(ShankOpt);
                mua_corr = cell(1,nshanks);
                for shid = 1:nshanks
                    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
                    parfor curr_depth = 1:length(unique_depths)
                        binned_spikes_depth(curr_depth,:) = histcounts(spikeTimes(spikeID& depth_group == unique_depths(curr_depth) & spikeShank==ShankOpt(shid)), corr_edges);
                    end
                    %     % Z-score
                    %     binned_spikes_depth = (binned_spikes_depth - nanmean(binned_spikes_depth(:)))./nanstd(binned_spikes_depth(:));

                    binned_spikes_depth(:,nansum(binned_spikes_depth,1)>quantile(nansum(binned_spikes_depth,1),0.95))=0;
                    mua_corr{shid} = smooth2a(corrcoef(binned_spikes_depth'),3);
                end
                mua_corr = cat(3,mua_corr{:});
                mua_corr = reshape(mua_corr,size(mua_corr,1),[]);
                limup = [quantile(mua_corr(:),0.1) quantile(mua_corr(:),0.95)];
                % Plot multiunit correlation
                figure('name',['MultiUnit ' MiceOpt{midx} ' ' thisdate ' ' thisprobe ])
                multiunit_ax = subplot(3,6,[1:5,7:11,13:17]);
                h=imagesc(1:length(depth_group_centers)*nshanks,depth_group_centers,mua_corr,limup);
                caxis([0,max(mua_corr(mua_corr~=1))]); colormap(hot);
                set(h,'Alphadata',~isnan(mua_corr))
                set(gca,'Color',[0.5 0.5 0.5])
                hold on
                for shid = 1:nshanks
                    line([size(mua_corr,1)*shid size(mua_corr,1)*shid],[startpoint endpoint],'color',[0.4 0.6 0],'LineStyle','--','LineWidth',3)
                end
                ylim([startpoint,endpoint]);
                xlim([1,length(depth_group_centers)*nshanks]);
                set(gca,'XTickLabel','')
                DChannels = endpoint-startpoint;


                set(multiunit_ax,'YDir','normal');
                title('MUA correlation');
                xlabel(multiunit_ax,'Multiunit depth X Probe');
                ylabel(multiunit_ax,'Multiunit depth');

                %% LFP
                if PlotLFP
                    freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};
                    FreqNames = {'Delta','Theta','Alpha','Beta','Gamma'};
                    LFPDir = fullfile(lfpD.folder,lfpD.name);
                    if ~isempty(strfind(LFPDir,'.lf')) % Saved out separately for NP1
                        % Get information from meta file
                        lfpD = dir(LFPDir);
                        [Imecmeta] = ReadMeta2(lfpD.folder,'lf');
                        lfpFs = str2num(Imecmeta.imSampRate);
                        nChansInFile = strsplit(Imecmeta.acqApLfSy,',');  % neuropixels phase3a, from spikeGLX
                        nChansInFile = str2num(nChansInFile{1})+1; %add one for sync

                        try
                            [lfpByChannel, allPowerEst, F, allPowerVar] = ...
                                lfpBandPower(fullfile(lfpD.folder,lfpD.name), lfpFs, nChansInFile, freqBands);
                            allPowerEst = allPowerEst(:,1:nChansInFile)'; % now nChans x nFreq
                            LFP_On =1;

                        catch ME
                            disp(ME)
                            LFP_On =0;
                        end

                    elseif ~isempty(strfind(LFPDir,'.ap')) % Saved out in .ap for NP1
                        [lfpByChannel, allPowerEst, F, allPowerVar] = lfpBandPowerNP2(LFPDir,freqBands);
                        LFP_On=1;

                    else
                        disp('No file found...')
                        LFP_On=0;
                    end

                    if LFP_On & PlotLFP
                        %normalize LFP per frequency
                        lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,1))./nanstd(lfpByChannel,[],1);
                        if size(channelpos,1)>size(lfpByChannel,1)
                            channelpos = channelpostmp; %Take channel opt from last session
                        end

                        %Channel to depth/xpos:
                        xposopt = unique(floor(channelpos(:,1)./100).*100);
                        yposopt = unique(round(channelpos(:,2)./100).*100);
                        Shank = floor(channelpos(:,1)./250);
                        ShankOpt = unique(Shank);
                        PwPerSh = arrayfun(@(Y) arrayfun(@(X) (nanmean(lfpByChannel(Shank==X & round(channelpos(:,2)./100).*100==Y,:),1)),ShankOpt,'UniformOutput',0),yposopt,'UniformOutput',0);
                        PwPerSh = cat(2,PwPerSh{:});
                        PwPerSh = cat(2,PwPerSh{:});
                        PwPerSh = reshape(PwPerSh,length(freqBands),length(ShankOpt),length(yposopt));


                        figure('name',['LFP ' MiceOpt{midx} ' ' thisdate ' ' thisprobe])
                        for shid = 1:length(ShankOpt)
                            subplot(1,length(ShankOpt),shid)
                            imagesc([],yposopt,squeeze(PwPerSh(:,shid,:))');
                            set(gca,'ydir','normal')
                            ylabel('Micron from tip')
                            xlabel('frequency')
                            colormap hot
                            title(['Shank ' num2str(ShankOpt(shid))])
                            set(gca,'XTick',1:length(FreqNames),'XTickLabel',FreqNames','XTickLabelRotation',45)
                            %                     title(sprintf('%d to %d Hz',round(freqBand{shid}(1)),round(freqBand{shid}(2))))
                        end
                    end

                    drawnow
                end

                %% Rastermap
                if 0
                    BinSizeOpt = [0.1, 0.3, 0.5, 1];
                    iPCOpt = [1:10];
                    sortedThis = nan(length(BinSizeOpt),length(iPCOpt),2,length(Good_IDx));
                    SortingSimilarityScore = nan(length(BinSizeOpt),length(iPCOpt));
                    GraphSimilarity = nan(length(BinSizeOpt),length(iPCOpt));
                    for bid = 1:length(BinSizeOpt)
                        for pid = 1:length(iPCOpt)

                            for cv = 1:2
                                % S = units X Time
                                %Spike rate / histogram
                                %                 binSize = 0.1;
                                binSize = BinSizeOpt(bid);
                                Edges = min(spikeTimes(:)):binSize:max(spikeTimes(:));

                                SpikeRatePerTP = arrayfun(@(Y) histcounts(spikeTimes(spikeCluster'== Y),...
                                    Edges)./binSize,cluster_id(Good_IDx),'UniformOutput',0);
                                SpikeRatePerTP = cat(1,SpikeRatePerTP{:});

                                DataThisProbe{probeid} = SpikeRatePerTP;

                                % Run rastermap

                                %                 [isort1, isort2, Sm] = mapTmap(csv);
                                %                 figure; imagesc(-csv(isort1,:),[-1 0])
                                %                 colormap gray
                                %                 ops.iPC = 50;
                                ops.iPC = 1:iPCOpt(pid);
                                if   any(ops.iPC>size(SpikeRatePerTP,1))
                                    continue
                                end
                                [sortedThis(bid,pid,cv,:), isort2, Sm{cv}] = mapTmap(SpikeRatePerTP(:,cv:2:end),ops);
                                %                             limup = quantile(SpikeRatePerTP(:),0.99);
                                %                             subplot(1,2,cv)
                                %                             imagesc(-SpikeRatePerTP(isort1,:),[-limup 0])
                                %                             colormap gray
                                %                             title(['Bin=' num2str(binSize) ', iPC=' num2str(ops.iPC) ', cv' num2str(cv)])
                                %                             makepretty
                                %                             drawnow
                                %                 subplot(1,3,2)
                                %                 imagesc(-SpikeRatePerTP(:,isort2),[-limup 0])
                                %                 colormap gray
                                %                 title('Sorted by time Rastermap')
                                %
                                %                 subplot(1,3,3)
                                %                 imagesc(-SpikeRatePerTP(isort1,isort2),[-limup 0])
                                %
                                %                 %                 imagesc(-Sm,[-quantile(Sm(:),0.99) 0])
                                %                 colormap gray
                                %                 title('Sorted by both Rastermap')
                                %
                                %
                                %                 zSpikeRate = (SpikeRatePerTP-nanmean(SpikeRatePerTP,2))./nanstd(SpikeRatePerTP,[],2);
                                %                 limup = quantile(abs(zSpikeRate(:)),0.95);
                                %                 figure; imagesc(zSpikeRate(SpikeRatePerTP,:),[-limup limup])
                                %                 colormap redblue

                                %                             drawnow
                                %
                            end
                            SortingSimilarityScore(bid,pid)=levenshtein_distance(squeeze(sortedThis(bid,pid,1,:))',squeeze(sortedThis(bid,pid,2,:))');
                            if size(Sm{1},2)>size(Sm{2},2) %More time points in SM1
                                Sm{1} = Sm{1}(:,1:size(Sm{2},2));
                            end
                            GraphSimilarity(bid,pid) = corr2(Sm{1},Sm{2});

                        end

                    end

                    % Best one:
                    [bid,pid] = find(GraphSimilarity==max(GraphSimilarity(:)));
                    [bid2,pid2] = find(SortingSimilarityScore==min(SortingSimilarityScore(:)));
                    bid = unique(bid);
                    pid = unique(pid);
                    bid2=unique(bid2);
                    pid2=unique(pid2);


                    if all(bid==bid2) && all(pid==pid2)
                        disp('Great, Levenshtein and correlation agree!')
                    elseif length(bid)>1
                        bid = find(ismember(bid2,bid));
                        disp('No agreement on bid, take best values based on both')
                    elseif length(pid)>1
                        pid = find(ismember(pid2,pid));
                        disp('No agreement on pid, take best values based on both')

                    end
                    binSize = BinSizeOpt(bid);
                    Edges = min(spikeTimes(:)):binSize:max(spikeTimes(:));

                    SpikeRatePerTP = arrayfun(@(Y) histcounts(spikeTimes(spikeCluster'== Y),...
                        Edges)./binSize,cluster_id(Good_IDx),'UniformOutput',0);
                    SpikeRatePerTP = cat(1,SpikeRatePerTP{:});

                    % Run rastermap

                    %                 [isort1, isort2, Sm] = mapTmap(csv);
                    %                 figure; imagesc(-csv(isort1,:),[-1 0])
                    %                 colormap gray
                    %                 ops.iPC = 50;
                    ops.iPC = 1:iPCOpt(pid);
                    [isort1, isort2, Sm_Best] = mapTmap(SpikeRatePerTP,ops);
                    SimilarityToStandard=levenshtein_distance(isort1,1:size(SpikeRatePerTP,1));
                    disp(['Similarity to standars order:' num2str(SimilarityToStandard)])


                    sortedThis(bid,pid,cv,:)=isort1;

                    figure('name', ['Bin=' num2str(BinSizeOpt(bid)) ', iPC=' num2str(iPCOpt(pid))])
                    limup = quantile(SpikeRatePerTP(:),0.99);
                    subplot(2,2,[2,4])
                    imagesc(-SpikeRatePerTP,[-limup 0])

                    %                 imagesc(-SpikeRatePerTP(isort1,:),[-limup 0])
                    colormap gray
                    set(gca,'XTickLabel',cellfun(@(X) num2str(str2num(X).*binSize),get(gca,'XTickLabel'),'UniformOutput',0))
                    xlabel('Time (s)')
                    ylabel('Sorted units')
                    title(['Bin=' num2str(binSize) ', iPC=' num2str(ops.iPC)])
                    makepretty
                    drawnow


                    subplot(2,2,[1])
                    h=imagesc(SortingSimilarityScore)
                    set(gca,'XTick',1:length(iPCOpt),'XTickLabel',arrayfun(@(X) num2str(X),iPCOpt,'UniformOutput',0))
                    set(gca,'YTick',1:length(BinSizeOpt),'YTickLabel',arrayfun(@(Y) num2str(Y),BinSizeOpt,'UniformOutput',0))
                    colorbar
                    colormap gray
                    xlabel('PCi')
                    ylabel('Binsize (s)')
                    title('Levenshtein Distance')
                    makepretty

                    subplot(2,2,[3])
                    h=imagesc(GraphSimilarity)
                    set(gca,'XTick',1:length(iPCOpt),'XTickLabel',arrayfun(@(X) num2str(X),iPCOpt,'UniformOutput',0))
                    set(gca,'YTick',1:length(BinSizeOpt),'YTickLabel',arrayfun(@(Y) num2str(Y),BinSizeOpt,'UniformOutput',0))
                    colorbar
                    colormap gray
                    xlabel('PCi')
                    ylabel('Binsize (s)')
                    title('Graph Similarity')
                    makepretty
                    drawnow

                    saveas(gcf,fullfile(myKsDir,'RasterMapSorting.fig'))
                    saveas(gcf,fullfile(myKsDir,'RasterMapSorting.bmp'))
                end
            end
        end
        %% number of good units across days
        figure('name','Units across days')
        nrdays = length(GoodUnits{midx});
        countday = 1;
        for didx = 1:nrdays
            tmpCell = GoodUnits{midx}{didx};
            if isempty(tmpCell)
                continue
            end
            nprobe = length(tmpCell);
            for pidx = 1:nprobe
                tmpUnits = tmpCell{pidx};
                depthGU = tmpUnits.depth(logical(tmpUnits.Good_ID));
                ShankGU = tmpUnits.Shank(logical(tmpUnits.Good_ID));
                if isempty(depthGU)
                    continue
                end

                scatter(countday-0.15+0.3*(pidx-1)+0.05*ShankGU,depthGU,8,'filled')
                hold on

                text(countday-0.15+0.3*(pidx-1),max(depthGU)+250,['n=' num2str(sum(tmpUnits.Good_ID))])

            end
            countday = countday+1;
        end
        ylim([0 max(get(gca,'ylim'))+300])
        xlim([0 countday])
        xlabel('Day/Shank')
        ylabel('Depth from probetip')
        title(['Good Units Across Days/Shank ' MiceOpt{midx}])
        makepretty

        saveas(gcf,fullfile(LocalDir,MiceOpt{midx},'UnitsAcrossDays.fig'))
        saveas(gcf,fullfile(LocalDir,MiceOpt{midx},'UnitsAcrossDays.bmp'))

    end
end

