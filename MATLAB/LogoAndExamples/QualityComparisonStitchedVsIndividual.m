%% User Input
%% Path information
DataDir = {'H:\MatchingUnits\RawData'};% Check DataDir2Use
SaveDir = 'H:\MatchingUnits\Output\'
tmpdatafolder = 'H:\MatchingUnits\Tmp'; % temporary folder
KilosortDir = 'H:\MatchingUnits\KilosortOutput';% 'E:\Data\KiloSortOutput';%
RunPyKSChronicStitchedOpt = [0,1];
%% Information on experiments
MiceOpt = {'AL032','AV008','CB016','EB019','JF067'};%{'AV009','AV015','EB014','EB019','AL032','AV008','JF067','CB016'}%,'AV008','JF067','CB016''EB019'}; %CB016 %AL032 'AV008' JF067 Add all mice you want to analyse
% nidq_sync_used = zeros(1,length(MiceOpt)); % Was an external nidq used for syncing (typically sync feeds directly into IMEC)
% nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1; % Except for these mice...
DataDir2Use = repmat(1,[1,length(MiceOpt)]); % In case you have multiple DataDir, index which directory is used for each mouse
RecordingType = repmat({'Acute'},1,length(MiceOpt)); % And whether recordings were acute (default)
RecordingType(ismember(MiceOpt,{'AL032','EB019','CB016','AV008','JF067'}))={'Chronic'}; %EB014', % Or maybe Chronic?
PrepareClusInfoparams.ReLoadAlways = 1; % If 1, SP & Clusinfo are always loaded from KS output

%% DateOpt
clear DateOpt
for idx = 1:length(DataDir)
    DateOpt{idx} = cellfun(@(X) dir(fullfile(DataDir{idx},X,'*-*')),MiceOpt(DataDir2Use==idx),'UniformOutput',0);
end
DateOpt = cat(2,DateOpt{:});
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);

%% Compare KS Stitched versus non-stitched
for midx = 1:length(MiceOpt)
    QMetricsAll = cell(1,2);
    MSEAll = cell(1,2);
    DepthVar = cell(1,2);
    SpikeRate = cell(1,2);
    for subid = 1:length(RunPyKSChronicStitchedOpt)
        %% Loading data from kilosort/phy easily
        myKsDir = fullfile(KilosortDir,MiceOpt{midx});
        subksdirs = dir(fullfile(myKsDir,'**','Probe*')); %This changed because now I suddenly had 2 probes per recording
        if length(subksdirs)<1
            clear subksdirs
            subksdirs.folder = myKsDir; %Should be a struct array
            subksdirs.name = 'Probe0';
        end

        ProbeOpt = (unique({subksdirs(:).name}));

        myKsDir = fullfile(KilosortDir,MiceOpt{midx});
        % Check for multiple subfolders?
        subsesopt = dir(fullfile(myKsDir,'**','channel_positions.npy'));
        subsesopt=arrayfun(@(X) subsesopt(X).folder,1:length(subsesopt),'Uni',0);
        if strcmp(RecordingType{midx},'Chronic')
            if ~RunPyKSChronicStitchedOpt(subid) %MatchUnitsAcrossDays
                disp('Unit matching in Matlab')
                subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),subsesopt,'UniformOutput',0))) = []; %Use separate days and match units via matlab script
            else
                disp('Using chronic pyks option')
                subsesopt = subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),subsesopt,'UniformOutput',0))); %Use chronic output from pyks
            end
        end

        % Copy file and then run pykilosort on it
        channelposition = cell(1,length(subsesopt));
        for id = 1:length(subsesopt)
            tmpfile = dir(fullfile(subsesopt{id},'channel_positions.npy'));
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

        %% Create saving directoryed
        clear params
        thisIMRO = '';
        thisdate = [];
        if ~exist(fullfile(SaveDir,MiceOpt{midx}))
            mkdir(fullfile(SaveDir,MiceOpt{midx}))
        end
        PrepareClusInfoparams.SaveDir = fullfile(SaveDir,MiceOpt{midx});

        %% Load prepared cluster information
        tmpprep = dir(fullfile(myKsDir,'**','PreparedData.mat'));
        if RunPyKSChronicStitchedOpt(subid)
            FolderIdx = find(cell2mat(cellfun(@(X) ~isempty(strfind(X,'Chronic')),{tmpprep(:).folder},'Uni',0)));
        else
            FolderIdx = find(~cell2mat(cellfun(@(X) ~isempty(strfind(X,'Chronic')),{tmpprep(:).folder},'Uni',0)));
        end


        % all units
        for fid = 1:length(FolderIdx)
            % Load Clusinfo
            tmp = load(fullfile(tmpprep(FolderIdx(fid)).folder,tmpprep(FolderIdx(fid)).name));
            OriIDAll = tmp.clusinfo.cluster_id;
            % Load quality metrics file
            d = dir(fullfile(tmpprep(FolderIdx(fid)).folder,'**','templates._bc_qMetrics.parquet'));
            for did = 1:length(d)
                qMetricsPath = d(did).folder;
                [~, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(qMetricsPath);

                qMclusterID = qMetric.clusterID-1; %0-index
                ThisIdx = ismember(qMclusterID,OriIDAll); % Find the same cluster IDs

                if fid==1 && did == 1
                    qMetricAllUnits = qMetric(ThisIdx,:); % Add sessions together, only take good units
                else
                    qMetricAllUnits = cat(1,qMetricAllUnits,qMetric(ThisIdx,:));
                end

            end

            % Amplitude of spikes
            tmp = load(fullfile(tmpprep(FolderIdx(fid)).folder,tmpprep(FolderIdx(fid)).name));
            clusinfo = tmp.clusinfo;
            sp = tmp.sp;

            sp.st = sp.st+(sp.SessionID-1)*max(sp.st);
            MSEtmp = nan(1,length(clusinfo.cluster_id));
            DepthVartmp = nan(1,length(clusinfo.cluster_id));
            SpikeRatetmp = nan(1,length(clusinfo.cluster_id));
            for clusid = 1:length(clusinfo.cluster_id)
                idx = sp.spikeTemplates == clusinfo.cluster_id(clusid);
                if sum(idx)<10
                    continue
                end
                Amplitudes = sp.spikeAmps(idx);
                Edges = [min(Amplitudes):2:max(Amplitudes)];
                if length(Edges)==1
                    keyboard
                end
                Vec = [min(Amplitudes)+1:2:max(Amplitudes)-1];
                h = histcounts(Amplitudes,Edges);
                h = (h-nanmin(h))./(nanmax(h)-nanmin(h));
                %                     figure; plot(Vec,h); hold on

                [m,s] = normfit(Amplitudes);
                y = normpdf(Vec,m,s);
                % Scale
                y = (y-nanmin(y(:)))./(nanmax(y(:))-nanmin(y(:)));
                %                     plot(Vec,y,'.')

                % Fit?
                MSEtmp(clusid) = nansum((y-h).^2);

                % Spike depths?
                DepthVartmp(clusid) = nanvar(sp.spikeDepths(idx));

                % Firing rate?
                SpikeRatetmp(clusid) = nanmean(histcounts(sp.st(idx),[0:1:max(sp.st)]));
            end
            if fid==1 
                MSEAllUnits = MSEtmp; % Add sessions together, only take good units
                DepthVarAllUnits = DepthVartmp;
                SpikeRateAllUnits = SpikeRatetmp;
            else
                MSEAllUnits = cat(2,MSEAllUnits,MSEtmp);
                DepthVarAllUnits = cat(2,DepthVarAllUnits,DepthVartmp);
                SpikeRateAllUnits = cat(2,SpikeRateAllUnits,SpikeRatetmp);
            end
        end


        QMetricsAll{subid} = qMetricAllUnits;
        MSEAll{subid} = MSEAllUnits;
        DepthVar{subid} = DepthVarAllUnits;
        SpikeRate{subid} = SpikeRateAllUnits;

    end
    cols = distinguishable_colors(2);

    %% Amplitude
    if midx==1
    ampfig = figure('name','Functional Scores');
    else
        figure(ampfig)
    end  
    subplot(length(MiceOpt),3,(midx-1)*3+1)
    clear tmp
    clear h
    lims = [floor(quantile(cat(2,MSEAll{:}),0.005)) ceil(quantile(cat(2,MSEAll{:}),0.995))];
    Edges = [lims(1):2:lims(2)];
    for hid = 1:2
        %             if numel(tmp{hid})<nmax % Fill up with nans
        %                 tmp{hid}(end+1:nmax)=nan;
        %             end

        h(hid) = histogram(MSEAll{hid},Edges,'FaceAlpha',0.5,'FaceColor',cols(hid,:));
        hold on
        line([nanmedian(MSEAll{hid}), nanmedian(MSEAll{hid})],get(gca,'ylim'),'color',cols(hid,:));
    end
    p = ranksum(MSEAll{1},MSEAll{2});
    if p<0.05
        if nanmedian(MSEAll{1})>nanmedian(MSEAll{2})
            disp(['Amplitude Distribution Fits significantly worse for non-stiched'])
        else
            disp(['Amplitude Distribution Fits significantly worse for stiched'])
        end
    end
    title([MiceOpt{midx} ', p=' num2str(round(p*1000)./1000)]) 
    xlim(lims)
    xlabel('MSE amplitude distribution fit')
    makepretty

    %% depth var
    subplot(length(MiceOpt),3,(midx-1)*3+2)
    clear tmp
    clear h
    lims = [floor(quantile(cat(2,DepthVar{:}),0.005)) ceil(quantile(cat(2,DepthVar{:}),0.995))];
    Edges = [lims(1):0.05:lims(2)];
    for hid = 1:2
        %             if numel(tmp{hid})<nmax % Fill up with nans
        %                 tmp{hid}(end+1:nmax)=nan;
        %             end

        h(hid) = histogram(DepthVar{hid},Edges,'FaceAlpha',0.5,'FaceColor',cols(hid,:));
        hold on
        line([nanmedian(DepthVar{hid}), nanmedian(DepthVar{hid})],get(gca,'ylim'),'color',cols(hid,:));
    end
    p = ranksum(DepthVar{1},DepthVar{2});
    if p<0.05
        if nanmedian(DepthVar{1})>nanmedian(DepthVar{2})
            disp(['Depth Variance significantly larger for non-stiched'])
        else
            disp(['Depth Variance significantly larger for stiched'])
        end
    end
    title([MiceOpt{midx} ', p=' num2str(round(p*1000)./1000)]) 
    xlim(lims)
    xlabel('Spike depth variance')
    makepretty

   %% SpikeRate
    subplot(length(MiceOpt),3,(midx-1)*3+3)
    clear tmp
    clear h
    lims = [floor(quantile(cat(2,SpikeRate{:}),0.005)) ceil(quantile(cat(2,SpikeRate{:}),0.995))];
    Edges = [lims(1):0.5:lims(2)];
    for hid = 1:2
        %             if numel(tmp{hid})<nmax % Fill up with nans
        %                 tmp{hid}(end+1:nmax)=nan;
        %             end

        h(hid) = histogram(SpikeRate{hid},Edges,'FaceAlpha',0.5,'FaceColor',cols(hid,:));
        hold on
        line([nanmedian(SpikeRate{hid}), nanmedian(SpikeRate{hid})],get(gca,'ylim'),'color',cols(hid,:));
    end
    p = ranksum(SpikeRate{1},SpikeRate{2});
    if p<0.05
        if nanmedian(SpikeRate{1})>nanmedian(SpikeRate{2})
            disp(['Spike Rate significantly larger for non-stiched'])
        else
            disp(['Spike Rate significantly larger for stiched'])
        end
    end
    title([MiceOpt{midx} ', p=' num2str(round(p*1000)./1000)]) 
    xlim(lims)
    xlabel('Spike Rate')
    makepretty


    if midx==length(MiceOpt)
        legend([h(:)],{'Non-stitched','Stitched'})
    end

    %% Quality metric differences?
    figure('name',[MiceOpt{midx} ' Qmetric distributions'])
    disp(MiceOpt{midx})
    nmax = max(cellfun(@height,QMetricsAll));
    for qmid = 1:size(QMetricsAll{1},2)
        subplot(ceil(sqrt(size(QMetricsAll{1},2))),round(sqrt(size(QMetricsAll{1},2))),qmid)

        clear tmp
        clear h
        for hid = 1:2
            tmp{hid} = table2array(QMetricsAll{hid}(:,qmid))';
%             if numel(tmp{hid})<nmax % Fill up with nans
%                 tmp{hid}(end+1:nmax)=nan;
%             end

            h(hid) = histogram(table2array(QMetricsAll{hid}(:,qmid)),'FaceAlpha',0.5,'FaceColor',cols(hid,:));
            hold on
            line([nanmedian(table2array(QMetricsAll{hid}(:,qmid))), nanmedian(table2array(QMetricsAll{hid}(:,qmid)))],get(gca,'ylim'),'color',cols(hid,:));
        end
        try
            %             alltmp = cat(1,tmp{:})';
            %             violinplot(alltmp);
            p = ranksum(tmp{1},tmp{2});
            if p<0.05
                if nanmedian(tmp{1})>nanmedian(tmp{2})
                    disp([QMetricsAll{1}.Properties.VariableNames{qmid} ' significant: larger for non-stiched'])
                else
                    disp([QMetricsAll{1}.Properties.VariableNames{qmid} ' significant: larger for stiched'])
                end
            end
            %             ylim([quantile(alltmp(:),0.001) quantile(alltmp(:),0.95)])

        catch ME
        end
        %         set(gca,'XTickLabel',{'Non-stitched','Stitched'})

        makepretty
        xlabel(QMetricsAll{1}.Properties.VariableNames{qmid})
        title(['p=' num2str(round(p*1000)./1000)])
    end
    legend([h(:)],{'Non-stitched','Stitched'})

  
end
