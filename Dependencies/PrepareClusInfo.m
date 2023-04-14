function [clusinfo,sp,AllChannelPos] = PrepareClusInfo(KiloSortPaths,Params,RawDataPaths)
% Prepares cluster information for subsequent analysis
%% Inputs:
% KiloSortPaths = List of directories pointing at kilosort output (same format as what you get when
% calling 'dir') for all sessions you want to analyze as one session (e.g.
% chronic recordings with same IMRO table: give a list with all sessions)

% Params: with:
% Params.loadPCs=1;
% Params.RunPyKSChronicStitched = 1
% Params.DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
% Params.RedoQM = 0; %if 1, redo quality metrics if it already exists
% Params.RunQualityMetrics = 1; % If 1, Run the quality metrics
% Params.InspectQualityMetrics =0; % Inspect the quality metrics/data set using the GUI
% Params.UnitMatch = 1; % Matching chronic recording using QM instead of using pyks chronic output
% Params.RedoUnitMatch = 0; % Redo unitmatch
% Params.SaveDir% Directory to save QM and UnitMatch results
% Params.tmpdatafolder %Directory to temporarily decompress data --> must
% be large enough!

% RawDataPaths = list of same directories pointing at the raw ephys data
% directory. If this input is missing, this code will try to find the raw
% ephys data in the params.py file, first line (dat_path).

%% Output
% Clusinfo: Struct with cluster information. Read in from (Py)KS but
% optionally Quality and Matches (oversplits/matches across recordings) are identified
% Sp: Struct with spike information for all recordings


%% Check inputs
if ~iscell(KiloSortPaths) %isstruct(KiloSortPaths) || isfield(KiloSortPaths(1),'name') || isfield(KiloSortPaths(1),'folder')
    error('This is not a cell... give correct input please')
end

try
    if nargin<2
        disp('No params given. Use default - although this is not advised...')
        Params.loadPCs=1;
        Params.RunPyKSChronicStitched = 1;
        Params.DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
        Params.RedoQM = 0; %if 1, redo quality metrics if it already exists
        Params.RunQualityMetrics = 1; % If 1, Run the quality metrics
        Params.InspectQualityMetrics =0; % Inspect the quality metrics/data set using the GUI
        Params.UnitMatch = 1; % Matching chronic recording using QM instead of using pyks chronic output
        Params.RedoUnitMatch = 0; % Redo unitmatch
        Params.SaveDir = KiloSortPaths(1); % Directory to save QM and UnitMatch results
        Params.tmpdatafolder = KiloSortPaths(1); %
        Params.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows
        Params.saveSp = 0;
    end
    if Params.RunQualityMetrics
        Params.loadPCs=1; %If you want to run QM you need this
    end

    if nargin<3
        disp('Finding raw ephys data using the params.py file from (py)kilosort output')
        UseParamsKS = 1;
    else
        UseParamsKS = 0;
    end
catch ME
    disp(ME)
    UseParamsKS = 1;
end

%% Initialize everything
channelmap=[];
channelpos = [];

AllKiloSortPaths = cell(1,0);
AllChannelPos = cell(1,0);
countid=1;
% figure;
cols = jet(length(KiloSortPaths));
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
        continue
    end
    %% initialize for new round
    Good_ID = [];
    Shank=[];
    recsesAll = [];
    channel = [];
    Label = [];
    depth = [];
    AllUniqueTemplates = [];
    cluster_id = [];
    recses = [];

    %% save data paths information
    if Params.RunPyKSChronicStitched %CALL THIS STITCHED --> Only works when using RunPyKS2_FromMatlab as well from this toolbox
        if UseParamsKS
            try
                spikeStruct = loadParamsPy(fullfile(KiloSortPaths{subsesid},'params.py'));
                rawD = spikeStruct.dat_path;
                rawD = strsplit(rawD,',');
                for rid=1:length(rawD)
                    rawD{rid} = rawD{rid}(strfind(rawD{rid},'"')+1:end);
                    rawD{rid} = rawD{rid}(1:strfind(rawD{rid},'"')-1);
                    rawD{rid} = dir(rawD{rid});
                    if isempty(rawD{rid})
                        rawD{rid} = dir(strrep(rawD{rid},'bin','cbin'));
                    end
                end
                rawD = cat(2,rawD{:});
            catch
                sessionsIncluded = dir(fullfile(KiloSortPaths{subsesid},'SessionsIncluded.mat'));
                sessionsIncluded = load(fullfile(sessionsIncluded.folder,sessionsIncluded.name));
                rawD = arrayfun(@(X) dir(sessionsIncluded.ThesePaths{X}),1:length(sessionsIncluded.ThesePaths),'UniformOutput',0);
                rawD = cat(2,rawD{:});
            end
        else
            rawD = dir(fullfile(RawDataPaths(subsesid).folder,RawDataPaths(subsesid).name));
        end

        RawDataPaths = rawD;
        AllKiloSortPaths = {AllKiloSortPaths{:} repmat(KiloSortPaths{subsesid},1,length(rawD))};
    else
        if UseParamsKS
            spikeStruct = loadParamsPy(fullfile(KiloSortPaths{subsesid},'params.py'));
            rawD = spikeStruct.dat_path;
            rawD = rawD(strfind(rawD,'"')+1:end);
            rawD = rawD(1:strfind(rawD,'"')-1);
            tmpdr = rawD;
            rawD = dir(rawD);
            if isempty(rawD)
                rawD = dir(strrep(tmpdr,'bin','cbin'));
            end


            % Save for later
            RawDataPaths(subsesid) = rawD;

            if isempty(rawD)
                disp('Bug...?')
                keyboard
            end
        else
            rawD = dir(fullfile(RawDataPaths(subsesid).folder,RawDataPaths(subsesid).name));
        end
        AllKiloSortPaths = {AllKiloSortPaths{:} KiloSortPaths{subsesid}};
    end
    DecompressionFlag = 0;

    %% Channel data
    myClusFile = dir(fullfile(KiloSortPaths{subsesid},'channel_map.npy'));
    channelmaptmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

    myClusFile = dir(fullfile(KiloSortPaths{subsesid},'channel_positions.npy'));
    channelpostmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    if length(channelmaptmp)<length(channelpostmp)
        channelmaptmp(end+1:length(channelpostmp))=length(channelmaptmp):length(channelpostmp)-1;
    end

    %% Is it correct channelpos though...? Check using raw data
    channelpostmpconv = ChannelIMROConversion(rawD(1).folder,1); % For conversion when not automatically done
    AllChannelPos{countid} = channelpostmpconv;

    %% Load existing?
    if exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat')) && ~Params.RedoQM && ~Params.ReLoadAlways
        % Check if parameters are the same, of not we have to redo it
        % anyway
        tmpparam = matfile(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'));
        tmpparam = tmpparam.Params;

        if tmpparam.RunQualityMetrics == Params.RunQualityMetrics && tmpparam.RunPyKSChronicStitched==Params.RunPyKSChronicStitched &&...
                tmpparam.separateIMRO == Params.separateIMRO
            disp(['Found existing data in ' KiloSortPaths{subsesid} ', Using this...'])
            %         load(fullfile(Params.SaveDir,'PreparedData.mat'))

            % Use UnitMatch Output if available
            %         if Params.UnitMatch
            %             disp('Using UnitMatch Clusters!')
            %             sp.clu = sp.UniqClu; %Temporary replace for rest of code
            %             clusinfo.cluster_id = clusinfo.UniqueID;
            %         end
            countid=countid+1;
            continue
        end
    end

    %% Load Spike Data
    sp = loadKSdir(fullfile(KiloSortPaths{subsesid}),Params); % Load Spikes with PCs
    [sp.spikeAmps, sp.spikeDepths, sp.templateDepths, sp.templateXpos, sp.tempAmps, sp.tempsUnW, sp.templateDuration, sp.waveforms] = templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.xcoords, sp.spikeTemplates, sp.tempScalingAmps); %from the spikes toolbox

    %% Remove noise; spikes across all channels'
    sp = RemoveNoiseAmplitudeBased(sp);


    %% Bombcell parameters
    % clear paramBC
    paramBC = bc_qualityParamValuesForUnitMatch(dir(strrep(fullfile(rawD(1).folder,rawD(1).name),'cbin','meta')),fullfile(Params.tmpdatafolder,strrep(rawD(1).name,'cbin','bin')));

    %% Load Cluster Info
    myClusFile = dir(fullfile(KiloSortPaths{subsesid},'cluster_info.tsv')); % If you did phy (manual curation) we will find this one... We can trust you, right?
    if isempty(myClusFile)
        disp('This data is not curated with phy! Hopefully you''re using automated quality metrics to find good units!')
        curratedflag=0;
        myClusFile = dir(fullfile(KiloSortPaths{subsesid},'cluster_group.tsv'));
        clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after phy)
        [clusinfo,sp] = ConvertTemplatesAfterPhy(clusinfo,sp);

        clusidtmp = clusinfo.cluster_id;
        cluster_id = cat(1,cluster_id,clusinfo.cluster_id);
        tmpLabel = char(length(clusinfo.cluster_id));
        KSLabelfile = tdfread(fullfile(KiloSortPaths{subsesid},'cluster_KSLabel.tsv'));
        tmpLabel(ismember(clusinfo.cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,clusinfo.cluster_id));
        Label = [Label,tmpLabel];
        totSpkNum = histc(sp.clu,sp.cids);
        Good_IDtmp = ismember(tmpLabel,'g') & totSpkNum'>paramBC.minNumSpikes; %Identify good clusters

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
            sp.spikeDepths(ismember(sp.clu,clusidtmp(clusid))) = depthtmp(clusid);

        end
        depth = cat(1,depth, depthtmp);
        channel = cat(1,channel,channeltmp);

    else
        disp('You did manual curation. You champion. If you have not enough time, maybe consider some automated algorithm...')
        CurationDone = 1;
        save(fullfile(KiloSortPaths{subsesid},'CuratedResults.mat'),'CurationDone')
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

        depth = [depth,depthtmp];

        channel = [channel,clusinfo.ch];
        Good_IDtmp = ismember(cellstr(clusinfo.group),'good');
        channeltmp = clusinfo.ch;
    end
    channelpostmp = channelpostmpconv;
    ypostmp = channelpostmp(:,2);
    xpostmp = channelpostmp(:,1);
    xposopt = (floor(xpostmp./250));% Assuming no new shank if not at least 100 micron apart
    Shanktmp = floor(xpostmp(channeltmp+1)./250);

    %     [~,minid] = arrayfun(@(X) (abs(floor(xpostmp(X)./250)-xposopt)),channeltmp+1,'UniformOutput',0);
    Shank = cat(1,Shank,Shanktmp);

    recses = cat(1, recses, repmat(countid,length(channeltmp),1));

    %     channelrecses = cat(2,channelrecses,repmat(countid,1,size(channelpostmp,1)));
    %     channelpos = cat(1,channelpos,channelpostmp);
    channelpos = channelpostmp; %only every need this here for mulitple sessions when chronic.
    channelmap = channelmaptmp;

    %     scatter(xtmp,depthtmp,10,cols(countid,:))
    %     hold on
    %     drawnow

    if Params.RunQualityMetrics
        theseuniqueTemplates = [];
        unitTypeAcrossRec= [];

        %% Quality metrics - Bombcell (https://github.com/Julie-Fabre/bombcell)
        for id = 1:length(rawD)
            ephysap_tmp = [];
            ephysap_path = fullfile(rawD(id).folder,rawD(id).name);
            if length(rawD)>1
                savePath = fullfile(myClusFile(1).folder,num2str(id));
            else
                savePath =myClusFile(1).folder;
            end
            qMetricsExist = ~isempty(dir(fullfile(savePath, '**', 'templates._bc_qMetrics.parquet'))); % ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) not used anymore?
            idx = sp.SessionID==id;
            InspectionFlag = 0;
            if ~qMetricsExist || Params.RedoQM
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(rawD(id).name,'cbin')) && Params.DecompressLocal
                    % detect whether data is compressed, decompress locally if necessary
                    if ~exist(fullfile(Params.tmpdatafolder,strrep(rawD(id).name,'cbin','bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                        decompDataFile = bc_extractCbinData(fullfile(rawD(id).folder,rawD(id).name),...
                            [], [], 0,  fullfile(Params.tmpdatafolder,strrep(rawD(id).name,'cbin','bin')));
                        copyfile(strrep(fullfile(rawD(id).folder,rawD(id).name),'cbin','meta'),strrep(fullfile(Params.tmpdatafolder,rawD(id).name),'cbin','meta'))
                    end
                    DecompressionFlag = 1;
                    if Params.InspectQualityMetrics
                        InspectionFlag=1;
                    end
                end

                paramBC.rawFile = fullfile(Params.tmpdatafolder,strrep(rawD(id).name,'cbin','bin'));
                paramBC.ephysMetaFile = (strrep(fullfile(Params.tmpdatafolder,rawD(id).name),'cbin','meta'));

                %             idx = ismember(sp.spikeTemplates,clusidtmp(Good_IDtmp)); %Only include good units
                %careful; spikeSites zero indexed
                [qMetric, unitType] = bc_runAllQualityMetrics(paramBC, sp.st(idx)*sp.sample_rate, sp.spikeTemplates(idx)+1, ...
                    sp.temps, sp.tempScalingAmps(idx),sp.pcFeat(idx,:,:),sp.pcFeatInd+1,channelpostmp, savePath); % Be careful, bombcell needs 1-indexed!

            else
                paramBC.rawFile = fullfile(Params.tmpdatafolder,strrep(rawD(id).name,'cbin','bin'));
                if exist('ephysap_tmp', 'var')
                    paramBC.tmpFolder = [ephysap_tmp, '/..'];
                end
                d = dir(fullfile(savePath, '**', 'templates._bc_qMetrics.parquet'));
                qMetricsPath = d.folder;
                [~, qMetric, fractionRPVs_allTauR] = bc_loadSavedMetrics(qMetricsPath);
                unitType = bc_getQualityUnitType(paramBC, qMetric);
%                 unitType(:) = 1; ???
                
                %{ 
                % Commmented by CB for now    
                spike_templates_0idx = readNPY([savePath filesep 'spike_templates.npy']);
                spikeTemplates = spike_templates_0idx + 1;
                uniqueTemplates = unique(spikeTemplates);
                % need to load forGUI.tempWv??
                bc_plotGlobalQualityMetric(qMetric, paramBC, unitType, uniqueTemplates, forGUI.tempWv);
                %}

                %                 load(fullfile(savePath, 'qMetric.mat'))
            end
            unitTypeAcrossRec{id} = unitType;
            theseuniqueTemplates{id} = unique(sp.spikeTemplates(idx));

            if InspectionFlag % Doesn't currently work: Julie will update bombcell
                bc_loadMetricsForGUI
                %% Inspect the results?
                % GUI guide:
                % left/right arrow: toggle between units
                % g : go to next good unit
                % m : go to next multi-unit
                % n : go to next noise unit
                % up/down arrow: toggle between time chunks in the raw data
                % u: brings up a input dialog to enter the unit you want to go to
                unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
                    param, probeLocation, unitType, loadRawTraces);
                %                 bc_unitQualityGUI(memMapData, ephysData, qMetric, rawWaveforms, paramBC,...
                %                     probeLocation, unitType, plotRaw);
                disp('Not continuing until GUI closed')
                while isvalid(unitQualityGuiHandle) && Params.InspectQualityMetrics
                    pause(0.01)
                end
                disp('GUI closed')
            end
            %% Apply QM findings:
            disp('Only use units that passed the quality metrics parameters')
            % This is necessary inside this loop when running stitched - don't change!!
            Good_IDtmp = nan(1,length(unitType)); % Replace with unitType
            Good_IDtmp(unitType ~= 1) = 0; % MUA
            Good_IDtmp(unitType == 1) = 1; % Good
            Good_IDtmp(unitType == 0) = 0;
            %         NoiseUnit = false(size(Good_IDtmp));
            %         NoiseUnit(unitType == 0)=1; % NOISE

           
            Good_ID = [Good_ID,Good_IDtmp]; %Identify good clusters

        end
       
        AllUniqueTemplates = cat(1,AllUniqueTemplates(:),cat(1,theseuniqueTemplates{:}));

    else
        AllUniqueTemplates = cat(1,AllUniqueTemplates,unique(sp.spikeTemplates));
        Good_ID = [Good_ID,Good_IDtmp]; %Identify good clusters

        %         NoiseUnit = false(size(Good_IDtmp));
    end
    addthis = nanmax(recsesAll);
    if isempty(addthis)
        addthis=0;
    end
    if exist('theseuniqueTemplates','var') && iscell(theseuniqueTemplates)
        recsesAlltmp = arrayfun(@(X) repmat(addthis+X,1,length(theseuniqueTemplates{X})),[1:length(theseuniqueTemplates)],'UniformOutput',0);
        recsesAll =cat(1,recsesAll(:), cat(2,recsesAlltmp{:})');
    else
        recsesAll = cat(1,recsesAll(:),repmat(addthis+1,1,length(Good_IDtmp))');
    end
    sp.RecSes = sp.SessionID+countid-1; %Keep track of recording session, as cluster IDs are not unique across sessions

    %% Compute cross-correlation matrices
    Good_Idx = find(Good_ID); % Only care about good units at this point
    
    % Define edges for this dataset
    edges = floor(min(sp.st))-Params.binsz/2:Params.binsz:ceil(max(sp.st))+Params.binsz/2;

    % bin data to create PSTH
    sr = nan(numel(Good_Idx),numel(edges)-1);
    for uid = 1:numel(Good_Idx)
        sr(uid,:) =  histcounts(sp.st(sp.spikeTemplates == sp.cids(Good_Idx(uid))),edges);
    end

    % Define folds (two halves)
    idx_fold1 = 1:floor(size(sr,2)./2);
    idx_fold2 = floor(size(sr,2)./2)+1:floor(size(sr,2)./2)*2;

    % Find cross-correlation in first and second half of session
    SessionCorrelations.fold1 = corr(sr(:,idx_fold1)',sr(:,idx_fold1)')';
    SessionCorrelations.fold2 = corr(sr(:,idx_fold2)',sr(:,idx_fold2)')';

    % Nan the diagonal
    SessionCorrelations.fold1(logical(eye(size(SessionCorrelations.fold1)))) = nan;
    SessionCorrelations.fold2(logical(eye(size(SessionCorrelations.fold2)))) = nan;

    %% Save out sp and clusinfo for this session in the correct folder
    ShankOpt = unique(Shank);
    ShankID = nan(size(Shank));
    for shankid = 1:length(ShankOpt)
        ShankID(Shank==ShankOpt(shankid))=shankid;
    end
    clusinfo.Shank = Shank;
    clusinfo.ShankID = ShankID;
    clusinfo.RecSesID = recsesAll;
    clusinfo.ch = channel;
    clusinfo.depth = depth;
    clusinfo.cluster_id = AllUniqueTemplates;
    clusinfo.group = Label;
    clusinfo.Good_ID = Good_ID;
    % clusinfo.Noise_ID = NoiseUnit;

    sp.sample_rate = sp.sample_rate(1);

    % Clean up unnecessary information in sp
    sp = rmfield(sp,'tempScalingAmps');
    sp = rmfield(sp,'temps');
    sp = rmfield(sp,'winv');
    sp = rmfield(sp,'pcFeat');
    sp = rmfield(sp,'pcFeatInd');
    sp = rmfield(sp,'tempAmps');
    sp = rmfield(sp,'tempsUnW');
    sp = rmfield(sp,'templateDuration');
    sp = rmfield(sp,'waveforms');

    save(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'clusinfo','Params','SessionCorrelations','-v7.3')
    if Params.saveSp
        save(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'sp','-append','-v7.3')
    end

    countid=countid+1;
    close all
end

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(KiloSortPaths));
countid=1;
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
        continue
    end
   
    disp(['Loading clusinfo for ' KiloSortPaths{subsesid}])
    tmp = matfile(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'));
    clusinfo{subsesid} = tmp.clusinfo;
    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = repmat(countid,size(clusinfo{subsesid}.RecSesID));
    countid=countid+1;
end

% Add all cluster information in one 'cluster' struct - can be used for further analysis
clusinfo = [clusinfo{:}];
clusinfoNew = struct;
fields = fieldnames(clusinfo(1));
for fieldid=1:length(fields)
    try
        eval(['clusinfoNew.' fields{fieldid} '= cat(1,clusinfo(:).' fields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.' fields{fieldid} '= cat(2,clusinfo(:).' fields{fieldid} ');'])
        catch ME
            keyboard
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

%% UnitMatch Parameters
% Use some bombcell parameters
if Params.RunQualityMetrics
    paramBC = bc_qualityParamValuesForUnitMatch;
    UMparam.sampleamount = paramBC.nRawSpikesToExtract; %500; % Nr. waveforms to include
    UMparam.spikeWidth =paramBC.spikeWidth; %82; % in sample space (time)
else
    UMparam.sampleamount = 200; %500; % Nr. waveforms to include
    UMparam.spikeWidth = 82; %82; % in sample space (time)
end

UMparam.RunPyKSChronicStitched = Params.RunPyKSChronicStitched;
UMparam.SaveDir = fullfile(Params.SaveDir,'UnitMatch');
if ~isdir(UMparam.SaveDir)
    mkdir(UMparam.SaveDir)
end
UMparam.ACGbinSize =  0.001;%paramBC.ACGbinSize;
UMparam.ACGduration = 1;%paramBC.ACGduration;

UMparam.UseBombCelRawWav = Params.RunQualityMetrics; % 1 by default
UMparam.KSDir = AllKiloSortPaths;
UMparam.channelpos = AllChannelPos;
UMparam.AllRawPaths = RawDataPaths;
UMparam.AllDecompPaths = arrayfun(@(X) fullfile(Params.tmpdatafolder,strrep(RawDataPaths(X).name,'cbin','bin')),1:length(RawDataPaths),'Uni',0);
UMparam.RedoExtraction = 0; % Only necessary if KS was redone!
UMparam.ProbabilityThreshold = 0.5;
UMparam.binsize = Params.binsz;
if Params.UnitMatch
    UnitMatchExist = dir(fullfile(UMparam.SaveDir,'UnitMatch.mat'));
    if ~isempty(UnitMatchExist) && ~Params.RedoUnitMatch
        load(fullfile(UMparam.SaveDir,'UnitMatch.mat'))
        clusinfo.UniqueID = UniqueID;
        clusinfo.MatchTable = MatchTable;
    else
        % Need to decompress if decompression wasn't done yet
        for id = 1:length(RawDataPaths)
            ephysap_tmp = [];
            ephysap_path = fullfile(RawDataPaths(id).folder,RawDataPaths(id).name);
            if (isempty(UnitMatchExist) || Params.RedoUnitMatch) && ~DecompressionFlag && ~UMparam.UseBombCelRawWav
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(RawDataPaths(id).name,'cbin')) && Params.DecompressLocal
                    if ~exist(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')

                        decompDataFile = bc_extractCbinData(fullfile(RawDataPaths(id).folder,RawDataPaths(id).name),...
                            [], [], 0, fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')));
                        paramBC.rawFile = decompDataFile;
                        % Decompression
%                         success = pyrunfile("MTSDecomp_From_Matlab.py","success",datapath = strrep(fullfile(RawDataPaths(id).folder,RawDataPaths(id).name),'\','/'),...
%                             JsonPath =  strrep(fullfile(RawDataPaths(id).folder,strrep(RawDataPaths(id).name,'cbin','ch')),'\','/'), savepath = strrep(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')),'\','/'));
%                         % Also copy metafile
                        copyfile(strrep(fullfile(RawDataPaths(id).folder,RawDataPaths(id).name),'cbin','meta'),strrep(fullfile(Params.tmpdatafolder,RawDataPaths(id).name),'cbin','meta'))
                    end
                    ephysap_tmp = fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin'));
%                     DecompressionFlag = 1;
                end
            end
        end

        % Run UnitMatch
        [UniqueIDConversion, MatchTable, WaveformInfo, AllSessionCorrelations, UMparam] = UnitMatch(clusinfo,UMparam);
        save(fullfile(UMparam.SaveDir,'UnitMatch.mat'),'UniqueIDConversion','MatchTable','WaveformInfo','AllSessionCorrelations','UMparam')
        clusinfo.UniqueID = UniqueIDConversion.UniqueID;
        clusinfo.MatchTable = MatchTable;
    end
elseif DecompressionFlag % You might want to at least save out averaged waveforms for every session to get back to later, if they were saved out by bomcell
    % Extract average waveforms
    ExtractAndSaveAverageWaveforms(clusinfo,UMparam)
end

%% Here we're going to actually load in all the sessions requested - sp
sp = cell(1,length(KiloSortPaths));
if Params.saveSp
    countid=1;
    for subsesid=1:length(KiloSortPaths)
        if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
            continue
        end

        disp(['Loading spike data for ' KiloSortPaths{subsesid}])
        tmp = matfile(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'));
        sp{subsesid} = tmp.sp;
        % Replace recsesid with subsesid
        sp{subsesid}.RecSes = repmat(countid,size(sp{subsesid}.RecSes));
        countid=countid+1;
    end
    % Add all spikedata in one spikes struct - can be used for further analysis
    sp = [sp{:}];
    spnew = struct;

    fields = fieldnames(sp(1));
    for fieldid=1:length(fields)
        try
            eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
        catch ME
            if strcmp(ME.message,'Out of memory.')
                eval(['spnew.' fields{fieldid} ' = sp(1).' fields{fieldid} ';'])
                for tmpid = 2:length(sp)
                    eval(['spnew.' fields{fieldid} ' = cat(1,spnew.' fields{fieldid} ', sp(tmpid).' fields{fieldid} ');'])
                end
            else
                eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])
            end
        end
    end
    sp = spnew;
    sp.sample_rate = sp.sample_rate(1);
    clear spnew

    % Save out sp unique cluster
    nclus = length(clusinfo.cluster_id);
    sp.UniqClu = sp.clu;
    if Params.UnitMatch
        for clusid=1:nclus
            sp.UniqClu(sp.clu==clusinfo.cluster_id(clusid) & sp.RecSes==clusinfo.RecSesID(clusid)) = clusinfo.UniqueID(clusid);
        end
    end
end

%% Remove temporary files
if 0%Params.DecompressLocal && DecompressionFlag
    clear memMapData
    clear ap_data
    try
        for id = 1:length(RawDataPaths)
            delete(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')))
            delete(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','meta')))

            % Bombcell takes up a lot of space: delete
            delete(fullfile(Params.tmpdatafolder,['RawWaveforms_' strrep(RawDataPaths(id).name,'cbin','bin')],'*'))
            rmdir(fullfile(Params.tmpdatafolder,['RawWaveforms_' strrep(RawDataPaths(id).name,'cbin','bin')]))
        end
    catch ME
        keyboard
    end
end

return

% % Find correct dataset index
% spikeTimes =cat(1,sp(:).st);
% spikeRecSes = cat(1,sp(:).RecSes);
% spikeSites = cat(1,sp(:).spikeTemplates);
% spikeCluster = cat(1,sp(:).clu);
% spikeAmps = cat(1,sp(:).spikeAmps);
% spikeDepths = cat(1,sp(:).spikeDepths);

% templateDepths = cat(1,sp(:).templateDepths);
% tempAmps = cat(1,sp(:).tempAmps);
% tempsUnW = cat(1,sp(:).tempsUnW);
% templateDuration = cat(1,sp(:).templateDuration);
% waveforms = cat(1,sp(:).waveforms);
% templateWaveforms = cat(1,sp(:).temps);
% templateAmplitudes = cat(1,sp(:).tempScalingAmps);
% pcFeatures = cat(1,sp(:).pcFeat);
% pcFeatureIdx = cat(1,sp(:).pcFeatInd);

