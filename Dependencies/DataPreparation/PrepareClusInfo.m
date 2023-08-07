function Params = PrepareClusInfo(KiloSortPaths, Params, RawDataPaths)
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
    if nargin < 2
        disp('No params given. Use default - although this is not advised...')
        Params.loadPCs = 1;
        Params.RunPyKSChronicStitched = 1;
        Params.DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
        Params.RedoQM = 0; %if 1, redo quality metrics if it already exists
        Params.RunQualityMetrics = 1; % If 1, Run the quality metrics
        Params.InspectQualityMetrics = 0; % Inspect the quality metrics/data set using the GUI
        Params.UnitMatch = 1; % Matching chronic recording using QM instead of using pyks chronic output
        Params.RedoUnitMatch = 0; % Redo unitmatch
        Params.SaveDir = KiloSortPaths(1); % Directory to save QM and UnitMatch results
        Params.tmpdatafolder = KiloSortPaths(1); %
        Params.binsz = 0.01; % Binsize in time (s) for the cross-correlation fingerprint. We recommend ~2-10ms time windows
        Params.saveSp = 0;
        Params.extractSync = 0;
        Params.deNoise = 1;
        Params.nSavedChans = 385;
        Params.nSyncChans = 1;
    end
    if Params.RunQualityMetrics
        Params.loadPCs = 1; %If you want to run QM you need this
    end

    if nargin < 3
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
channelmap = [];
channelpos = [];

AllKiloSortPaths = cell(1, 0);
AllChannelPos = cell(1, 0);
countid = 1;
% figure;
cols = jet(length(KiloSortPaths));
for subsesid = 1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid}, '*.npy')))
        continue
    end

    %% initialize for new round
    Good_ID = [];
    Shank = [];
    ProbeAll = [];
    recsesAll = [];
    channel = [];
    Label = [];
    depth = [];
    AllUniqueTemplates = [];
    cluster_id = [];
    recses = [];
    if any(strfind(KiloSortPaths{subsesid}, 'Probe'))
        probeid = str2num(KiloSortPaths{subsesid}(strfind(KiloSortPaths{subsesid}, 'Probe') + 5));
    else
        disp('Probe ID unknown')
        probeid = 0;
    end

    %% save data paths information
    if Params.RunPyKSChronicStitched %CALL THIS STITCHED --> Only works when using RunPyKS2_FromMatlab as well from this toolbox
        if UseParamsKS
            try
                spikeStruct = loadParamsPy(fullfile(KiloSortPaths{subsesid}, 'params.py'));
                rawD = spikeStruct.dat_path;
                rawD = strsplit(rawD, ',');
                for rid = 1:length(rawD)
                    rawD{rid} = rawD{rid}(strfind(rawD{rid}, '"') + 1:end);
                    rawD{rid} = rawD{rid}(1:strfind(rawD{rid}, '"') - 1);
                    rawD{rid} = dir(rawD{rid});
                    if isempty(rawD{rid})
                        rawD{rid} = dir(strrep(rawD{rid}, 'bin', 'cbin'));
                    end
                end
                rawD = cat(2, rawD{:});
            catch
                sessionsIncluded = dir(fullfile(KiloSortPaths{subsesid}, 'SessionsIncluded.mat'));
                sessionsIncluded = load(fullfile(sessionsIncluded.folder, sessionsIncluded.name));
                rawD = arrayfun(@(X) dir(sessionsIncluded.ThesePaths{X}), 1:length(sessionsIncluded.ThesePaths), 'UniformOutput', 0);
                rawD = cat(2, rawD{:});
            end
        else
            rawD = dir(fullfile(RawDataPaths(subsesid).folder, RawDataPaths(subsesid).name));
        end

        RawDataPaths = rawD;
        AllKiloSortPaths = cat(2, AllKiloSortPaths{:}, repmat(KiloSortPaths(subsesid), 1, length(rawD)));
    else
        if UseParamsKS
            spikeStruct = loadParamsPy(fullfile(KiloSortPaths{subsesid}, 'params.py'));
            rawD = spikeStruct.dat_path;
            rawD = rawD(strfind(rawD, '"')+1:end);
            rawD = rawD(1:strfind(rawD, '"')-1);
            tmpdr = rawD;
            rawD = dir(rawD);
            if isempty(rawD)
                rawD = dir(strrep(tmpdr, 'bin', 'cbin'));
            end


            % Save for later
            RawDataPaths(subsesid) = rawD;

            if isempty(rawD)
                disp('Bug...?')
                keyboard
            end
        else
            if isstruct(RawDataPaths)
                rawD = dir(fullfile(RawDataPaths(subsesid).folder, RawDataPaths(subsesid).name));

            else
                rawD = dir(fullfile(RawDataPaths{subsesid}));
            end

        end
        AllKiloSortPaths = {AllKiloSortPaths{:}, KiloSortPaths{subsesid}};
    end
    DecompressionFlag = 0;

    %% Channel data
    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'channel_map.npy'));
    channelmaptmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));

    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'channel_positions.npy'));
    channelpostmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));
    if length(channelmaptmp) < length(channelpostmp)
        channelmaptmp(end+1:length(channelpostmp)) = length(channelmaptmp):length(channelpostmp) - 1;
    end
   

    %% Is it correct channelpos though...? Check using raw data
    channelpostmpconv = ChannelIMROConversion(rawD(1).folder, 0); % For conversion when not automatically done
    AllChannelPos{countid} = channelpostmpconv;

    %% Load existing?
    if exist(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat')) && ~Params.RedoQM && ~Params.ReLoadAlways
        % Check if parameters are the same, of not we have to redo it
        % anyway
        tmpparam = matfile(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'));
        tmpparam = tmpparam.Params;

        if tmpparam.RunQualityMetrics == Params.RunQualityMetrics && tmpparam.RunPyKSChronicStitched == Params.RunPyKSChronicStitched && ...
                tmpparam.separateIMRO == Params.separateIMRO
            disp(['Found existing data in ', KiloSortPaths{subsesid}, ', Using this...'])
            %         load(fullfile(Params.SaveDir,'PreparedData.mat'))

            % Use UnitMatch Output if available
            %         if Params.UnitMatch
            %             disp('Using UnitMatch Clusters!')
            %             sp.clu = sp.UniqClu; %Temporary replace for rest of code
            %             clusinfo.cluster_id = clusinfo.UniqueID;
            %         end
            countid = countid + 1;
            continue
        end
    end

    %% Load histology if available
    tmphisto = dir(fullfile(KiloSortPaths{subsesid}, 'HistoEphysAlignment.mat'));
    clear Depth2AreaPerUnit
    if ~isempty(tmphisto)
        histodat = load(fullfile(tmphisto.folder, tmphisto.name));
        Depth2AreaPerUnit = histodat.Depth2AreaPerUnit;
    end

    %% Load Spike Data
    sp = loadKSdir(fullfile(KiloSortPaths{subsesid}), Params); % Load Spikes with PCs
    [sp.spikeAmps, sp.spikeDepths, sp.templateDepths, sp.templateXpos, sp.tempAmps, sp.tempsUnW, sp.templateDuration, sp.waveforms] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.xcoords, sp.spikeTemplates, sp.tempScalingAmps); %from the spikes toolbox
    templateWaveforms = sp.temps;
    %% Remove noise; spikes across all channels'
    if ~isfield(Params,'deNoise')
        Params.deNoise = 1;
    end
    if Params.deNoise == 1
        sp = RemoveNoiseAmplitudeBased(sp);
    end

    %% Bombcell parameters
    % clear paramBC
    paramBC = bc_qualityParamValuesForUnitMatch(dir(strrep(fullfile(rawD(1).folder, rawD(1).name), 'cbin', 'meta')), fullfile(Params.tmpdatafolder, strrep(rawD(1).name, 'cbin', 'bin')));

    %% Load Cluster Info
    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'cluster_info.tsv')); % If you did phy (manual curation) we will find this one... We can trust you, right?
    if isempty(myClusFile)
        disp('This data is not curated with phy! Hopefully you''re using automated quality metrics to find good units!')
        curratedflag = 0;
        myClusFile = dir(fullfile(KiloSortPaths{subsesid},'cluster_group.tsv'));
        if isempty(myClusFile)
            clusinfo = tdfread(fullfile(KiloSortPaths{subsesid}, 'cluster_KSLabel.tsv'));
        else
            clusinfo = tdfread(fullfile(myClusFile(1).folder, myClusFile(1).name));
        end
        
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after splitting/merging)
        [clusinfo,sp,emptyclus] = RemovingEmptyClusters(clusinfo,sp);

        %clusidtmp = clusinfo.cluster_id;
        clusinfo.cluster_id = unique(sp.spikeTemplates);
        clusidtmp = clusinfo.cluster_id;
        tmpLabel = char(length(clusinfo.cluster_id));
        KSLabelfile = tdfread(fullfile(KiloSortPaths{subsesid}, 'cluster_KSLabel.tsv'));
        tmpLabel(ismember(clusinfo.cluster_id, KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id, clusinfo.cluster_id));
        Label = [Label, tmpLabel];
        totSpkNum = histc(sp.clu, sp.cids);
        Good_IDtmp = ismember(tmpLabel, 'g') & totSpkNum' > paramBC.minNumSpikes; %Identify good clusters

        % Find depth and channel
        depthtmp = nan(length(clusinfo.cluster_id), 1);
        xtmp = nan(length(clusinfo.cluster_id), 1);
        channeltmp = nan(length(clusinfo.cluster_id), 1);
        for clusid = 1:length(depthtmp)
            % Depth according to PyKS2 output
            depthtmp(clusid) = round(sp.templateDepths(clusid)); %round(nanmean(sp.spikeDepths(find(sp.clu==clusid-1))));
            xtmp(clusid) = sp.templateXpos(clusid);
            [~, minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1, channelpostmp(X, :), [xtmp(clusid), depthtmp(clusid)]), 'euclidean'), 1:size(channelpostmp, 1), 'UniformOutput', 0)));
            try
                channeltmp(clusid) = channelmaptmp(minidx);
                depthtmp(clusid) = channelpostmpconv(minidx, 2);
                xtmp(clusid) = channelpostmpconv(minidx, 1);
            catch
                channeltmp(clusid) = channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmpconv(minidx-1, 2);
                xtmp(clusid) = channelpostmpconv(minidx-1, 1);
            end
            sp.spikeDepths(ismember(sp.clu, clusidtmp(clusid))) = depthtmp(clusid);

        end
        depth = cat(1, depth, depthtmp);
        channel = cat(1, channel, channeltmp);

    else
        disp('You did manual curation. You champion. If you have not enough time, maybe consider some automated algorithm...')
        CurationDone = 1;
        save(fullfile(KiloSortPaths{subsesid}, 'CuratedResults.mat'), 'CurationDone')
        clusinfo = tdfread(fullfile(myClusFile(1).folder, myClusFile(1).name));
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after  splitting/merging)
        [clusinfo, sp, emptyclus] = RemovingEmptyClusters(clusinfo,sp);
       
        curratedflag = 1;
        if isfield(clusinfo, 'id')
            clusidtmp = clusinfo.id;
            cluster_id = [cluster_id, clusinfo.id];
        elseif isfield(clusinfo, 'cluster_id')
            clusidtmp = clusinfo.cluster_id;

            cluster_id = [cluster_id, clusinfo.cluster_id];
        else
            keyboard
            disp('Someone thought it was nice to change the name again...')
        end
        KSLabel = clusinfo.KSLabel;
        Label = [Label, clusinfo.group]; % You want the group, not the KSLABEL!
        % Find depth and channel
        depthtmp = nan(length(clusidtmp), 1);
        xtmp = nan(length(clusidtmp), 1);
        channeltmp = nan(length(clusidtmp), 1);
        for clusid = 1:length(depthtmp)
            % Depth according to PyKS2 output
            depthtmp(clusid) = round(sp.templateDepths(clusid)); %round(nanmean(sp.spikeDepths(find(sp.clu==clusid-1))));
            xtmp(clusid) = sp.templateXpos(clusid);
            [~, minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1, channelpostmp(X, :), [xtmp(clusid), depthtmp(clusid)]), 'euclidean'), 1:size(channelpostmp, 1), 'UniformOutput', 0)));
            try
                channeltmp(clusid) = channelmaptmp(minidx);
                depthtmp(clusid) = channelpostmpconv(minidx, 2);
                xtmp(clusid) = channelpostmpconv(minidx, 1);
            catch
                channeltmp(clusid) = channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmpconv(minidx-1, 2);
                xtmp(clusid) = channelpostmpconv(minidx-1, 1);
            end
            sp.spikeDepths(ismember(sp.clu, clusidtmp(clusid))) = depthtmp(clusid);

        end

        depth = [depth, depthtmp];

        channel = [channel, clusinfo.ch];
        Good_IDtmp = ismember(cellstr(clusinfo.group), 'good');
        channeltmp = clusinfo.ch;
    end
    channelpostmp = channelpostmpconv;
    ypostmp = channelpostmp(:, 2);
    xpostmp = channelpostmp(:, 1);
    xposopt = (floor(xpostmp./250)); % Assuming no new shank if not at least 100 micron apart
    Shanktmp = floor(xpostmp(channeltmp+1)./250);

    %     [~,minid] = arrayfun(@(X) (abs(floor(xpostmp(X)./250)-xposopt)),channeltmp+1,'UniformOutput',0);
    Shank = cat(1, Shank, Shanktmp);

    recses = cat(1, recses, repmat(countid, length(channeltmp), 1));

    %     channelrecses = cat(2,channelrecses,repmat(countid,1,size(channelpostmp,1)));
    %     channelpos = cat(1,channelpos,channelpostmp);
    channelpos = channelpostmp; %only every need this here for mulitple sessions when chronic.
    channelmap = channelmaptmp;

    %     scatter(xtmp,depthtmp,10,cols(countid,:))
    %     hold on
    %     drawnow

    if Params.RunQualityMetrics
        theseuniqueTemplates = [];
        unitTypeAcrossRec = [];

        %% Quality metrics - Bombcell (https://github.com/Julie-Fabre/bombcell)
        for id = 1:length(rawD)
            ephysap_tmp = [];

            ephysap_path = fullfile(rawD(id).folder, rawD(id).name);

                savePath = fullfile(KiloSortPaths{subsesid});

            qMetricsExist = ~isempty(dir(fullfile(savePath, '**', 'templates._bc_qMetrics.parquet'))); % ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) not used anymore?
            idx = sp.SessionID == 1;
            InspectionFlag = 0;
            if ~isempty(dir(fullfile(savePath, '**', 'RawWaveforms'))) % if raw waveforms have not been extract, decompress data for extraction
                disp('Extracting sync file...')
                % detect whether data is compressed, decompress locally if necessary
                if ~exist(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')))
                    disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                    decompDataFile = bc_extractCbinData(fullfile(rawD(id).folder, rawD(id).name), ...
                        [], [], 0, fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')));
                    status = copyfile(strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'), strrep(fullfile(Params.tmpdatafolder, rawD(id).name), 'cbin', 'meta')); %QQ doesn't work on linux
                end
                DecompressionFlag = 1;
                 if status == 0 %could not copy meta file - use original meta file 
                     [Imecmeta] = ReadMeta2(fullfile(rawD(id).folder, strrep(rawD(id).name, 'cbin', 'meta')), 'ap');
                 else
                    [Imecmeta] = ReadMeta2(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'meta')), 'ap');
                end
                nchan = strsplit(Imecmeta.acqApLfSy, ',');
                nChansInFile = str2num(nchan{1}) + str2num(nchan{3});

                syncDatImec = extractSyncChannel(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')), nChansInFile, nChansInFile); %Last channel is sync
                status = copyfile(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, '.cbin', '_sync.dat')), fullfile(rawD(id).folder, strrep(rawD(id).name, '.cbin', '_sync.dat'))); %QQ doesn't work on linux

            end
            if ~qMetricsExist || Params.RedoQM
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(rawD(id).name, 'cbin')) && Params.DecompressLocal
                    % detect whether data is compressed, decompress locally if necessary
                    if ~exist(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                        decompDataFile = bc_extractCbinData(fullfile(rawD(id).folder, rawD(id).name), ...
                            [], [], 0, fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')));
                        status = copyfile(strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'), strrep(fullfile(Params.tmpdatafolder, rawD(id).name), 'cbin', 'meta')); %QQ doesn't work on linux
                    end

                    DecompressionFlag = 1;
                    if Params.InspectQualityMetrics
                        InspectionFlag = 1;
                    end
                end

                paramBC.rawFile = fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'));
                paramBC.ephysMetaFile = (strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'));

                %             idx = ismember(sp.spikeTemplates,clusidtmp(Good_IDtmp)); %Only include good units
                %careful; spikeSites zero indexed
                [qMetric, unitType] = bc_runAllQualityMetrics(paramBC, sp.st(idx)*sp.sample_rate, sp.spikeTemplates(idx)+1, ...
                    templateWaveforms, sp.tempScalingAmps(idx), sp.pcFeat(idx, :, :), sp.pcFeatInd+1, channelpostmp, savePath); % Be careful, bombcell needs 1-indexed!

            else
                paramBC.rawFile = fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'));
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
            qMetricclusterID = qMetric.clusterID;

            unitTypeAcrossRec{id} = unitType;
            theseuniqueTemplates{id} = unique(sp.spikeTemplates);

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
            Good_IDtmp = nan(1, length(Good_IDtmp)); % Replace with unitType
            Good_IDtmp(qMetricclusterID(unitType ~= 1)) = 0; % MUA
            Good_IDtmp(qMetricclusterID(unitType == 1)) = 1; % Good
            Good_IDtmp(qMetricclusterID(unitType == 0)) = 0; % MUA
            Good_IDtmp(isnan(Good_IDtmp)) = 0; % Noise
            %         NoiseUnit = false(size(Good_IDtmp));
            %         NoiseUnit(unitType == 0)=1; % NOISE


            Good_ID = [Good_ID, Good_IDtmp]; %Identify good clusters

        end

        AllUniqueTemplates = cat(1, AllUniqueTemplates(:), cat(1, theseuniqueTemplates{:}));

    else
        AllUniqueTemplates = cat(1, AllUniqueTemplates, unique(sp.spikeTemplates));
        Good_ID = [Good_ID, Good_IDtmp]; %Identify good clusters

        %         NoiseUnit = false(size(Good_IDtmp));
    end
    addthis = nanmax(recsesAll);
    if isempty(addthis)
        addthis = 0;
    end
    if exist('theseuniqueTemplates', 'var') && iscell(theseuniqueTemplates)
        recsesAlltmp = arrayfun(@(X) repmat(addthis+X, 1, length(theseuniqueTemplates{X})), [1:length(theseuniqueTemplates)], 'UniformOutput', 0);
        recsesAll = cat(1, recsesAll(:), cat(2, recsesAlltmp{:})');
        ProbeAll = cat(1, ProbeAll(:), repmat(probeid, 1, length(cat(2, recsesAlltmp{:})))');
    else
        recsesAll = cat(1, recsesAll(:), repmat(addthis+1, 1, length(Good_IDtmp))');
        ProbeAll = cat(1, ProbeAll(:), repmat(probeid, 1, length(Good_IDtmp))');
    end
    sp.RecSes = sp.SessionID + countid - 1; %Keep track of recording session, as cluster IDs are not unique across sessions

    %% Save out sp and clusinfo for this session in the correct folder
    ShankOpt = unique(Shank);
    ShankID = nan(size(Shank));
    for shankid = 1:length(ShankOpt)
        ShankID(Shank == ShankOpt(shankid)) = shankid;
    end
    clusinfo.Shank = Shank;
    clusinfo.ProbeID = ProbeAll;
    clusinfo.ShankID = ShankID;
    clusinfo.RecSesID = recsesAll;
    clusinfo.ch = channel;
    clusinfo.depth = depth;
    clusinfo.cluster_id = AllUniqueTemplates;
    clusinfo.group = Label;
    if length(Good_ID)>length(recsesAll)
        NonEmptyIdx = true(1,length(recsesAll));
        NonEmptyIdx(emptyclus) = false;
        clusinfo.Good_ID = Good_ID(NonEmptyIdx);
    else
        clusinfo.Good_ID = Good_ID;
    end

    if exist('Depth2AreaPerUnit', 'var') % Add area information
        Idx = cell2mat(arrayfun(@(X) find(Depth2AreaPerUnit.Cluster_ID-1 == X), clusinfo.cluster_id, 'Uni', 0));
        clusinfo.Area = Depth2AreaPerUnit.Area(Idx);
        clusinfo.Coordinates = Depth2AreaPerUnit.Coordinates(Idx);
        if length(clusinfo.Coordinates) ~= length(clusinfo.cluster_id)
            clusinfo.Coordinates = []; % Doesn't work out
            clusinfo.Area = [];
        end
    end
    % clusinfo.Noise_ID = NoiseUnit;

    sp.sample_rate = sp.sample_rate(1);

    % Clean up unnecessary information in sp
    sp = rmfield(sp, 'tempScalingAmps');
    sp = rmfield(sp, 'temps');
    sp = rmfield(sp, 'winv');
    sp = rmfield(sp, 'pcFeat');
    sp = rmfield(sp, 'pcFeatInd');
    sp = rmfield(sp, 'tempAmps');
    sp = rmfield(sp, 'tempsUnW');
    sp = rmfield(sp, 'templateDuration');
    sp = rmfield(sp, 'waveforms');
    save(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'), 'clusinfo', 'Params', '-v7.3')
    if Params.saveSp %QQ not using savePaths
        try
            save(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'), 'sp', '-append')
        catch ME
            disp(ME)
        end
    end

    countid = countid + 1;
end

Params.AllChannelPos = AllChannelPos;
Params.RawDataPaths = RawDataPaths;
Params.DecompressionFlag = DecompressionFlag;

%% Remove temporary files
if 0 %Params.DecompressLocal && DecompressionFlag
    clear memMapData
    clear ap_data
    try
        for id = 1:length(RawDataPaths)
            delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths(id).name, 'cbin', 'bin')))
            delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths(id).name, 'cbin', 'meta')))
        end
    catch ME
        keyboard
    end
end

return
