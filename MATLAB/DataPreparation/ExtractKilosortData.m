function Params = ExtractKilosortData(KiloSortPaths, Params, RawDataPathsInput)

% Prepares cluster information for subsequent analysis
%% Inputs:
% KiloSortPaths = List of directories pointing at kilosort output (same format as what you get when
% calling 'dir') for all sessions you want to analyze as one session (e.g.
% chronic recordings with same IMRO table: give a list with all sessions)

% Params: default parameters available using DefaultParametersExtractKSData
% Params.loadPCs=0;
% Params.RunPyKSChronicStitched = 1
% Params.DecompressLocal = 1; %if 1, uncompress data first if it's currently compressed
% Params.RedoQM = 0; %if 1, redo quality metrics if it already exists
% Params.RunQualityMetrics = 1; % If 1, Run the quality metrics
% Params.InspectQualityMetrics =0; % Inspect the quality metrics/data set using the GUI
% Params.tmpdatafolder %Directory to temporarily decompress data --> must
% be large enough!

% RawDataPaths = list of same directories pointing at the raw ephys data
% directory. If this input is missing, this code will try to find the raw
% ephys data in the params.py file, first line (dat_path).

%% Check inputs
if ~iscell(KiloSortPaths) %isstruct(KiloSortPaths) || isfield(KiloSortPaths(1),'name') || isfield(KiloSortPaths(1),'folder')
    error('This is not a cell... give correct input please')
end
if nargin < 2
    disp('No params given. Use default - although this is not advised...')
    Params = struct;
end
Params = DefaultParametersExtractKSData(Params,KiloSortPaths);

if nargin>2 && exist('RawDataPathsInput')
    Params.RawDataPaths = RawDataPathsInput;
end
try
    if ~isfield(Params,'RawDataPaths') || isempty(Params.RawDataPaths) || length(Params.RawDataPaths)~=length(KiloSortPaths)
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

AllKiloSortPaths = cell(1, length(KiloSortPaths));
AllChannelPos = cell(1, length(KiloSortPaths));
AllProbeSN = cell(1, length(KiloSortPaths));
RawDataPaths = cell(1, length(KiloSortPaths));
IncludeThese = true(1,length(KiloSortPaths));

countid = 1;
% figure;
cols = jet(length(KiloSortPaths));
for subsesid = 1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
        IncludeThese(subsesid) = false;
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
    ExtractChannelMapThenContinue = 0;

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
            if isstruct(RawDataPathsInput)
                rawD = dir(fullfile(RawDataPathsInput(subsesid).folder, RawDataPathsInput(subsesid).name));
            else
                rawD = dir(fullfile(RawDataPathsInput{subsesid}));
            end
        end

        RawDataPaths{subsesid} = rawD; % Definitely save as cell
        AllKiloSortPaths{subsesid} = repmat(KiloSortPaths(subsesid), 1, length(rawD));
    else
        if UseParamsKS
            spikeStruct = loadParamsPy(fullfile(KiloSortPaths{subsesid}, 'params.py'));
            rawD = spikeStruct.dat_path;
            if any(strfind(rawD, '"'))
                rawD = rawD(strfind(rawD(1:10), '"')+1:end);
                rawD = rawD(1:strfind(rawD, '"')-1);
            end
            if any(strfind(rawD,'../'))
                rawD = rawD(strfind(rawD,'../')+3:end);
            end
            tmpdr = rawD;
            rawD = dir(strrep(rawD, '/', filesep));
            if isempty(rawD)
                rawD = dir(strrep(tmpdr, 'bin', 'cbin'));
            end
            if isempty(rawD)
                rawD = dir(strrep(tmpdr, 'cbin', 'bin'));
            end
            % Try another way
            if isempty(rawD)
                FolderParts = strsplit(KiloSortPaths{subsesid},{'pyKS','PyKS','kilosort2'});
                rawD = dir(fullfile(FolderParts{1},'**','*bin'));
                if length(rawD)>1
                    rawD = rawD(1);
                end
                %                 rawD = fullfile(rawD.folder,rawD.name);
            end
        else
            try
                if isstruct(Params.RawDataPaths)
                    rawD = Params.RawDataPaths(subsesid);
                elseif iscell((Params.RawDataPaths))
                    rawD = Params.RawDataPaths{subsesid};
                    if ~isstruct(rawD)
                        rawD = dir(rawD);
                    end
                else
                    rawD = dir(fullfile(Params.RawDataPaths{subsesid}));
                end
            catch ME
                keyboard
            end
        end
        if isempty(rawD)
            disp(['No raw data available for ' KiloSortPaths{subsesid}])
            rawD = [];
        end

        RawDataPaths{subsesid} = rawD; % Definitely save as cell
        AllKiloSortPaths{subsesid} = KiloSortPaths{subsesid};
    end
    Params.DecompressionFlag = 0;

    %% Channel data
    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'channel_map.npy'));
    channelmaptmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));

    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'channel_positions.npy'));
    channelpostmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));

    % Some KS versions do not store empty channels?
    channelpos = nan(Params.nSavedChans-Params.nSyncChans,2);
    channelpos(channelmaptmp+1,:) = channelpostmp;

    % Also make channelmap appropriate size
    channeltmpnew = nan(Params.nSavedChans-Params.nSyncChans,1);
    channeltmpnew(channelmaptmp+1) = channelmaptmp;

    channelmaptmp = channeltmpnew; % Replace
    % Usually can just fill the missing channels easily
    channelmaptmp = round(fillmissing(channelmaptmp,'linear'));
    try
        channelpostmp = fillmissingprobe(channelpos);
    catch ME
        warning('Your channel map has gaps, please check')
        channelpostmp = channelpos;
    end


    %% Load existing?
    if exist(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat')) && Params.ReLoadAlways~=1
        % Check if parameters are the same, of not we have to redo it
        % anyway
        dateViolationflag = 0;
        FileInfo = dir(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'));
        if isfield(Params,'FromDate') && datetime(FileInfo(1).date) < Params.FromDate && Params.ReLoadAlways==2
            dateViolationflag = 1;
        end

        try
            tmpparam = load(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'),'SessionParams');
            tmpparam = tmpparam.SessionParams;
            overwrite = 0;
        catch  % Old way of storing things, less flexible so overwrite
            tmpparam = load(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'),'Params');
            tmpparam = tmpparam.Params;
            overwrite = 1;
        end

        if tmpparam.RunQualityMetrics == Params.RunQualityMetrics && tmpparam.RunPyKSChronicStitched == Params.RunPyKSChronicStitched && ~dateViolationflag
            if ~isfield(tmpparam,'RecordingDuration') || tmpparam.RecordingDuration < Params.MinRecordingDuration
                if ~isempty(rawD) & ~contains(rawD(1).name,'.dat') & ~contains(rawD(1).name,'.raw')
                    [channelpostmpconv, probeSN, recordingduration] = ChannelIMROConversion(rawD(1).folder, 0); % For conversion when not automatically done
                    Params.RecordingDuration = recordingduration;
                    if recordingduration<Params.MinRecordingDuration
                        disp([KiloSortPaths{subsesid} ' recording too short, skip...'])
                        continue
                    end
                end
            end
            disp(['Found existing data in ', KiloSortPaths{subsesid}, ', Using this...'])

            if isfield(tmpparam,'AllProbeSN')
                if exist('channelpostmpconv','var') && (size(channelpostmp,1) ~= size(tmpparam.AllChannelPos{1},1) || any(tmpparam.AllChannelPos{1}(:) ~= channelpostmp(:)))

                    figure; scatter(tmpparam.AllChannelPos{1}(:,1),tmpparam.AllChannelPos{1}(:,2))
                    hold on; scatter(channelpostmpconv(:,1),channelpostmpconv(:,2))
                    scatter(channelpostmp(:,1),channelpostmp(:,2))
                    title('Blue = stored in prepared data, red = correct map, yellow = used by KS')

                    SpacingsKS = max(unique(diff(channelpostmp(:,1))));
                    SpacingsReal = max(unique(diff(channelpostmpconv(:,1))));
                    if (SpacingsReal - SpacingsKS) == 50
                        disp('Kilosort used different shank spacing for IMRO table. Not a big problem')
                    else
                        warning('Kilosort probably used wrong IMRO table. Consider reanalyzing')
                    end
                    tmpparam.AllChannelPos = {channelpostmp}; % Using KS map now
                    tmpparam.AllProbeSN = {probeSN};
                    overwrite = 1;
                end
                AllChannelPos{subsesid} = tmpparam.AllChannelPos{1};
                AllProbeSN{subsesid} = tmpparam.AllProbeSN{1};
                countid = countid + 1;
                if ~overwrite
                    continue
                end
            else
                ExtractChannelMapThenContinue = 1;
            end
        end
        if overwrite
            SessionParams = tmpparam; % But only store for this specific session - otherwise confusion!
            SessionParams.KSDir = KiloSortPaths(subsesid);
            SessionParams.RawDataPaths = RawDataPaths(subsesid);
            save(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'),'SessionParams', '-append')
            if ~ExtractChannelMapThenContinue
                continue
            end
        end    
    end


    %% Is it correct channelpos though...? Check using raw data. While reading this information, also extract recording duration and Serial number of probe
    if ~isempty(rawD) & (~contains(rawD(1).name,'.dat') & ~contains(rawD(1).name,'.raw'))
        [channelpostmpconv, probeSN, recordingduration] = ChannelIMROConversion(rawD(1).folder, 0); % For conversion when not automatically done
        if recordingduration<Params.MinRecordingDuration
            disp([KiloSortPaths{subsesid} ' recording too short, skip...'])
            continue
        end
        AllProbeSN{subsesid} = probeSN;
    else
        probeSN = '000000';
        AllProbeSN{subsesid} = probeSN;
    end
    AllChannelPos{subsesid} = channelpostmp; % Use what was used for KS regardless.

    if ExtractChannelMapThenContinue % Version compatibility
        countid = countid + 1;
        ExtractChannelMapThenContinue = 0;
        continue
    end


    %% Load histology if available
    tmphisto = dir(fullfile(KiloSortPaths{subsesid}, 'HistoEphysAlignment.mat'));
    clear Depth2Area
    if ~isempty(tmphisto)
        histodat = load(fullfile(tmphisto.folder, tmphisto.name));
        Depth2Area = histodat.Depth2Area;
    end

    %% Load Spike Data
    sp = loadKSdir_UM(KiloSortPaths{subsesid}, Params); % Load Spikes with PCs
    [sp.spikeAmps, sp.spikeDepths, sp.templateDepths, sp.templateXpos, sp.tempAmps, sp.tempsUnW, sp.templateDuration, sp.waveforms] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.xcoords, sp.spikeTemplates, sp.tempScalingAmps); %from the spikes toolbox
    % Align with bombcell
    templateWaveforms = zeros(size(sp.temps));
    for t = 1:size(templateWaveforms, 1)
        templateWaveforms(t, :, :) = squeeze(sp.temps(t, :, :)) * sp.winv;
    end

    %% Check for spike hole bug in kilosort
    stderiv = [diff(sp.st)];
    Idx = (stderiv>7/1000); % Retain spike times >7ms
    tmp = diff(sp.st(Idx))./2.18689567; % time difference in batches
    figure; hh=histogram(tmp,500,'binLimits',[0 10]); % batch size
    values = hh.Values;
    Edges = hh.BinEdges;
    Ratio = nansum(values(logical(Edges(2:end)==round(Edges(2:end)))))./nansum(values(~logical(Edges(2:end)==round(Edges(2:end)))));
    title(['Spike hole ratio batch edges versus non-batch edges: ' num2str(Ratio)])
    xlabel('Batch number')
    ylabel('Count of spike times void>7ms')

    if Ratio > 1
        warning('Redo Kilosort for this session -- spike holes detected')
        saveas(gcf,fullfile(KiloSortPaths{subsesid},['SpikeHoles.bmp']))
    end


    %% Remove noise; spikes across all channels'
    if ~isfield(Params, 'deNoise')
        Params.deNoise = 1;
    end
    if Params.deNoise == 1
        sp = RemoveNoiseAmplitudeBased(sp);
    end

    %% Bombcell parameters
    % clear paramBC
    if Params.RunQualityMetrics
        if any(strfind(KiloSortPaths{subsesid},'KS4')) || any(strfind(KiloSortPaths{subsesid},'kilosort4')) % Used KS4?
            if ~isempty(rawD)
                paramBC = bc.qm.qualityParamValuesForUnitMatch(dir(strrep(fullfile(rawD(1).folder, rawD(1).name), 'cbin', 'meta')), fullfile(Params.tmpdatafolder, strrep(rawD(1).name, 'cbin', 'bin')),KiloSortPaths{subsesid},[],4);
            else
                paramBC = bc.qm.qualityParamValuesForUnitMatch([], [],KiloSortPaths{subsesid},[],4);
            end
        else
            if ~isempty(rawD)
                paramBC = bc.qm.qualityParamValuesForUnitMatch(dir(strrep(fullfile(rawD(1).folder, rawD(1).name), 'cbin', 'meta')), fullfile(Params.tmpdatafolder, strrep(rawD(1).name, 'cbin', 'bin')),KiloSortPaths{subsesid},[],2);

            else
                paramBC = bc.qm.qualityParamValuesForUnitMatch([],[],KiloSortPaths{subsesid},[],2);
            end
        end
    else
        paramBC.minNumSpikes = 300;
    end
    %% Load Cluster Info
    myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'cluster_info.tsv')); % If you did phy (manual curation) we will find this one... We can trust you, right?
    if isempty(myClusFile) || Params.RunQualityMetrics
        if ~Params.RunQualityMetrics
            disp('This data is not curated with phy! Hopefully you''re using automated quality metrics to find good units!')
        end
        curratedflag = 0;
        myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'cluster_group.tsv'));
        if isempty(myClusFile)
            clusinfo = tdfread(fullfile(KiloSortPaths{subsesid}, 'cluster_KSLabel.tsv'));
        else
            clusinfo = tdfread(fullfile(myClusFile(1).folder, myClusFile(1).name));
        end

        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after splitting/merging)
        [clusinfo, sp, emptyclus] = RemovingEmptyClusters(clusinfo, sp);

        %clusidtmp = clusinfo.cluster_id;
        sp.clu = int32(sp.clu);
        % clusinfo.cluster_id = unique(sp.spikeTemplates);
        clusidtmp = int32(clusinfo.cluster_id);
        tmpLabel = char(length(clusinfo.cluster_id));
        KSLabelfile = tdfread(fullfile(KiloSortPaths{subsesid}, 'cluster_KSLabel.tsv'));
        tmpLabel(ismember(clusinfo.cluster_id, KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id, clusinfo.cluster_id));
        Label = [Label, tmpLabel];
        totSpkNum = histc(sp.clu, clusinfo.cluster_id);
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
                depthtmp(clusid) = channelpostmp(minidx, 2);
                xtmp(clusid) = channelpostmp(minidx, 1);
            catch
                channeltmp(clusid) = channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmp(minidx-1, 2);
                xtmp(clusid) = channelpostmp(minidx-1, 1);
            end
            sp.spikeDepths(ismember(sp.clu, clusidtmp(clusid))) = depthtmp(clusid);

        end
        depth = cat(1, depth, depthtmp);
        channel = cat(1, channel, channeltmp);

    else
        disp('You did manual curation. You champion. If you have not enough time, maybe consider some automated algorithm...')
        CurationDone = 1;
        save(fullfile(KiloSortPaths{subsesid}, 'CuratedResults.mat'), 'CurationDone')
        myClusFile = dir(fullfile(KiloSortPaths{subsesid}, 'cluster_group.tsv'));
        clusinfo = tdfread(fullfile(myClusFile(1).folder, myClusFile(1).name));
        if ~isfield(clusinfo,'cluster_id')
            clusinfo.cluster_id = clusinfo.id;
        end
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after  splitting/merging)
        [clusinfo, sp, emptyclus] = RemovingEmptyClusters(clusinfo, sp);
        if ~any(~isnan(clusinfo.group))
            disp('clusinfo.group is empty, taking KS labels')
            clusinfo.group = clusinfo.KSLabel;
        end
        curratedflag = 1;
        if isfield(clusinfo, 'cluster_id')
            clusidtmp = clusinfo.cluster_id;
            cluster_id = [cluster_id, clusinfo.cluster_id];
        else
            keyboard
            disp('Someone thought it was nice to change the name again...')
        end


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
                depthtmp(clusid) = channelpostmp(minidx, 2);
                xtmp(clusid) = channelpostmp(minidx, 1);
            catch
                channeltmp(clusid) = channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmp(minidx-1, 2);
                xtmp(clusid) = channelpostmp(minidx-1, 1);
            end
            sp.spikeDepths(ismember(sp.clu, clusidtmp(clusid))) = depthtmp(clusid);

        end
        depth = [depth, depthtmp];
        channel = [channel, channeltmp];
        Good_IDtmp = ismember(cellstr(clusinfo.group), 'good');
    end

    if any(emptyclus) % TEmporarily
        paramBC.reextractRaw = 1;
    end

    % channelpostmp = channelpostmpconv;
    ypostmp = channelpostmp(:, 2);
    xpostmp = channelpostmp(:, 1);
    xdiffs = unique(abs(diff(unique(xpostmp))));
    xdiffs(xdiffs<50) = []; % Assume differences smaller than 50 micron means they're on the same shank?
    xposopt = (floor(xpostmp./min(xdiffs)+1)); % Sort of hacky
    Shanktmp = floor(xpostmp(channeltmp+1)./(min(xdiffs)+1));
    if isempty(xposopt)
        xposopt = ones(length(channeltmp),1);
        Shanktmp = zeros(length(channeltmp),1);
    end
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
        Donotinclude = 0;
        theseuniqueTemplates = [];
        unitTypeAcrossRec = [];

        %% Quality metrics - Bombcell (https://github.com/Julie-Fabre/bombcell)
        for id = 1:length(rawD)
            statusCopy = 0;
            ephysap_tmp = [];

            %ephysap_path = fullfile(rawD(id).folder, rawD(id).name);

            if length(rawD) > 1 % DO NOT DELETE!
                savePath = fullfile(myClusFile(1).folder, num2str(id));
                kilosortPath = myClusFile(1).folder;
                idx = sp.SessionID == id;
            else
                savePath = fullfile(KiloSortPaths{subsesid});
                kilosortPath = savePath;
                idx = sp.SessionID == 1;
            end

            qMetricsExist = ~isempty(dir(fullfile(savePath, '**', 'templates._bc_qMetrics.parquet'))); % ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) not used anymore?

            if ~qMetricsExist & ~Params.ExtractNewDataNow
                disp('No new extractions now... continue')
                Donotinclude = 1;
                IncludeThese(subsesid) = false;

                break
            end

            InspectionFlag = 0;
            if isempty(dir(fullfile(rawD(id).folder, strrep(rawD(id).name, '.cbin', '_sync.dat')))) 
                disp('Extracting sync file...')
                % detect whether data is compressed, decompress locally if necessary
                if ~exist(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'))) && ~contains(rawD(1).name,'.dat') && ~contains(rawD(1).name,'.raw')
                    disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                    decompDataFile = bc.dcomp.extractCbinData(fullfile(rawD(id).folder, rawD(id).name), ...
                        [], [], 0, fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')));
                    statusCopy = copyfile(strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'), strrep(fullfile(Params.tmpdatafolder, rawD(id).name), 'cbin', 'meta')); %QQ doesn't work on linux
                end
                Params.DecompressionFlag = 1;
                [~,~,currext] = fileparts(rawD(id).name);

                if ~exist('statusCopy','var') ||statusCopy == 0 %could not copy meta file - use original meta file
                    [Imecmeta] = ReadMeta2(fullfile(rawD(id).folder, strrep(rawD(id).name, currext, '.meta')), 'ap');
                else
                    [Imecmeta] = ReadMeta2(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, currext, '.meta')), 'ap');
                end
                nchan = strsplit(Imecmeta.acqApLfSy, ',');
                nChansInFile = str2num(nchan{1}) + str2num(nchan{3});

                syncDatImec = extractSyncChannel(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, currext, '.bin')), nChansInFile, nChansInFile); %Last channel is sync (function of spikes toolbox)
                statusCopy = copyfile(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, currext, '_sync.dat')), fullfile(rawD(id).folder, strrep(rawD(id).name, '.cbin', '_sync.dat'))); %QQ doesn't work on linux

            end
            if ~qMetricsExist || Params.RedoQM
                % Check if all waveforms are extracted
                if ~isempty(dir(fullfile(savePath, 'templates._bc_rawWaveforms_kilosort_format.npy')))
                    rawWaveformsFull = readNPY(fullfile(savePath, 'templates._bc_rawWaveforms_kilosort_format.npy'));

                    % Find empty rows (where all elements are zero)
                    row_sums = squeeze(sum(sum(abs(rawWaveformsFull), 2), 3));
                    empty_row_indices = find(isnan(row_sums));
                else
                    empty_row_indices = [];
                end
             
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(rawD(id).name, 'cbin')) && Params.DecompressLocal && (isempty(dir(fullfile(savePath, '**', 'RawWaveforms'))) || (Params.RedoQM && paramBC.reextractRaw) || any(empty_row_indices))  % if raw waveforms have not been extract, decompress data for extraction
                    % detect whether data is compressed, decompress locally if necessary
                    if ~exist(fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'))) && ~contains(rawD(1).name,'.dat') && ~contains(rawD(1).name,'.raw')
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                        decompDataFile = bc.dcomp.extractCbinData(fullfile(rawD(id).folder, rawD(id).name), ...
                            [], [], 0, fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin')));
                        statusCopy = copyfile(strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'), strrep(fullfile(Params.tmpdatafolder, rawD(id).name), 'cbin', 'meta')); %QQ doesn't work on linux
                    end

                    Params.DecompressionFlag = 1;
                    if Params.InspectQualityMetrics
                        InspectionFlag = 1;
                    end
                end

                paramBC.rawFile = fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'));
                paramBC.ephysMetaFile = (strrep(fullfile(rawD(id).folder, rawD(id).name), 'cbin', 'meta'));

                %             idx = ismember(sp.spikeTemplates,clusidtmp(Good_IDtmp)); %Only include good units
                %careful; spikeSites zero indexed
                if isempty(sp.pcFeat)
                    spfeat = [];
                else
                    spfeat = sp.pcFeat(idx, :, :);
                end
                [qMetric, unitType] = bc.qm.runAllQualityMetrics(paramBC, sp.st(idx)*sp.sample_rate, sp.clu(idx)+1, ...
                    templateWaveforms, sp.tempScalingAmps(idx), spfeat, sp.pcFeatInd+1, channelpostmp, savePath); % Be careful, bombcell needs 1-indexed!

            else
                paramBC.rawFile = fullfile(Params.tmpdatafolder, strrep(rawD(id).name, 'cbin', 'bin'));
                if exist('ephysap_tmp', 'var')
                    paramBC.tmpFolder = [ephysap_tmp, '/..'];
                end
                d = dir(fullfile(savePath, '**', 'templates._bc_qMetrics.parquet'));
                qMetricsPath = d.folder;
                [~, qMetric, fractionRPVs_allTauR] = bc.load.loadSavedMetrics(qMetricsPath);
                unitType = bc.qm.getQualityUnitType(paramBC, qMetric);
                %                 unitType(:) = 1; ???


                % Commmented by CB for now
                spike_templates_0idx = readNPY([kilosortPath filesep 'spike_clusters.npy']); % changed back 20230920 JF

                spikeTemplates = spike_templates_0idx + 1;
                uniqueTemplates = unique(spikeTemplates);
                % need to load forGUI.tempWv??
                try
                    tmpfile = dir(fullfile(savePath,'**','templates.qualityMetricDetailsforGUI.mat'));
                    tmpGUI = load(fullfile(tmpfile.folder,tmpfile.name));
                    bc.qm.plotGlobalQualityMetric(qMetric, paramBC, unitType, uniqueTemplates, tmpGUI.forGUI.tempWv);
                catch ME
                    disp(ME)
                end

                %                 load(fullfile(savePath, 'qMetric.mat'))
            end
            theseuniqueTemplates{id} = unique(sp.clu);
            qMetricclusterID = qMetric.clusterID;

            unitTypeAcrossRec{id} = unitType;

            if InspectionFlag % Doesn't currently work: Julie will update bombcell
                bc.load.loadMetricsForGUI

                %% Inspect the results?
                % GUI guide:
                % left/right arrow: toggle between units
                % g : go to next good unit
                % m : go to next multi-unit
                % n : go to next noise unit
                % up/down arrow: toggle between time chunks in the raw data
                % u: brings up a input dialog to enter the unit you want to go to
                unitQualityGuiHandle = bc.viz.unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
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

        if Donotinclude
            continue
        end
        AllUniqueTemplates = cat(1, AllUniqueTemplates(:), cat(1, theseuniqueTemplates{:}));

    else
        AllUniqueTemplates = cat(1, AllUniqueTemplates, clusinfo.cluster_id);
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
        ProbeAll = cat(1, ProbeAll(:), repmat(probeSN, 1, length(cat(2, recsesAlltmp{:})))'); % Replace to serial number
    else
        recsesAll = cat(1, recsesAll(:), repmat(addthis+1, 1, length(Good_IDtmp))');
        ProbeAll = cat(1, ProbeAll(:), repmat(probeSN, 1, length(Good_IDtmp))'); % Replace to serial numer
    end
    sp.RecSes = sp.SessionID + countid - 1; %Keep track of recording session, as cluster IDs are not unique across sessions

    %% Save out sp and clusinfo for this session in the correct folder
    ShankOpt = unique(Shank);
    ShankID = nan(size(Shank));
    for shankid = 1:length(ShankOpt)
        ShankID(Shank == ShankOpt(shankid)) = shankid;
    end
    if numel(Shank)~=numel(AllUniqueTemplates)
        keyboard
    end
    clusinfo.Shank = Shank;
    clusinfo.ProbeID = ProbeAll;
    clusinfo.ShankID = ShankID;
    clusinfo.RecSesID = recsesAll;
    clusinfo.ch = channel;
    clusinfo.depth = depth;
    clusinfo.cluster_id = AllUniqueTemplates;
    clusinfo.group = Label;
    if length(Good_ID) > length(recsesAll)
        NonEmptyIdx = true(1, length(Good_ID));
        NonEmptyIdx(emptyclus) = false;
        clusinfo.Good_ID = Good_ID(NonEmptyIdx);
    else
        clusinfo.Good_ID = Good_ID;
    end

    if exist('Depth2Area', 'var') % Add area information
        Idx = cell2mat(arrayfun(@(X) find(Depth2Area.Cluster_ID-1 == X), clusinfo.cluster_id, 'Uni', 0));
        clusinfo.Area = Depth2Area.Area(Idx);
        clusinfo.Coordinates = Depth2Area.Coordinates(Idx);
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

    SessionParams = Params; % But only store for this specific session - otherwise confusion!
    SessionParams.KSDir = KiloSortPaths(subsesid);
    SessionParams.RawDataPaths = RawDataPaths(subsesid);
    SessionParams.AllChannelPos = {channelpostmp};
    SessionParams.AllProbeSN = {probeSN};

    save(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'), 'clusinfo', 'SessionParams', '-v7.3')
    if Params.saveSp %QQ not using savePaths
        try
            save(fullfile(KiloSortPaths{subsesid}, 'PreparedData.mat'), 'sp', '-append')
        catch ME
            disp(ME)
        end
        % Store all figures
        if Params.RunQualityMetrics
            figHandles = findall(0,'Type','figure');
            for figid = 1:numel(figHandles)
                try
                saveas(figHandles(figid),fullfile(KiloSortPaths{subsesid},['Fig' num2str(figid) '.fig']))
                saveas(figHandles(figid),fullfile(KiloSortPaths{subsesid},['Fig' num2str(figid) '.bmp']))
                catch
                end
            end
        end


    end
    close all
    countid = countid + 1;
end

Params.KSDir = KiloSortPaths(IncludeThese);
Params.AllChannelPos = AllChannelPos(IncludeThese);
Params.AllProbeSN = AllProbeSN(IncludeThese);
Params.RawDataPaths = RawDataPaths(IncludeThese);

%% Remove temporary files
if isstruct(RawDataPaths)
    if any(cellfun(@(X) strcmp(X,Params.tmpdatafolder), {RawDataPaths(:).folder}))
        Params.CleanUpTemporary = 1;
    end
end

CleanUpCheckFlag = 1; % Put to 1 is own responsibility! Make sure not to delete stuff from the server directly!
if Params.DecompressLocal && Params.CleanUpTemporary
    try
        if isnan(CleanUpCheckFlag) && any(cellfun(@(X) exist(fullfile(Params.tmpdatafolder, strrep(X.name, 'cbin', 'bin'))),RawDataPaths(find(~cellfun(@isempty,RawDataPaths)))))

            answer = questdlg(['Automatically remove data from ', Params.tmpdatafolder, '?'], ...
                'REMOVING -- CHECK!!!', ...
                'YES', 'NO', 'YES');
            if strcmpi(answer, 'YES')
                CleanUpCheckFlag = 1;
            else
                CleanUpCheckFlag = 0;
            end
        elseif isnan(CleanUpCheckFlag)
            CleanUpCheckFlag = 0;
        end
    catch
        CleanUpCheckFlag = 0;
    end

    if CleanUpCheckFlag
        clear memMapData
        clear ap_data
        try
            for id = 1:length(RawDataPaths)
                if ~isempty(RawDataPaths{id})
                    delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths{id}.name, 'cbin', 'bin')))
                    delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths{id}.name, 'cbin', 'meta')))
                    delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths{id}.name, '.ap.cbin', '_kilosortChanMap.mat')))
                    delete(fullfile(Params.tmpdatafolder, strrep(RawDataPaths{id}.name, '.cbin', '_sync.dat')))
                    delete(fullfile(Params.tmpdatafolder, RawDataPaths{id}.name))
                end
            end
        catch ME
            keyboard
        end
    end
end

return
