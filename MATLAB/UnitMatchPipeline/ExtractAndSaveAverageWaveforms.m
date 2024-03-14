function Path4UnitNPY = ExtractAndSaveAverageWaveforms(clusinfo,param)
%% Called by UnitMatch, but can also be used on its own to save two averaged waveforms per unit per session

%% Read in from param
RedoExtraction = param.RedoExtraction; % Raw waveform and parameter extraction
AllDecompPaths = param.AllDecompPaths;
sampleamount = param.sampleamount; %500; % Nr. waveforms to include
spikeWidth = param.spikeWidth; % in sample space (time) - number of samples
halfWidth = floor(spikeWidth/2);
%% Extract all cluster info
AllClusterIDs = clusinfo.cluster_id;
AllRecSes = clusinfo.RecSesID;
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')
end
GoodRecSesID = clusinfo.RecSesID(Good_Idx);

% Define day stucture
nclus = length(Good_Idx);
UseMemMap = 1;

%% Actual extraction
dataTypeNBytes = numel(typecast(cast(0, 'uint16'), 'uint8')); % Define datatype

% Initialize
Path4UnitNPY = cell(1,nclus);

timercounter = tic;
fprintf(1,'Extracting raw waveforms. Progress: %3d%%',0)
Currentlyloaded = 0;
for uid = 1:nclus
    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    if length(param.KSDir)>1 || length(param.RawDataPaths) == 1
        tmppath = dir(fullfile(param.KSDir{GoodRecSesID(uid)},'**','RawWaveforms*'));
    else %Stitched KS
        tmppath = dir(fullfile(param.KSDir{1},'**','RawWaveforms*'));
    end
    if length(tmppath)>1
        % Probably stitched:
        tmppath = tmppath(GoodRecSesID(uid));
    end
    if isempty(tmppath)
        mkdir(fullfile(param.KSDir{GoodRecSesID(uid)},'RawWaveforms'));
        tmppath = dir(fullfile(param.KSDir{GoodRecSesID(uid)},'**','RawWaveforms*'));
    end
    Path4UnitNPY{uid} = fullfile(tmppath.folder,tmppath.name,['Unit' num2str(AllClusterIDs(Good_Idx(uid))) '_RawSpikes.npy']); %0-indexed

    if exist(Path4UnitNPY{uid}) && ~RedoExtraction
        continue
    else       
        if ~(GoodRecSesID(uid) == Currentlyloaded) % Only load new memmap if not already loaded

            % Map the data
            if ~UseMemMap & exist('fid')
                fclose(fid);
            end
            clear memMapData
            if ~isstruct(AllDecompPaths{GoodRecSesID(uid)})
                spikeFile = dir(AllDecompPaths{GoodRecSesID(uid)});
            else
                spikeFile = AllDecompPaths{GoodRecSesID(uid)};
                AllDecompPaths{GoodRecSesID(uid)} = fullfile(AllDecompPaths{GoodRecSesID(uid)}.folder, AllDecompPaths{GoodRecSesID(uid)}.name);
            end
            % Check memory
            [userview,systemview] = memory;
            if spikeFile.bytes>userview.MaxPossibleArrayBytes
                disp('Cannot load raw data to memory, serial loading instead')
                UseMemMap = 0;
            end
            if UseMemMap
                try %hacky way of figuring out if sync channel present or not
                    n_samples = spikeFile.bytes / (param.nChannels * dataTypeNBytes);
                    nChannels = param.nChannels - param.nSyncChans; % Last channel is sync, ignore for now
                    ap_data = memmapfile(AllDecompPaths{GoodRecSesID(uid)}, 'Format', {'int16', [param.nChannels, n_samples], 'data'});
                catch
                    ap_data = memmapfile(AllDecompPaths{GoodRecSesID(uid)}, 'Format', {'int16', [nChannels, n_samples], 'data'});
                end
                memMapData = ap_data.Data.data;
            else
                n_samples = spikeFile.bytes / (param.nSavedChans * dataTypeNBytes);
                nChannels = param.nSavedChans - param.nSyncChans; % Last channel is sync, ignore for now
                fid = fopen(AllDecompPaths{GoodRecSesID(uid)},'r');
            end
            Currentlyloaded = GoodRecSesID(uid);
        end

        %load sp
        if length(param.KSDir)>1
            tmp = matfile(fullfile(param.KSDir{GoodRecSesID(uid)},'PreparedData.mat'));
        else %Stitched Kilosort
            tmp = matfile(fullfile(param.KSDir{1},'PreparedData.mat'));
        end
        sp = tmp.sp;
        tmpclusinfo = tmp.clusinfo; %Load original clusinfo

        % Spike samples
        if length(param.KSDir)>1
            idx1=(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid)) & sp.SessionID == tmpclusinfo.RecSesID(Good_Idx(uid)- find(AllRecSes==GoodRecSesID(uid),1,'first')+1)).*round(sp.sample_rate));  % Spike times in samples;
        else
            idx1=(sp.st(sp.spikeTemplates == AllClusterIDs(Good_Idx(uid)) & sp.SessionID == tmpclusinfo.RecSesID(Good_Idx(uid))).*round(sp.sample_rate));  % Spike times in samples;
        end
        %Extract raw waveforms on the fly - % Unit uid
        if sampleamount<length(idx1)
            spksubi = round(linspace(1, numel(idx1), sampleamount))';
            spikeIndicestmp = sort(idx1(spksubi));
        else
            spikeIndicestmp = sort(idx1);
        end
        spikeMap = nan(spikeWidth,nChannels,length(spikeIndicestmp));
        for iSpike = 1:length(spikeIndicestmp)
            thisSpikeIdx = spikeIndicestmp(iSpike);

            if UseMemMap
                if thisSpikeIdx > halfWidth && (thisSpikeIdx + halfWidth) < size(memMapData,2) % check that it's not out of bounds
                    tmp = smoothdata(double(memMapData(1:param.nSavedChans,thisSpikeIdx-halfWidth:thisSpikeIdx+halfWidth)),2,'gaussian',5);
                    tmp = (tmp - mean(tmp(:,1:20),2))';
                    tmp(:,end+1:nChannels) = nan(size(tmp,1),nChannels-size(tmp,2));
                    % Subtract first 10 samples to level spikes
                    spikeMap(:,:,iSpike) = tmp(1:spikeWidth,1:nChannels);
                end
            else
                if ((thisSpikeIdx - halfWidth) * param.nSavedChans) * dataTypeNBytes > halfWidth &&...
                        (thisSpikeIdx + halfWidth) * param.nSavedChans * dataTypeNBytes < spikeFile.bytes %
                    bytei = ((thisSpikeIdx - halfWidth) * param.nSavedChans) * dataTypeNBytes;
                    fseek(fid, bytei, 'bof');
                    data0 = fread(fid, param.nSavedChans*spikeWidth, 'int16=>int16'); % read individual waveform from binary file
                    frewind(fid);
                    data = reshape(data0, param.nSavedChans, []);
                    tmp = smoothdata(double(data(1:nChannels,:)),2,'gaussian',5);
                    tmp = (tmp - mean(tmp(:,1:20),2))';
                    tmp(:,end+1:nChannels) = nan(size(tmp,1),nChannels-size(tmp,2));
                    spikeMap(:,:,iSpike) = tmp(1:spikeWidth,1:nChannels);
                end
            end
        end
        %Actual number of wavefroms
        nwavs = sum(sum(~isnan(nanmean(spikeMap,2)),1) == spikeWidth); % Actual number of waves
        for cv = 1:2
            if cv==1
                wavidx = floor(1:nwavs/2);
            else
                wavidx = floor(nwavs/2+1:nwavs);
            end
            spikeMapAvg(:,:,cv) = nanmedian(spikeMap(:,:,wavidx),3);
        end
        spikeMap = spikeMapAvg;
        clear spikeMapAvg
        writeNPY(spikeMap, Path4UnitNPY{uid})

    end
end
if ~UseMemMap & exist('fid')
    fclose(fid);
else
    clear memMapData
end

fprintf('\n')
disp(['Extracting raw waveforms took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])
