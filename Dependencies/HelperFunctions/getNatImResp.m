function [spikeData,proc] = getNatImResp(spikesAll,expFolders,binFileRef,proc)
    %%% Load natural images responses
    % spikesAll.times: spike times
    % spikesAll.clusters: cluster ID for each spike
    % expFolders: folder for each experiment to load for that spike train
    % binFileRef: associated ephys file

    % Get processing parameters
    if ~exist('proc','var')
        proc.window = [-0.3 0.5 ... % around onset
            0.0 0.5]; % around offset
        proc.binSize = 0.002; % in ms
        proc.smoothSize = 5; % PSTH smoothing filter
        gw = gausswin(proc.smoothSize,3);
        proc.smWin = gw./sum(gw);
    end
    nBins = int64((proc.window(2) - proc.window(1) + proc.window(4) - proc.window(3))/proc.binSize);

    clusters.IDs = spikesAll.clusterIDs;

    % Loop through exp
    for ee = 1:numel(expFolders)
        expFolder = expFolders{ee};
        blockFile = dir(fullfile(expFolder,'*Block*'));
        load(fullfile(blockFile.folder,blockFile.name));
        
        if contains(block.rigName,'zelda')
            % Zelda exp -- already aligned.

            % subject = regexp(expFolder,'[A-Z]{2,}\d{3}','match'); subject = subject{1};
            % expDate = regexp(expFolder,'\d{4}-\d{2}-\d{2}','match'); expDate = expDate{1};
            % [~,expNum] = fileparts(expFolder);

            % Get events
            % expInfo = csv.queryExp(subject=subject,expDate=expDate,expNum=expNum);
            % events = csv.loadData(expInfo,dataType = 'eventsFull');
            evPQTPath = dir(fullfile(expFolder,'ONE_preproc','events',['_av_trials.table.*.pqt' ]));
            if ~isempty(evPQTPath)
                events = table2struct(parquetread(fullfile(evPQTPath.folder,evPQTPath.name)),"ToScalar",1);
            else
                error('Missing events. Check?')
            end

            imageOnsetTimes = events.imageOnsetTimes;
            imageOffsetTimes = events.imageOffsetTimes;
            imageIDs = events.imageIDs;

            % align ephys
            alignmentFile = dir(fullfile(expFolder,'*alignment.mat'));
            alignment = load(fullfile(alignmentFile.folder,alignmentFile.name));
            probeNum = [];
            for pp = 1:numel(alignment)
                binFileAlign = dir(fullfile(alignment.ephys(pp).ephysPath,'*ap.*bin'));
                binFileRefDir = dir(binFileRef);
                if strcmp(binFileAlign.name, binFileRefDir.name)
                    probeNum = pp;
                end
            end

            if isempty(probeNum)
                error('Problem with ephys ref, cannot find its alignment')
            end

            % Get spikes time -- in timeline time
            spikes.times = interp1(alignment.ephys(probeNum).originTimes, alignment.ephys(probeNum).timelineTimes, spikesAll.times, 'linear', nan);
            nanIdx = isnan(spikes.times);

            if any(nanIdx)
                refOffsets = alignment.ephys(probeNum).timelineTimes(:)-alignment.ephys(probeNum).originTimes(:);
                offSetsPerPoint = interp1(alignment.ephys(probeNum).originTimes, refOffsets, spikesAll.times, 'nearest', 'extrap');

                spikes.times(nanIdx) = spikesAll.times(nanIdx)+offSetsPerPoint(nanIdx);
            end
            
            expLength = block.duration;
            spk2keep = (spikes.times>0) & (spikes.times<expLength);
            spikes.times = spikes.times(spk2keep);
            spikes.clusters = spikesAll.clusters(spk2keep);
        
        elseif contains(block.rigName,'zgood')
            % Kilotrode
            
            spikes = spikesAll;
            alignmentFolder = fullfile(fileparts(expFolder),'alignments');
            [~,expNum] = fileparts(expFolder);

            % Get block alignment
            alignmentFileBlock = dir(fullfile(alignmentFolder, sprintf('correct_block_%s_to_timeline_*.npy',expNum)));
            % Get timeline ref
            timeRef = extractAfter(alignmentFileBlock.name,"timeline_");
            timeRef = timeRef(1:end-4);

            % Get ephys alignment
            alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_*.npy',timeRef))); % should be the only one, has already been aligned
            if ~(numel(alignmentFileEphys) == 1)
                error('Couldn''t find unique alignement file for ephys. Check?')
            end

            % Get event times -- in ephys time
            %%% SHOULD BE TAKEN FROM TIMELINE
            al = readNPY(fullfile(alignmentFileBlock.folder,alignmentFileBlock.name));
            bTLtoMaster = readNPY(fullfile(alignmentFileEphys.folder,alignmentFileEphys.name));
            eventTimes = block.stimWindowUpdateTimes;
            eventTimes = applyCorrection(eventTimes, al);
            eventTimes = applyCorrection(eventTimes, bTLtoMaster);
            imageOnsetTimes = eventTimes(1:2:end);
            imageOffsetTimes = eventTimes(2:2:end);
            %%% HACK 
            imageOffsetTimes = [imageOffsetTimes; imageOnsetTimes(end)+median(imageOffsetTimes-imageOnsetTimes(1:end-1))];
            imageIDs = block.events.numValues;
        end

        nClusters = numel(clusters.IDs);
        nTrials = numel(imageOnsetTimes);
        baSmtmp = zeros(nTrials, nBins, nClusters);

        % get all onset-centered psths
        for c = 1:nClusters
            temp = clusters.IDs(c);
            st = spikes.times(spikes.clusters == temp);

            % get psth
            if diff(proc.window(1:2)) > 0
                [~, ~, ~, ~, ~, baOn] = psthAndBA(st, imageOnsetTimes, proc.window(1:2), proc.binSize);
            else
                baOn = [];
            end
            if diff(proc.window(3:4)) > 0
                [~, ~, ~, ~, ~, baOff] = psthAndBA(st, imageOffsetTimes, proc.window(3:4), proc.binSize);
            else
                baOff = [];
            end

            % smooth ba
            ba = cat(2,baOn,baOff);
            baSmtmp(:,:,c) = conv2(proc.smWin,1,ba', 'same')'./proc.binSize;
        end

        % not optimal here?
        trials = imageIDs(1:nTrials);
        trialid = unique(trials);
        baSm{ee} = nan(numel(trialid), nBins, nClusters, ceil(nTrials/numel(trialid)));
        for tt = 1:numel(trialid)
            idxrep = trials == trialid(tt);
            baSm{ee}(tt,:,:,1:sum(idxrep)) = permute(baSmtmp(idxrep,:,:),[2 3 1]);
        end
    end
    spikeData = cat(4,baSm{:});
end
   