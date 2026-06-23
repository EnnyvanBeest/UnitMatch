function [imageOnsetTimes, imageOffsetTimes, imageIDs] = extractStimTimings(expFolder, binFileRef)
    blockFile = dir(fullfile(expFolder,'*Block*'));
    load(fullfile(blockFile.folder,blockFile.name),'block');

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

        % Get image onsets / offsets in spike time
        imageOnsetTimes = interp1(alignment.ephys(probeNum).timelineTimes, alignment.ephys(probeNum).originTimes, imageOnsetTimes, 'linear', nan);
        imageOffsetTimes = interp1(alignment.ephys(probeNum).timelineTimes, alignment.ephys(probeNum).originTimes, imageOffsetTimes, 'linear', nan);


    elseif contains(block.rigName,'zgood')
        % Kilotrode
        
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

    nanIdx = isnan(imageOnsetTimes) | isnan(imageOffsetTimes);
    imageOnsetTimes(nanIdx) = [];
    imageOffsetTimes(nanIdx) = [];

    nTrials = numel(imageOnsetTimes);

    % not optimal here?
    imageIDs = imageIDs(1:nTrials);
end