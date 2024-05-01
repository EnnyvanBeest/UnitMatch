function spAll = getSpikesFromPrepData(AllKSDir)

    % Load SP
    nKSFiles = length(AllKSDir);
    spAll = cell(1, nKSFiles);
    for did = 1:nKSFiles
        % Load prepared data
        if exist(fullfile(AllKSDir{did}, 'PreparedData.mat'))
            try
                load(fullfile(AllKSDir{did}, 'PreparedData.mat'), 'sp', 'SessionParams');
            catch % Old way of saving?
                load(fullfile(AllKSDir{did}, 'PreparedData.mat'), 'sp', 'Params');
                SessionParams = Params;
            end
            if ~exist('SessionParams')
                load(fullfile(AllKSDir{did}, 'PreparedData.mat'), 'Params');
                SessionParams = Params;
            end

        else
            SessionParams.RunPyKSChronicStitched = 0;
            warning('No PreparedData.mat found... loading in directly from KS.. have not checked for empty clusters... Consider using ExtractKilosortData.m')

            %% Load Spike Data
            sp = loadKSdir(AllKSDir{did}); % Load Spikes with PCs
            [sp.spikeAmps, sp.spikeDepths, sp.templateDepths, sp.templateXpos, sp.tempAmps, sp.tempsUnW, sp.templateDuration, sp.waveforms] = ...
                templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.xcoords, sp.spikeTemplates, sp.tempScalingAmps); %from the spikes toolbox
        end

        % Only keep parameters used
        spAll{did}.st = sp.st;
        spAll{did}.spikeTemplates = sp.spikeTemplates;
        spAll{did}.spikeAmps = sp.spikeAmps;
        spAll{did}.spikeDepths = sp.spikeDepths;
        % Replace recsesid with subsesid
        if SessionParams.RunPyKSChronicStitched
            spAll{did}.RecSes = sp.RecSes;
        else
            spAll{did}.RecSes = repmat(did, size(spAll{did}.st));
        end
    end
    
    % Add all spikedata in one spikes struct - can be used for further analysis
    spAll = [spAll{:}];
    spnew = struct;
    fields = fieldnames(spAll(1));
    for fieldid = 1:length(fields)
        try
            spnew.(fields{fieldid}) = cat(1,spAll(:).(fields{fieldid}));
        catch ME
            if strcmp(ME.message, 'Out of memory.')
                spnew.(fields{fieldid}) = cat(1,spAll(1).(fields{fieldid}));
                for tmpid = 2:length(spAll)
                    spnew.(fields{fieldid}) = cat(1,spnew.(fields{fieldid}),spAll(tmpid).(fields{fieldid}));
                end
            else
                spnew.(fields{fieldid}) = cat(1,spAll(:).(fields{fieldid}));
                eval(['spnew.', fields{fieldid}, '= cat(2,sp(:).', fields{fieldid}, ');'])
            end
        end
    end
    spAll = spnew;

end