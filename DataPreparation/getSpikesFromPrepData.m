function spAll = getSpikesFromPrepData(AllKSDir)

    % Load SP
    nKSFiles = length(AllKSDir);
    spAll = cell(1, nKSFiles);
    for did = 1:nKSFiles
        % Load prepared data
        load(fullfile(AllKSDir{did}, 'PreparedData.mat'), 'sp', 'Params');
    
        % Only keep parameters used
        spAll{did}.st = sp.st;
        spAll{did}.spikeTemplates = sp.spikeTemplates;
        spAll{did}.spikeAmps = sp.spikeAmps;
        spAll{did}.spikeDepths = sp.spikeDepths;
        % Replace recsesid with subsesid
        if Params.RunPyKSChronicStitched
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