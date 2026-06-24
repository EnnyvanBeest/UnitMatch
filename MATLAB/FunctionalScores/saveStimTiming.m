function saveStimTiming(UMfolder, recompute)

    if nargin < 2 || isempty(recompute)
        recompute = 0;
    end

    d = dir(fullfile(UMfolder,'**','UnitMatch.mat'));
    for dd = 1:numel(d)
        SaveDir = d(dd).folder;
        load(fullfile(SaveDir, 'UnitMatch.mat'), 'UMparam');

        %%% TO REMOVE
        if ~isfield(UMparam,'AllRawPaths') % For now no to disturb Célian
            UMparam.AllRawPaths = UMparam.RawDataPaths;
        end
        
        for ss = 1:numel(UMparam.AllRawPaths)
        
            % Get the original binFile (also for stitched?)
            Flag = 1;
            if ~isempty(UMparam.AllRawPaths{ss}) % When no raw data is available
                if iscell(UMparam.AllRawPaths{ss})
                    binFileRef = fullfile(UMparam.AllRawPaths{s});
                elseif isstr(UMparam.AllRawPaths{ss})
                    binFileRef = UMparam.AllRawPaths{ss};
                else
                    binFileRef = fullfile(UMparam.AllRawPaths{ss}.folder,UMparam.AllRawPaths{ss}.name);
                end
            else
                Flag = 0;
            end
        
            if ~exist(fullfile(fileparts(binFileRef), 'trial.onsetTimes.npy'), 'file') || recompute
        
                % Find the associated experiments
                if Flag
                    exp2keep = getNatImExpRef(binFileRef);
        
                    if ~isempty(exp2keep)
        
                        % Loop through exp
                        imageOnsetTimesAll = [];
                        imageOffsetTimesAll = [];
                        imageIDsAll = [];
                        for ee = 1:numel(exp2keep)
                            expFolder = exp2keep{ee};
                            [imageOnsetTimes, imageOffsetTimes, imageIDs] = extractStimTimings(expFolder, binFileRef);
        
                            imageOnsetTimesAll = cat(1, imageOnsetTimesAll, imageOnsetTimes(:));
                            imageOffsetTimesAll = cat(1, imageOffsetTimesAll, imageOffsetTimes(:));
                            imageIDsAll = cat(1, imageIDsAll, imageIDs(:));
                        end
        
                        expPath = fileparts(binFileRef);
                        writeNPY(imageOnsetTimesAll, fullfile(expPath, 'trial.onsetTimes.npy'));
                        writeNPY(imageOffsetTimesAll, fullfile(expPath, 'trial.offsetTimes.npy'));
                        writeNPY(imageIDsAll, fullfile(expPath, 'trial.imageIDs.npy'));
                    end
                end
            end
        end
    end
end

