function exp2keep = getNatImExpRef(binFile)
    %%% Find the experiments associated with a binFile

    if ~isstruct(binFile)
        binFile = dir(binFile);
    end

    % Get info -- hardcoded but shouldn't change...?
    info = split(binFile.folder,'\');
    server = ['\\' info{3} '\' info{4}];
    subject = info{5};
    expDate = info{6};

    exp2keep = [];
    expFolders = dir(fullfile(server,subject,expDate,'*'));
    expFolders = expFolders(cell2mat(cellfun(@(x) ~isempty(regexp(x,'\d*','match')), {expFolders.name}, 'uni', 0)));

    % Get associated bin file (non-stitched)
    for ee = 1:numel(expFolders)
        expFolder = fullfile(expFolders(ee).folder,expFolders(ee).name);
        blockFile = dir(fullfile(expFolder,'*Block*'));
        if isempty(blockFile)
            continue
        end
        load(fullfile(blockFile.folder,blockFile.name));

        if contains(block.expDef,'imageWorld')
            if contains(block.rigName,'zelda')

                % On a zelda rig
                alignmentFile = dir(fullfile(expFolder,'*alignment.mat'));
                alignment = load(fullfile(alignmentFile.folder,alignmentFile.name));

                % Check that that recording was aligned
                for probeNum = 1:numel(alignment)
                    if strcmp(alignment.ephys(probeNum).ephysPath,binFile.folder) %%% SHOULD BE ENOUGH TO IDENTIFY RECORDING?
                        exp2keep = cat(1,exp2keep,{expFolder});
                    end
                end

            elseif contains(block.rigName,'zgood')

                % On kilotrode
                alignmentFolder = fullfile(fileparts(expFolder),'alignments');
                if exist(alignmentFolder,'dir')

                    % Check if block is associated with a timeline
                    [~,expNum] = fileparts(expFolder);
                    alignmentFileBlock = dir(fullfile(alignmentFolder, sprintf('correct_block_%s_to_timeline_*.npy',expNum)));
                    if ~isempty(alignmentFileBlock)
                        % Get timeline ref
                        timeRef = extractAfter(alignmentFileBlock.name,"timeline_");
                        timeRef = timeRef(1:end-4);
                    else
                        warning('Block wasn''t aligned to timeline.')
                        continue
                    end

                    % Check if ephys is associated with the same
                    % timeline -- could do it the other way around with
                    % ephys before but never sure will be in kilotrode
                    [~,ephysFolder] = fileparts(binFile.folder);
                    alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_%s.npy',timeRef,ephysFolder)));
                    if ~isempty(alignmentFileEphys)
                        exp2keep = cat(1,exp2keep,{expFolder});
                    else
                        warning('Ephys wasn''t aligned to timeline.')
                        continue
                    end
                end
            end
        end
    end

    exp2keep = unique(exp2keep);

    if isempty(exp2keep)
        warning('Couldn''t find any associated Natural Images experiment.')
    end
end