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
        if numel(blockFile)>1
            blockFile = blockFile(contains({blockFile.name},expDate));
        end
        if isempty(blockFile)
            continue
        end
        load(fullfile(blockFile.folder,blockFile.name));
        
        if numel(block.stimWindowUpdateTimes) < 223
            warning('Incomplete block. Skip.\n')
            continue
        end

        if contains(block.expDef,'imageWorld')
            if contains(block.rigName,'zelda')

                % On a zelda rig
                alignmentFile = dir(fullfile(expFolder,'*alignment.mat'));
                if isempty(alignmentFile)
                    continue
                end
                alignment = load(fullfile(alignmentFile.folder,alignmentFile.name));

                % Check that that recording was aligned
                for probeNum = 1:numel(alignment)
                    if ~isa(alignment.ephys(probeNum),'double')
                        dAlign = dir(fullfile(alignment.ephys(probeNum).ephysPath,'*cbin'));
                        if ~isempty(dAlign)
                            if strcmp(dAlign.name,binFile.name) %%% SHOULD BE ENOUGH TO IDENTIFY RECORDING?
                                exp2keep = cat(1,exp2keep,{expFolder});
                            end
                        else
                            % can happen when there was an alignment error
                            fprintf('Alignment issue in the pinkrigs. Skip.\n',ee)
                        end
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
                        fprintf('Block %d wasn''t aligned to timeline.\n',ee)
                        continue
                    end

                    % Check if ephys is associated with the same
                    % timeline -- could do it the other way around with
                    % ephys before but never sure will be in kilotrode
                    [~,ephysFolder] = fileparts(binFile.folder);
                    alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_%s.npy',timeRef,ephysFolder)));
                    if isempty(alignmentFileEphys)
                        %%% HACK BECAUSE TAGS ARE WEIRD SOMETIMES
                        % Get ephys tag
                        dEphys = dir(fullfile(server, subject, expDate, 'ephys*', '*')); 
                        dEphys(ismember({dEphys.name},{'.','..'})) = [];
                        dEphys(~[dEphys.isdir]) = [];
                        tags = {[]};
                        if numel(dEphys)>=1
                            for q = 1:numel(dEphys)
                                tags{q} = dEphys(q).name(numel(subject)+2:end);
                            end
                        end
                        tags = tags(cellfun(@(x) ~isempty(x), tags));
                        tag = tags(cell2mat(cellfun(@(x) contains(ephysFolder,x),tags,'uni',0)));
                        if numel(tag)>1
                            error('Too many ephys tags corresponding?')
                        end
                        if ~isempty(tag)
                            alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_ephys_%s.npy',timeRef,tag{1})));
                            if isempty(alignmentFileEphys)
                                % Sometime the "g0" has been removed at the end
                                alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_ephys_%s.npy',timeRef,regexprep(tag{1},'_g\d',''))));
                            end
                        end
                        % Can deal with that case only if contains natim in
                        % its name.
                        if isempty(alignmentFileEphys) && contains(lower(binFile.name),'natim')
                            % Check if several "natim" ephys... If yes,
                            % check that currently looking at the first. 
                            % If all good, consider that it's ok. 
                            alignmentFileEphys = dir(fullfile(alignmentFolder, sprintf('correct_timeline_%s_to_ephys_.npy',timeRef)));
                            if ~isempty(alignmentFileEphys)
                                fprintf('tricky situation -- there''s an undefined alignment file for this timeline...\n')
                                if contains(lower(binFile.name),'natim2')
                                    printf('File called natim2. Check if there''s another one?\n')
                                    if numel(unique({dEphys.folder}))==1
                                        cbinFilesForDate = dir(fullfile(dEphys(1).folder,'**','*.ap.cbin'));
                                    else
                                        cbinFilesForDate = dir(fullfile(server, subject, expDate,'**','*.ap.cbin')); % can't take ephys folder because there can be several...
                                    end
                                    if numel(cbinFilesForDate)>1
                                        if sum(contains({cbinFilesForDate.name},'natim','IgnoreCase',1))>1
                                            warning('Several natim, won''t align the current one. Check?')
                                            alignmentFileEphys = [];
                                        end
                                    end
                                else
                                    fprintf('Assume it''s the correct one.\n')
                                end
                            end
                        end
                    end
                    if ~isempty(alignmentFileEphys)
                        exp2keep = cat(1,exp2keep,{expFolder});
                    else
                        fprintf('Ephys wasn''t aligned to timeline %d.\n', timeRef)
                        continue
                    end
                end
            end
        end
    end

    exp2keep = unique(exp2keep);

    if isempty(exp2keep)
        warning('Couldn''t find any associated Natural Images experiment for %s.', binFile.name)
        if contains(lower(binFile.name),'natim')
            warning('But there''s ''natim'' in the name of the file. Double-check?')
        end
    end
end