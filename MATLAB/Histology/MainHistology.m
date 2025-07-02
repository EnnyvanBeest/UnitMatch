%% User Input
NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
RedoAfterClustering = 0;
RedoUserInput = 0;
UseLFP = 0;
% directory of reference atlas files
ann = 10; %Steps in micron that's used
annotation_volume_location = 'annotation_volume_10um_by_index.npy';
structure_tree_location = 'structure_tree_safe_2017.csv';
% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% probe insertion direction 'down' (i.e. from the dorsal surface, downward -- most common!)
% or 'up' (from a ventral surface, upward)
probe_insertion_direction = 'down';

% show a table of regions that the probe goes through, in the console
show_region_table = true;

% black brain?
black_brain = false;
AllenCCFPath = fullfile(GithubDir,'allenCCF');
InclNrecordings = 8; % we typically have 8 mapping sessions at the start
%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0);
DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
NewHistologyNeededOri = NewHistologyNeeded;
RedoUserInputOri = RedoUserInput;

% load the reference brain annotations
if ~exist('av','var') || ~exist('st_allen','var')
    disp('loading reference atlas...')
    av = readNPY(fullfile(AllenCCFPath,annotation_volume_location));
    st_allen = readtable(fullfile(AllenCCFPath, structure_tree_location));
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
end

nInsertionsIncluded = 0;
nMiceIncluded = false(length(MiceOpt),1);
fwireframe = [];
% create a new figure with wireframe
fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on;
fwireframe.InvertHardcopy = 'off';

ProbeColors = distinguishable_colors(length(MiceOpt)+1);

for midx = 1:length(MiceOpt)

    figs = findall(0, 'Type', 'figure'); % Get all figures
    figs(figs == fwireframe) = []; % Remove the figure you want to keep
    close(figs); % Close the rest

    %% which probes?
    myKsDir = fullfile(KilosortDir,MiceOpt{midx});
    subksdirs = dir(fullfile(myKsDir,'*','Probe*')); %This changed because now I suddenly had 2 probes per recording
    multidate = cellfun(@(X) strsplit(X,'\'),{subksdirs(:).folder},'UniformOutput',0);
    multidate = cellfun(@(X) X{end},multidate,'UniformOutput',0);
    % if length(unique(multidate))>1 && strcmp(RecordingType{midx},'Chronic')
    %     multidate=1;
    % else
    multidate=0;
    % end
    if NewHistologyNeededOri || RedoUserInput
        thefirstperprobe = zeros(1,2);
    else
        thefirstperprobe = ones(1,2);
    end
    % For every date a different dataset
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)

        try
            if ~multidate
                NewHistologyNeeded = NewHistologyNeededOri;
                RedoUserInput = RedoUserInputOri;
            end

            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};

            tmpephysdir = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys',['*' MiceOpt{midx} '*']));
            if exist('IgnoreTheseFiles','var')
                for id = 1:length(IgnoreTheseFiles)
                    % Check if there's an ephys folder, if so run pyks2
                    tmpephysdir(find(cell2mat(cellfun(@(X) any(strfind(X,IgnoreTheseFiles{id})),{tmpephysdir(:).name},'UniformOutput',0))))=[];
                end
            end
            if isempty(tmpephysdir)
                continue
            end

            %% Loading data from kilosort/phy easily
            myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate);
            subksdirs = dir(fullfile(myKsDir,'Probe*')); %This changed because now I suddenly had 2 probes per recording
            if length(subksdirs)<1
                clear subksdirs
                subksdirs.folder = myKsDir; %Should be a struct array
                subksdirs.name = 'Probe0';
            end
            for probeid = 1:length(subksdirs)
                myKsDir = fullfile(subksdirs(probeid).folder,subksdirs(probeid).name)
                if ~isdir(myKsDir)
                    continue
                end
                RedoUserInput = RedoUserInputOri;
                NewHistologyNeeded = NewHistologyNeededOri;

                % Copy to other folders with same probe SN
                tmp = dir(fullfile(subksdirs(probeid).folder,subksdirs(probeid).name,'**','PreparedData.mat'));
                if ~isempty(tmp)
                    load(fullfile(tmp(1).folder,tmp(1).name),'SessionParams')
                    if isfield(SessionParams,'AllProbeSN')
                        SN = SessionParams.AllProbeSN{1};
                    else
                        rawD = SessionParams.RawDataPaths{1};
                        [channelpostmpconv, SN, recordingduration] = ChannelIMROConversion(rawD(1).folder, 0); % For conversion when not automatically done
                    end
                else
                    tmpfil = dir(fullfile(subksdirs(probeid).folder,subksdirs(probeid).name,'**','Params.py'));
                    spikeStruct = loadParamsPy(fullfile(tmpfil(1).folder,tmpfil(1).name));
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
                    [channelpostmpconv, SN, recordingduration] = ChannelIMROConversion(rawD(1).folder, 0); % For conversion when not automatically done
                end


                if multidate && thefirstperprobe(probeid)
                    RedoUserInput = 0;
                end

                %Saving directory
                thisprobe = subksdirs(probeid).name
                if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment_Auto.mat']))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat')))
                    if RedoUserInput
                        tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment_Auto.fig']));
                        answer = questdlg(['Would you like to redo or keep this alignment? '  MiceOpt{midx} ' ' thisdate ' Probe' thisprobe], ...
                            'Redo Alignment', ...
                            'Redo','Keep','Keep');
                        close(tmpfig)
                        if strcmp(answer,'Redo')
                            NewHistologyNeeded=1;
                        else
                            disp('Skip')
                            DrawProbeInBrain
                            nInsertionsIncluded = nInsertionsIncluded + 1;
                            nMiceIncluded(midx) = true;
                            thefirstperprobe(probeid) = 1;
                            continue
                        end
                    elseif ~multidate


                        disp('Skip')
                        try
                            DrawProbeInBrain
                        catch ME

                        end
                        nInsertionsIncluded = nInsertionsIncluded + 1;
                        nMiceIncluded(midx) = true;


                        continue
                    end
                elseif RedoAfterClustering || NewHistologyNeeded || ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment_Auto.mat']))
                    myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate,thisprobe);
                    myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                    if isempty(myClusFile)
                        disp([MiceOpt{midx} ' ' thisdate 'is not yet curated with phy!!'])
                        if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment_Auto.mat'])))
                            if RedoUserInput
                                tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment_Auto.fig']));
                                answer = questdlg(['Would you like to redo or keep this alignment? '  MiceOpt{midx} ' ' thisdate ' Probe' thisprobe], ...
                                    'Redo Alignment', ...
                                    'Redo','Keep','Keep');
                                close(tmpfig)
                                if strcmp(answer,'Redo')
                                    NewHistologyNeeded=1;
                                else
                                    disp('Skip')
                                    DrawProbeInBrain
                                    nInsertionsIncluded = nInsertionsIncluded + 1;
                                    nMiceIncluded(midx) = true;

                                    continue
                                end
                            else
                                disp('Skip')
                                DrawProbeInBrain
                                nInsertionsIncluded = nInsertionsIncluded + 1;
                                nMiceIncluded(midx) = true;

                                continue
                            end

                        else
                            NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
                        end
                    end
                end


                if multidate & thefirstperprobe(probeid) & NewHistologyNeeded
                    histfile = dir(fullfile(SaveDir,MiceOpt{midx},'**',[num2str(SN) '_HistoEphysAlignment_Auto.mat']));
                    if ~isempty(histfile)
                        if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
                            mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
                        end
                        copyfile(fullfile(histfile(1).folder,histfile(1).name),fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,histfile(1).name))
                        continue
                    end
                elseif multidate & thefirstperprobe(probeid) & ~NewHistologyNeeded
                    continue
                elseif multidate & ~thefirstperprobe(probeid) & ~NewHistologyNeeded
                    disp('Skip')
                    try
                        DrawProbeInBrain
                    catch ME

                    end
                    nInsertionsIncluded = nInsertionsIncluded + 1;
                    nMiceIncluded(midx) = true;
                    thefirstperprobe(probeid) = 1;

                    continue
                end
                thefirstperprobe(probeid) = 1;

                %% Get cluster information
                PipelineParams.thisdate = thisdate;
                PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe);
                if multidate
                    % Find all recording sessions
                    myKsDir = dir(fullfile(KilosortDir,MiceOpt{midx},thisdate,'**','spike_clusters.npy'));
                    IncludeThese = false(1,numel(myKsDir));

                    for ksid = 1:numel(myKsDir) % Only include the ones with the same SN for the probe
                        tmp = dir(fullfile(myKsDir(ksid).folder,'PreparedData.mat'));
                        if ~isempty(tmp)
                            tmp = load(fullfile(tmp(1).folder,tmp(1).name),'SessionParams');
                            if tmp.SessionParams.AllProbeSN{1} == SN
                                IncludeThese(ksid) = 1;
                            end
                        elseif any(strfind(myKsDir(ksid).folder,thisdate)) && any(strfind(myKsDir(ksid).folder,thisprobe))
                            IncludeThese(ksid) = 1;
                        end

                    end
                    myKsDir = myKsDir(IncludeThese);
                else
                    myKsDir = dir(fullfile(myKsDir,'**','spike_clusters.npy'));
                end
                myKsDir = arrayfun(@(X) X.folder,myKsDir,'UniformOutput',0);
                if numel(myKsDir)>InclNrecordings
                    myKsDir = myKsDir(1:InclNrecordings);
                end

                try
                    [clusinfo, sp, Params]  = LoadPreparedClusInfo(myKsDir,PipelineParams);
                catch ME
                    disp(ME)
                    PipelineParams = ExtractKilosortData(myKsDir,PipelineParams);
                    [clusinfo, sp, Params]  = LoadPreparedClusInfo(myKsDir,PipelineParams);
                end
                % This extracts the parameters within clusinfo and sp
                % struct for further analysis
                ExtractFields({sp,clusinfo})
                Good_IDx = find(Good_ID');
                if isempty(Good_IDx)
                    disp(['No good units for ' MiceOpt{midx} ' ' thisdate ' ' thisprobe ', continue'])
                    continue
                end

                %% Get LFP?
                if UseLFP
                    myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,'ephys');
                    lfpD = dir(fullfile([myLFDir '*'], '**\*.lf.*bin')); % lf file from spikeGLX specifically
                    if isempty(lfpD)
                        disp('No LFP data found')
                    elseif length(lfpD)~=length(subksdirs)
                        disp('Should be a different amount of probes?')
                        keyboard
                    else
                        lfpD = lfpD(probeid);
                    end

                    if isempty(lfpD)
                        % get AP data
                        lfpD = dir(fullfile([myLFDir '*'], '**\*.ap.*bin')); % lf file from spikeGLX specifically
                        if isempty(lfpD)
                            disp('No LFP data found')
                        elseif length(lfpD)~=length(subksdirs)
                            disp('Should be a different amount of probes?')
                            lfpD = lfpD(end);

                        else
                            lfpD = lfpD(end);
                        end
                    end
                else
                    lfpD = [];
                end


                %% Get Histology output
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts


                if ~histoflag
                    disp([MiceOpt{midx} ' ' thisdate ' ' thisprobe 'No histology data, skip...'])
                    continue
                end



                %% Plot in atlas space
                try
                    DrawProbeInBrain
                    nInsertionsIncluded = nInsertionsIncluded + 1;
                    nMiceIncluded(midx) = true;

                catch ME
                    disp(ME)
                end
            end
        catch ME
            disp(ME)
        end
    end
    %     end
end
% Save image
figure(fwireframe)
saveas(gcf,fullfile(SaveDir,['3DProbes2.fig']))
saveas(gcf,fullfile(SaveDir,['3DProbes2.svg']))

set(fwireframe,'Units','normalized','Position',[0 0 1 1])

disp([num2str(nInsertionsIncluded) ' insertions included in ' num2str(sum(nMiceIncluded)) ' mice'])

spinningGIF(fullfile(SaveDir,'ProbeInsertionsAcrossMice.gif'))



