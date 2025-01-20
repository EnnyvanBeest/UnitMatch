%% User Input
NewHistologyNeeded = 0; %Automatically to 1 after RedoAfterClustering
RedoAfterClustering=0;
RedoUserInput = 0;

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
    av = readNPY(annotation_volume_location);
    st_allen = loadStructureTree(structure_tree_location);
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

ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0];
if length(ProbeColors)<length(MiceOpt)
    ProbeColors = distinguishable_colors(length(MiceOpt));
end

for midx = 1:length(MiceOpt)
    %% which probes?
    myKsDir = fullfile(KilosortDir,MiceOpt{midx});
    subksdirs = dir(fullfile(myKsDir,'*','Probe*')); %This changed because now I suddenly had 2 probes per recording
    multidate = cellfun(@(X) strsplit(X,'\'),{subksdirs(:).folder},'UniformOutput',0);
    multidate = cellfun(@(X) X{end},multidate,'UniformOutput',0);
    if length(unique(multidate))>1 && strcmp(RecordingType{midx},'Chronic')
        multidate=1;
    else
        multidate=0;
    end
    % if NewHistologyNeededOri || RedoUserInput
        thefirstperprobe = zeros(1,2);
    % else

        % thefirstperprobe = ones(1,2);
    % end
    NewHistologyNeeded = NewHistologyNeededOri;
    RedoUserInput = RedoUserInputOri;
    % For every date a different dataset
    Dates4Mouse = DateOpt{midx};
    for didx = 1:length(Dates4Mouse)
        
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
            if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat']))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat')))
                if RedoUserInput
                    tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.fig']));
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
            elseif RedoAfterClustering || NewHistologyNeeded || ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat']))
                myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate,thisprobe);
                myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                if isempty(myClusFile)
                    disp([MiceOpt{midx} ' ' thisdate 'is not yet curated with phy!!'])
                    if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat'])))
                        if RedoUserInput
                            tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.fig']));
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


            if multidate & NewHistologyNeeded & thefirstperprobe(probeid)
                histfile = dir(fullfile(SaveDir,MiceOpt{midx},'**',[num2str(SN) '_HistoEphysAlignment.mat']));
                if ~isempty(histfile)
                    if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
                        mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
                    end
                    copyfile(fullfile(histfile(1).folder,histfile(1).name),fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,histfile(1).name))
                    continue
                end

            elseif multidate & ~thefirstperprobe(probeid)
                disp('Skip')
                try
                    DrawProbeInBrain
                catch ME

                end
                nInsertionsIncluded = nInsertionsIncluded + 1;
                nMiceIncluded(midx) = true;
                thefirstperprobe(probeid) = 1;

                continue
            elseif multidate & thefirstperprobe(probeid)
                continue
            end
            thefirstperprobe(probeid) = 1;

            %% Get cluster information
            PipelineParams.thisdate = thisdate;
            PipelineParams.SaveDir = fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe);
            if multidate
                myKsDir = dir(fullfile(KilosortDir,MiceOpt{midx},'**',thisprobe,'**','spike_clusters.npy'));
            else
                myKsDir = dir(fullfile(myKsDir,'**','spike_clusters.npy'));
            end
            myKsDir = arrayfun(@(X) X.folder,myKsDir,'UniformOutput',0);
            if numel(myKsDir)>30
                myKsDir = myKsDir(1:30);
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

            %% Get LFP?
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

          
            %% Get Histology output
            if numel(myKsDir) == 1
                myKsDir = myKsDir{1};
            elseif multidate
                disp('Multidate')
            else
                keyboard
            end
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



