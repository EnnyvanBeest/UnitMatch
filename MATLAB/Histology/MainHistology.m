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

SaveFigs = SaveDir;
%% Automated
% Load all data
% Find available datasets (always using dates as folders)
clear DateOpt
DateOpt = arrayfun(@(X) dir(fullfile(DataDir{DataDir2Use(X)},MiceOpt{X},'*-*')),1:length(MiceOpt),'UniformOutput',0);

DateOpt = cellfun(@(X) X([X.isdir]),DateOpt,'UniformOutput',0);
DateOpt = cellfun(@(X) {X.name},DateOpt,'UniformOutput',0);
NewHistologyNeededOri = NewHistologyNeeded;

 % load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
end


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
    if length(unique(multidate))>1
        multidate=1;
    else
        multidate=0;
    end
%     if strcmp(RecordingType{midx},'Chronic') %& ~multidate  % These are my chronic mice, one dataset per mouse
%         %% Loading data from kilosort/phy easily
%         if ~MatchUnitsAcrossDays % Don't use pykilosort
%             myKsDir = fullfile(LocalDir,MiceOpt{midx},'Chronic');
%         else
%             myKsDir = fullfile(LocalDir,MiceOpt{midx});
%         end
%         subksdirs = dir(fullfile(myKsDir,'*','Probe*')); %This changed because now I suddenly had 2 probes per recording
%         if length(subksdirs)<1
%             clear subksdirs
%             subksdirs.folder = myKsDir; %Should be a struct array
%             subksdirs.name = 'Probe0';
%         end
%         ProbeOpt = unique({subksdirs(:).name});
%         for probeid = 1:length(ProbeOpt)
%             thisprobe = ProbeOpt{probeid}
%             myKsDir = (fullfile(subksdirs(1).folder,ProbeOpt{probeid}))
%             if ~isdir(myKsDir)
%                 continue
%             end
%             
%             %Saving directory
%             if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisprobe,'HistoEphysAlignment.mat'))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisprobe,'CuratedResults.mat')))
%                 disp([MiceOpt{midx} ' already aligned... '])
%                 if RedoUserInput
%                     tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisprobe,'HistoEphysAlignment.fig'));
%                     answer = questdlg(['Would you like to redo or keep this alignment? '  MiceOpt{midx} ' Probe' thisprobe], ...
%                         'Redo Alignment', ...
%                         'Redo','Keep','Keep');
%                     close(tmpfig)
% 
%                     if strcmp(answer,'Redo')
%                         NewHistologyNeeded=1;
%                     else
%                         disp('Skip')
%                         DrawProbeInBrain
%                         continue
%                     end
%                     
%                 else
%                     disp('Skip')
%                     DrawProbeInBrain
%                     continue
% 
%                 end
%             elseif RedoAfterClustering || NewHistologyNeeded
%                 myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
%                 if isempty(myClusFile)
%                     disp([MiceOpt{midx} ' is not yet curated with phy!!'])
%                 end
%                 NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering
%             end
%             
%             %% Get cluster information
%             myKsDir = fullfile(LocalDir,MiceOpt{midx},'*',ProbeOpt{probeid});
%             clear params
%             params.loadPCs=true;
%             thisdate = [];
%             PrepareClusInfo
%             
%             %% Get LFP?
%             myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},'*','ephys');
%             lfpD = dir(fullfile(myLFDir,'*','*','*.lf.*bin')); % ap file from spikeGLX specifically
%             if isempty(lfpD)
%                 disp('No LFP data found, maybe it lives in .ap file?')
%                 lfpD = dir(fullfile(myLFDir,'*','*','*.ap.*bin')); % ap file from spikeGLX specifically
%             end
%             if isempty(lfpD)
%                 disp('No, really no data found..')
%             elseif length(lfpD)>length(subksdirs)
%                 disp('Just take data from the last recording')
%                 lfpD = lfpD(end);
%             elseif length(lfpD)<length(subksdirs)
%                 disp('Should be a different amount of probes?')
%                 keyboard
%             else
%                 lfpD = lfpD(probeid);
%             end
%             
%             %% Get Histology output
%             GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
%             if ~histoflag
%                 disp([MiceOpt{midx} ' ' thisprobe 'No histology data, skip...'])
%                 continue
%             end
% 
%             %% Plot in atlas space
%            DrawProbeInBrain
% 
% 
%         end
%     else
        % For every date a different dataset
        Dates4Mouse = DateOpt{midx};
        for didx = 1:length(Dates4Mouse)
            % Within folders, look for 'RF mapping sessions'
            thisdate = Dates4Mouse{didx};
           
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
                
                %Saving directory
                thisprobe = subksdirs(probeid).name
                if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat'))) && (~RedoAfterClustering || exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat')))
                    if RedoUserInput
                        tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.fig'));
                        answer = questdlg(['Would you like to redo or keep this alignment? '  MiceOpt{midx} ' ' thisdate ' Probe' thisprobe], ...
                            'Redo Alignment', ...
                            'Redo','Keep','Keep');
                        close(tmpfig)
                        if strcmp(answer,'Redo')
                            NewHistologyNeeded=1;
                        else
                            disp('Skip')
                            DrawProbeInBrain
                            continue
                        end
                    else
                      

                        disp('Skip')
                        DrawProbeInBrain
                        continue
                    end
                elseif RedoAfterClustering || NewHistologyNeeded
                    myKsDir = fullfile(KilosortDir,MiceOpt{midx},thisdate,thisprobe);
                    myClusFile = dir(fullfile(myKsDir,'cluster_info.tsv'));
                    if isempty(myClusFile)
                        disp([MiceOpt{midx} ' ' thisdate 'is not yet curated with phy!!'])
                        if (~NewHistologyNeeded && exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.mat')))
                            if RedoUserInput
                                tmpfig = open(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'HistoEphysAlignment.fig'));
                                answer = questdlg(['Would you like to redo or keep this alignment? '  MiceOpt{midx} ' ' thisdate ' Probe' thisprobe], ...
                                    'Redo Alignment', ...
                                    'Redo','Keep','Keep');
                                close(tmpfig)
                                if strcmp(answer,'Redo')
                                    NewHistologyNeeded=1;
                                else
                                    disp('Skip')
                                    DrawProbeInBrain
                                    continue
                                end
                            else
                                disp('Skip')
                                DrawProbeInBrain
                                continue
                            end

                        else
                            NewHistologyNeeded = 1; %Automatically to 1 after RedoAfterClustering

                        end
                    end
                end
                
                
                %% Get cluster information
                clear params
                params.loadPCs=true;
                params.thisdate = thisdate;
                PrepareClusInfo
                
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
                 
                %% Get Histology output
                GetHistologyOutput      %I created an extra function to have one line of code in the different scripts
                if ~histoflag
                    disp([MiceOpt{midx} ' ' thisdate ' ' thisprobe 'No histology data, skip...'])
                    continue
                end
                
                NewHistologyNeeded = NewHistologyNeededOri;

                %% Plot in atlas space
              
                try
                DrawProbeInBrain
                catch ME
                    disp(ME)
                end
            end
            if strcmp(RecordingType{midx},'Chronic')
                break
            end
        end
%     end
end
% Save image
figure(fwireframe)
saveas(gcf,fullfile(SaveFigs,['3DProbes2.fig']))
saveas(gcf,fullfile(SaveFigs,['3DProbes2.svg']))

% set(fwireframe,'Units','normalized','Position',[0 0 1 1])

% spinningGIF(fullfile(SaveDir,'ProbeInsertionsAcrossMice.gif'))




