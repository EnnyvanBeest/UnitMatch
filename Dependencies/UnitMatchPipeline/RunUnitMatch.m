function UMparam = RunUnitMatch(AllKiloSortPaths,Params)
%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(AllKiloSortPaths));
addthis=0;
for subsesid=1:length(AllKiloSortPaths)
    if isempty(dir(fullfile(AllKiloSortPaths{subsesid},'*.npy')))
        continue
    end

    disp(['Loading clusinfo for ' AllKiloSortPaths{subsesid}])
    tmp = matfile(fullfile(AllKiloSortPaths{subsesid},'PreparedData.mat'));
    clusinfo{subsesid} = tmp.clusinfo;

    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = clusinfo{subsesid}.RecSesID+addthis;
    addthis=max(clusinfo{subsesid}.RecSesID);
end

% Add all cluster information in one 'cluster' struct - can be used for further analysis
clusinfo = [clusinfo{:}];
clusinfoNew = struct;
fields = fieldnames(clusinfo(1));
for fieldid=1:length(fields)
    try
        eval(['clusinfoNew.' fields{fieldid} '= cat(1,clusinfo(:).' fields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.' fields{fieldid} '= cat(2,clusinfo(:).' fields{fieldid} ');'])
        catch ME
            keyboard
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

%% UnitMatch Parameters
% Use some bombcell parameters
if Params.RunQualityMetrics
    paramBC = bc_qualityParamValuesForUnitMatch;
    UMparam.sampleamount = paramBC.nRawSpikesToExtract; %500; % Nr. waveforms to include
    UMparam.spikeWidth =paramBC.spikeWidth; %82; % in sample space (time)
else
    UMparam.sampleamount = 1000; % 
    UMparam.spikeWidth = 82; %82; % in sample space (time)
end

UMparam.RunPyKSChronicStitched = Params.RunPyKSChronicStitched;
UMparam.SaveDir = fullfile(Params.SaveDir,'UnitMatch');
if ~isdir(UMparam.SaveDir)
    mkdir(UMparam.SaveDir)
end
UMparam.ACGbinSize =  0.001;%paramBC.ACGbinSize;
UMparam.ACGduration = 1;%paramBC.ACGduration;
UMparam.UseBombCelRawWav = Params.RunQualityMetrics; % 1 by default
UMparam.KSDir = AllKiloSortPaths;
UMparam.channelpos = Params.AllChannelPos;
UMparam.AllRawPaths = Params.RawDataPaths;
UMparam.AllDecompPaths = arrayfun(@(X) fullfile(Params.tmpdatafolder,strrep(Params.RawDataPaths(X).name,'cbin','bin')),1:length(Params.RawDataPaths),'Uni',0);
UMparam.RedoExtraction = 0; % Only necessary if KS was redone!
UMparam.ProbabilityThreshold = 0.8; %Typically where noise splits from real matches
UMparam.binsize = Params.binsz;
UMparam.Scores2Include = Params.Scores2Include; %
UMparam.ApplyExistingBayesModel = Params.ApplyExistingBayesModel; %If 1, use probability distributions made available by us
UMparam.MakePlotsOfPairs = Params.MakePlotsOfPairs; % Plot all pairs
UMparam.AssignUniqueID = Params.AssignUniqueID; %Assign Unique ID
UMparam.GoodUnitsOnly = Params.GoodUnitsOnly;
if Params.UnitMatch
    UnitMatchExist = dir(fullfile(UMparam.SaveDir,'UnitMatch.mat'));
    if ~isempty(UnitMatchExist) && ~Params.RedoUnitMatch
        load(fullfile(UMparam.SaveDir,'UnitMatch.mat'))
        clusinfo.UniqueID = UniqueIDConversion.UniqueID;
        clusinfo.MatchTable = MatchTable;
    else
        % Need to decompress if decompression wasn't done yet
        for id = 1:length(Params.RawDataPaths)
            ephysap_tmp = [];
            ephysap_path = fullfile(Params.RawDataPaths(id).folder,Params.RawDataPaths(id).name);
            if (isempty(UnitMatchExist) || Params.RedoUnitMatch) && ~Params.DecompressionFlag && ~UMparam.UseBombCelRawWav
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(Params.RawDataPaths(id).name,'cbin')) && Params.DecompressLocal
                    if ~exist(fullfile(Params.tmpdatafolder,strrep(Params.RawDataPaths(id).name,'cbin','bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')

                        decompDataFile = bc_extractCbinData(fullfile(Params.RawDataPaths(id).folder,Params.RawDataPaths(id).name),...
                            [], [], 0, fullfile(Params.tmpdatafolder,strrep(Params.RawDataPaths(id).name,'cbin','bin')));
                        paramBC.rawFile = decompDataFile;
                        % Decompression
                        %                         success = pyrunfile("MTSDecomp_From_Matlab.py","success",datapath = strrep(fullfile(RawDataPaths(id).folder,RawDataPaths(id).name),'\','/'),...
                        %                             JsonPath =  strrep(fullfile(RawDataPaths(id).folder,strrep(RawDataPaths(id).name,'cbin','ch')),'\','/'), savepath = strrep(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')),'\','/'));
                        %                         % Also copy metafile
                        copyfile(strrep(fullfile(Params.RawDataPaths(id).folder,Params.RawDataPaths(id).name),'cbin','meta'),strrep(fullfile(Params.tmpdatafolder,Params.RawDataPaths(id).name),'cbin','meta'))
                    end
                    ephysap_tmp = fullfile(Params.tmpdatafolder,strrep(Params.RawDataPaths(id).name,'cbin','bin'));
                    %                     DecompressionFlag = 1;
                end
            end
        end

        %% Run UnitMatch
        GlobalUnitMatchClock = tic;
        [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo,UMparam);
        UMrunTime = toc(GlobalUnitMatchClock);
        save(fullfile(UMparam.SaveDir,'UnitMatch.mat'),'UniqueIDConversion','MatchTable','WaveformInfo','UMparam','UMrunTime')
        DataSizeParam = CalculateDuration(Params.SaveDir);
        disp(['UnitMatch took ' num2str(round(UMrunTime/60*10)/10) 'minutes to run'])

    end
elseif Params.DecompressionFlag % You might want to at least save out averaged waveforms for every session to get back to later, if they were saved out by bomcell
    % Extract average waveforms
    ExtractAndSaveAverageWaveforms(clusinfo,UMparam)
end

return
