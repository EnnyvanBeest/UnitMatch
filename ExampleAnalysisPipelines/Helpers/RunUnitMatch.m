function [UMparam, UniqueIDConversion, MatchTable, WaveformInfo] = RunUnitMatch(Params)


%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
AllKiloSortPaths = Params.KSDir;
clusinfo = cell(1, length(AllKiloSortPaths));
addthis = 0;
for subsesid = 1:length(AllKiloSortPaths)
    if isempty(dir(fullfile(AllKiloSortPaths{subsesid}, '*.npy')))
        continue
    end

    disp(['Loading clusinfo for ', AllKiloSortPaths{subsesid}])
    tmp = matfile(fullfile(AllKiloSortPaths{subsesid}, 'PreparedData.mat'));
    clusinfo{subsesid} = tmp.clusinfo;

    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = clusinfo{subsesid}.RecSesID + addthis;
    addthis = max(clusinfo{subsesid}.RecSesID);
end

% Add all cluster information in one 'cluster' struct - can be used for further analysis
% Apparently sometimes clusinfo does not have the same fields (KS version?)
fields = fieldnames(clusinfo{1});
IncludeField = true(1,length(fields));
for sesid=2:length(clusinfo)
    IncludeField(~ismember(fields,fieldnames(clusinfo{sesid}))) = false;
end
for sesid=1:length(clusinfo)
    thesefields = fieldnames(clusinfo{sesid});
    clusinfo{sesid} = rmfield(clusinfo{sesid},thesefields(ismember(thesefields,fields(find(~IncludeField))))); % REmove fields
end

% Continue adding them together
clusinfo = [clusinfo{:}];
clusinfoNew = struct;
fields = {fields{find(IncludeField)}}';
for fieldid = 1:length(fields)
    try
        eval(['clusinfoNew.', fields{fieldid}, '= cat(1,clusinfo(:).', fields{fieldid} ');'])
    catch ME
        try
            eval(['clusinfoNew.', fields{fieldid}, '= cat(2,clusinfo(:).', fields{fieldid}, ');'])
        catch ME
            keyboard
        end
    end
end
clusinfo = clusinfoNew;
clear clusinfoNew

if sum(clusinfo.Good_ID)<25
    disp('Less than 25 neurons.. skip')
    UMparam.Error = 'Too Little Units'
    return
end


%% UnitMatch Parameters
% Use some bombcell parameters
if Params.RunQualityMetrics
    paramBC = bc_qualityParamValuesForUnitMatch;
    UMparam.sampleamount = paramBC.nRawSpikesToExtract; %500; % Nr. waveforms to include
    UMparam.spikeWidth = paramBC.spikeWidth; %82; % in sample space (time)
else
    UMparam.sampleamount = 1000; %
    UMparam.spikeWidth = 82; %82; % in sample space (time)
end

UMparam.RunPyKSChronicStitched = Params.RunPyKSChronicStitched;
UMparam.SaveDir = fullfile(Params.SaveDir, 'UnitMatch');
if ~isdir(UMparam.SaveDir)
    mkdir(UMparam.SaveDir)
end
UMparam.ACGbinSize = 0.001; %paramBC.ACGbinSize;
UMparam.ACGduration = 1; %paramBC.ACGduration;
UMparam.UseBombCelRawWav = Params.RunQualityMetrics; % 1 by default
UMparam.KSDir = AllKiloSortPaths;
UMparam.channelpos = Params.AllChannelPos;
UMparam.AllRawPaths = Params.RawDataPaths;
UMparam.UseHistology = Params.UseHistology;
if isstruct(Params.RawDataPaths{1})
    if length(Params.RawDataPaths)==1
        UMparam.AllDecompPaths = arrayfun(@(X) fullfile(Params.tmpdatafolder, strrep(X.name, 'cbin', 'bin')), Params.RawDataPaths{1}, 'Uni', 0);;
    else
        UMparam.AllDecompPaths = cellfun(@(X) fullfile(Params.tmpdatafolder, strrep(X.name, 'cbin', 'bin')), Params.RawDataPaths, 'Uni', 0);
    end
else
    allDecompPaths_dirs = arrayfun(@(X) dir(Params.RawDataPaths{X}), 1:length(Params.RawDataPaths), 'Uni', 0);
    UMparam.AllDecompPaths = arrayfun(@(X) fullfile(Params.tmpdatafolder, strrep(allDecompPaths_dirs{X}.name, 'cbin', 'bin')), 1:length(allDecompPaths_dirs), 'Uni', 0);
end
UMparam.RedoExtraction = 0; % Only necessary if KS was redone!
UMparam.ProbabilityThreshold = 0.5; %Increase a bit for being more strict
UMparam.binsize = Params.binsz;
UMparam.Scores2Include = Params.Scores2Include; %
UMparam.ApplyExistingBayesModel = Params.ApplyExistingBayesModel; %If 1, use probability distributions made available by us
UMparam.MakePlotsOfPairs = Params.MakePlotsOfPairs; % Plot all pairs
if isfield(Params,'GUI')
    UMparam.GUI = Params.GUI; % Open GUI for manual curation of pairs
else
    UMparam.GUI = 0;
end

UMparam.AssignUniqueID = Params.AssignUniqueID; %Assign Unique ID
UMparam.GoodUnitsOnly = Params.GoodUnitsOnly;
if isfield(Params, 'RecType') && strcmpi(Params.RecType, 'Acute')
    UMparam.ExpectMatches = 0; %We do actually not expect any matches. Change prior accordingly
else
    UMparam.ExpectMatches = 1;
end
if Params.UnitMatch
    UnitMatchExist = dir(fullfile(UMparam.SaveDir, 'UnitMatch.mat'));
    if ~isempty(UnitMatchExist) && ~Params.RedoUnitMatch
        load(fullfile(UMparam.SaveDir, 'UnitMatch.mat'))
        clusinfo.UniqueID = UniqueIDConversion.UniqueID;
        clusinfo.MatchTable = MatchTable;
    else
        % Need to decompress if decompression wasn't done yet
        for id = 1:length(Params.RawDataPaths)
            ephysap_tmp = [];
            if isstruct(Params.RawDataPaths{id})
                ephysap_path = fullfile(Params.RawDataPaths{id}.folder, Params.RawDataPaths{id}.name);
            else
                ephysap_path = fullfile(Params.RawDataPaths{id});
            end
            if (isempty(UnitMatchExist) || Params.RedoUnitMatch) && ~Params.DecompressionFlag && ~UMparam.UseBombCelRawWav
                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(Params.RawDataPaths(id).name, 'cbin')) && Params.DecompressLocal
                    if ~exist(fullfile(Params.tmpdatafolder, strrep(Params.RawDataPaths(id).name, 'cbin', 'bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')

                        decompDataFile = bc_extractCbinData(fullfile(Params.RawDataPaths(id).folder, Params.RawDataPaths(id).name), ...
                            [], [], 0, fullfile(Params.tmpdatafolder, strrep(Params.RawDataPaths(id).name, 'cbin', 'bin')));
                        paramBC.rawFile = decompDataFile;
                        % Decompression
                        %                         success = pyrunfile("MTSDecomp_From_Matlab.py","success",datapath = strrep(fullfile(RawDataPaths(id).folder,RawDataPaths(id).name),'\','/'),...
                        %                             JsonPath =  strrep(fullfile(RawDataPaths(id).folder,strrep(RawDataPaths(id).name,'cbin','ch')),'\','/'), savepath = strrep(fullfile(Params.tmpdatafolder,strrep(RawDataPaths(id).name,'cbin','bin')),'\','/'));
                        %                         % Also copy metafile
                        copyfile(strrep(fullfile(Params.RawDataPaths(id).folder, Params.RawDataPaths(id).name), 'cbin', 'meta'), strrep(fullfile(Params.tmpdatafolder, Params.RawDataPaths(id).name), 'cbin', 'meta'))
                    end
                    ephysap_tmp = fullfile(Params.tmpdatafolder, strrep(Params.RawDataPaths(id).name, 'cbin', 'bin'));
                    %                     DecompressionFlag = 1;
                end
            end
        end

        %% Run UnitMatch
        [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
        tmpfile = dir(fullfile(UMparam.SaveDir,'UnitMatch.mat'));
        if UMparam.AssignUniqueID
            AssignUniqueID(fullfile(tmpfile.folder,tmpfile.name));
        end

    end
elseif Params.DecompressionFlag % You might want to at least save out averaged waveforms for every session to get back to later, if they were saved out by bomcell
    % Extract average waveforms
    ExtractAndSaveAverageWaveforms(clusinfo, UMparam)
end

return