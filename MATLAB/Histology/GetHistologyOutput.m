histodone=0;
histoflag = 0;
Depth2Area = [];
AlignmentFeatures = [];
removenoise = 0;
AllenCCFPath = fullfile(GithubDir,'allenCCF');
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end
tmphistfile = dir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,['*HistoEphysAlignment_Auto.mat']));
if isempty(tmphistfile)
    tmphistfile = dir(fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdate,thisprobe,['*HistoEphysAlignment_Auto.mat']));
end
if ~isempty(tmphistfile) && ~NewHistologyNeeded
    tmpfile = load(fullfile(tmphistfile.folder,tmphistfile.name));
    try
        Depth2Area = tmpfile.Depth2Area;
        histodone=1;
        disp('Data already aligned to histology')
        histoflag=1;
    catch ME
        histodone=0;
    end
end
if ~histodone %Just in case it's not done yet
    if iscell(myKsDir)
        myKsDir = myKsDir{1};
    end
    histofile = dir(fullfile(myKsDir,'channel_locations.json'));
    if isempty(histofile)
        histofile = dir(fullfile(myKsDir,'*.csv'));
        if isempty(histofile)
            % Try the histology folder
            histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',thisdate,thisprobe,'*.csv'));
        end
        if isempty(histofile)
            % Maybe chronic?
            try
            myLFDir = dir(Params.RawDataPaths{1});
            MetaFileDir = dir(fullfile(myLFDir.folder,'*.meta'));
            ImecMeta =ReadMeta2(fullfile(MetaFileDir(1).folder,MetaFileDir(1).name));
            ProbeSN = ImecMeta.imDatPrb_sn;

            histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',['*' ProbeSN],'*.csv')); % Use Probe Serial Number is safest!
            if isempty(histofile)
                histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',thisprobe,'*.csv'));
            end
            catch ME
                disp(ME)

            end
        end
        if length(histofile)>1
            disp('Detecting multiple files... npix2 probe?')
        end
        if isempty(histofile)
            histofile = dir(fullfile(myKsDir,'probe_ccf.mat'));
            if isempty(histofile)
                disp('No channel locations known yet, do histology first')
                histoflag=0;
            else
                disp('This doesn''t work yet')
                keyboard
                % alignment with petersen probe output
                histoflag = 1;
                histinfo =load(fullfile(histofile(1).folder,histofile(1).name));
                fullfile(fullfile(histofile(1).folder,histofile(1).name))

                % Align ephys data with probe
                 [Depth2Area, AlignmentFeatures]  = alignatlasdata_automated(histinfo,AllenCCFPath,sp,clusinfo,1,1,2);

            end
        else % Automatic alignment with brain globe output
            histoflag = 1;
            histinfo = arrayfun(@(X) readtable(fullfile(histofile(X).folder,histofile(X).name),'ReadVariableNames',1,'Delimiter',','),1:length(histofile),'UniformOutput',0);

            % If available; find actual depth in brain
            % corresponding to this:
            trackcoordinates = arrayfun(@(X) readNPY(fullfile(histofile(X).folder,strrep(histofile(X).name,'.csv','.npy'))),1:length(histofile),'UniformOutput',0);
            for shid = 1:numel(histinfo)
                trackcoordinates{shid} = trackcoordinates{shid}(1:size(histinfo{shid},1),:);
            end
            % Align ephys data with probe
            if ~exist('HistologyAlignmentMethod','var'), HistologyAlignmentMethod = 'constrained'; end
            if strcmpi(HistologyAlignmentMethod,'automated')
                [Depth2Area, AlignmentFeatures] = alignatlasdata_automated(histinfo,AllenCCFPath,sp,clusinfo,removenoise,trackcoordinates);
                AlignmentFit = [];
            else
                % Constrained alignment (scale/offset with the track end as tip, limited local warp),
                % with LFP landmarks from all recordings used here (they may cover different shanks)
                LFP = [];
                RawForLFP = Params.RawDataPaths; if ~iscell(RawForLFP), RawForLFP = {RawForLFP}; end
                for rfi = 1:numel(RawForLFP)
                    rf = RawForLFP{rfi}; if isstruct(rf), rf = fullfile(rf(1).folder,rf(1).name); end
                    try
                        % running vs still periods from the mpep sessions in this recording (if saved)
                        if exist('SaveRFDir','var'), RunDir = SaveRFDir; else, RunDir = SaveDir; end
                        % running speed of all sessions in this recording (each its own Timeline):
                        % synced here if not saved yet (MainAnalysis_SpatialNavigation pipeline)
                        if exist('ExtractRunTiming','file') == 2 && exist('nidq_sync_used','var') && exist('DataDir2Use','var')
                            try
                                ExtractRunTiming(rf,DataDir,DataDir2Use,MiceOpt,midx,thisdate,thisprobe,nidq_sync_used,RunDir);
                            catch ME
                                warning('GetHistologyOutput:RunTiming','No running speed for %s: %s',rf,ME.message)
                            end
                        end
                        Run = RunningOnProbe(fullfile(RunDir,MiceOpt{midx},thisdate,thisprobe),rf);
                        L = LFPDepthFeatures(rf,struct('Run',Run));
                        fprintf('LFP %s: brain state from %s\n',rf,L.State.Source)
                        if isempty(LFP), LFP = L; elseif ~isempty(L)
                            LFP.Shank = [LFP.Shank; L.Shank]; LFP.Depth = [LFP.Depth; L.Depth]; LFP.Values = [LFP.Values; L.Values];
                            if isempty(LFP.RippleCSD), LFP.RippleCSD = L.RippleCSD; elseif ~isempty(L.RippleCSD), LFP.RippleCSD = [LFP.RippleCSD L.RippleCSD]; end
                            for fc = {'ThetaCSD','DentateCSD','SlowCSD'} % per-shank CSD structs of each recording
                                if isempty(LFP.(fc{1})), LFP.(fc{1}) = L.(fc{1}); elseif ~isempty(L.(fc{1})), LFP.(fc{1}) = [LFP.(fc{1}) L.(fc{1})]; end
                            end
                        end
                    catch ME
                        warning('GetHistologyOutput:LFP','No LFP features from %s: %s',rf,ME.message)
                    end
                end
                % Tuning landmarks from the passive stimulation analyses of this probe/day
                % (LoadSpikesRFBoot, Generalmpep_Analysis: run these before MainHistology)
                if exist('SaveRFDir','var'), TuningDir = SaveRFDir; else, TuningDir = SaveDir; end
                try
                    TuningLM = TuningLandmarks(fullfile(TuningDir,MiceOpt{midx},thisdate,thisprobe));
                catch ME
                    warning('GetHistologyOutput:Tuning','No tuning landmarks: %s',ME.message)
                    TuningLM = [];
                end
                % Visually evoked CSD (early sink: VIS layer 4 / superficial SC), from the mpep
                % trial onsets on the probe clock (StimTiming in MPEPData.mat)
                EvokedCSD = [];
                try
                    [CSDLM,EvokedCSD] = EvokedCSDLandmarks(fullfile(TuningDir,MiceOpt{midx},thisdate,thisprobe));
                    if ~isempty(CSDLM), TuningLM = [TuningLM(:); CSDLM(:)]'; end
                catch ME
                    warning('GetHistologyOutput:CSD','No visual CSD: %s',ME.message)
                end
                [Depth2Area, AlignmentFeatures, AlignmentFit] = alignatlasdata_constrained(histinfo,trackcoordinates,sp,clusinfo,AllenCCFPath, ...
                    struct('LFP',LFP,'Landmarks',TuningLM,'CSD',EvokedCSD));
            end
        end
    else
        histoflag=1;
        histinfo = fileread(fullfile(histofile(1).folder,histofile(1).name)); %Read json file
        histinfo = jsondecode(histinfo);% Decode json text

        % Make new type table of this
        AllCh = fieldnames(histinfo);
        histinfonew = cell(length(AllCh)-1,4);
        for id = 1:length(AllCh)
            if strcmp(AllCh{id},'origin')
                continue
            end
            eval(['histinfonew(id,1) = {' num2str(id) '};'])
            eval(['histinfonew(id,2) = {histinfo.' AllCh{id} '.brain_region_id};'])
            eval(['histinfonew(id,3) = {histinfo.' AllCh{id} '.brain_region};'])
            eval(['histinfonew(id,4) = {histinfo.' AllCh{id} '.brain_region};'])
        end
        histinfonew = cell2table(histinfonew);
        histinfonew.Properties.VariableNames = {'Position','RegionID','RegionAcronym','RegionName'};

        % Align ephys data with probe
        [Depth2Area, AlignmentFeatures] = alignatlasdata_automated(histinfonew,AllenCCFPath,sp,clusinfo,removenoise,0,fullfile(lfpD.folder,lfpD.name),2);
    end
end
if histoflag && ~histodone && isempty(Depth2Area)
    warning('GetHistologyOutput:NoAlignment','%s %s %s: automatic alignment failed, nothing saved',MiceOpt{midx},thisdate,thisprobe)
    histoflag = 0;
elseif histoflag && ~histodone
    OutBase = fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment']);
    % Keep the previous (possibly manually curated) alignment as *_beforeRedo: the new automatic one
    % replaces it. Only the first backup is kept, so running again does not overwrite the curated version.
    for ext = {'.mat','.fig','.bmp'}
        if exist([OutBase '_Auto' ext{1}],'file') && ~exist([OutBase '_beforeRedo' ext{1}],'file')
            movefile([OutBase '_Auto' ext{1}],[OutBase '_beforeRedo' ext{1}])
        end
    end
    set(gcf,'name',sprintf('%s %s %s',MiceOpt{midx},thisdate,thisprobe))
    saveas(gcf,[OutBase '_Auto.fig'])
    saveas(gcf,[OutBase '_Auto.bmp'])

    ManuallyCurated  = false;
    FlaggedForReview = false;
    if ~exist('AlignmentFit','var'), AlignmentFit = []; end
    % A tuning/LFP landmark without its areas on the track: check the track in the curation GUI
    if isstruct(AlignmentFit) && isfield(AlignmentFit,'LandmarkWarnings') && ~isempty(AlignmentFit.LandmarkWarnings)
        FlaggedForReview = true;
    end
    save([OutBase '_Auto.mat'],'Depth2Area','AlignmentFeatures','AlignmentFit','ManuallyCurated','FlaggedForReview')
end

% Add information to clusinfo

if exist('tmp') && isstruct(tmp) && isfield(tmp,'VRDat') && ~isfield(tmp,'clusinfo')
    myKsDir = dir(fullfile(SaveDir,tmp.VRDat.Mouse(1,:),tmp.VRDat.Date(1,:),tmp.VRDat.Session(1,:),thisprobe,'SpikeData.mat'));
    tmpload = load(fullfile(myKsDir.folder,myKsDir.name));
    clusinfo = tmpload.clusinfo;

    myKsDir = dir(fullfile(KilosortDir,tmp.VRDat.Mouse(1,:),tmp.VRDat.Date(1,:),thisprobe,'**','PreparedData.mat'));
    if numel(myKsDir)>1
        mostlikely = zeros(1,numel(myKsDir));
        for kid = 1:numel(myKsDir)
            tmpspl = strsplit(myKsDir(kid).folder,filesep);
            if strcmp(tmpspl{end},'s')
                mostlikely(kid) = 1;
            end                
        end
        if sum(mostlikely)==1
            myKsDir = myKsDir(find(mostlikely));
        elseif sum(mostlikely) == 0
            for ksid = 1:numel(myKsDir)
                tmpload = load(fullfile(myKsDir(ksid).folder,myKsDir(ksid).name),'SessionParams');
                PipelineParams.SaveDir = tmpload.SessionParams.SaveDir;
                [clusinfo, ~, ~]  = LoadPreparedClusInfo({myKsDir(ksid).folder},PipelineParams);

                if all(ismember(tmp.VRDat.clusid,clusinfo.cluster_id(logical(clusinfo.Good_ID))))
                    mostlikely(ksid)=sum(ismember(clusinfo.cluster_id(logical(clusinfo.Good_ID)),tmp.VRDat.clusid))./numel(clusinfo.cluster_id(logical(clusinfo.Good_ID)));
                end
            end
            [~,maxid] = max(mostlikely);
            myKsDir = myKsDir(maxid);
        end
    end
    tmpload = load(fullfile(myKsDir.folder,myKsDir.name),'SessionParams');
    PipelineParams.SaveDir = tmpload.SessionParams.SaveDir;
    [clusinfo, ~, ~]  = LoadPreparedClusInfo({myKsDir.folder},PipelineParams);
end
if exist('clusinfo','var') & ~isempty(Depth2Area)
    rgbColors = cell2mat(cellfun(@hex2rgb, Depth2Area.Color, 'Uni', 0));
    grayMask = std(rgbColors, 0, 2) < 0.05;  % threshold can be adjusted

    % Calculate all distances
    if all(isnan(clusinfo.AdjustedDepth))
        clusinfo.AdjustedDepth(isnan(clusinfo.AdjustedDepth)) = clusinfo.depth; % Just take normal depth
    end
    % Shank number (0-based, as clusinfo.Shank) of each Depth2Area row. Depth2Area.Shank is the position
    % of the shank among the shanks with units in the recordings the alignment was made with
    % (alignatlasdata_automated), not the shank number; newer alignments store ShankAbs. All older
    % alignments checked (EB036-EB055, 16 probes) were made with units on all shanks, so there
    % shank number = Shank - 1. (Comparing Shank directly gave units on shanks 1-3 the areas of the
    % neighbouring shank.)
    if ismember('ShankAbs',Depth2Area.Properties.VariableNames)
        D2AShank = Depth2Area.ShankAbs;
    else
        D2AShank = Depth2Area.Shank - 1;
    end
    NoAlignment = setdiff(unique(clusinfo.Shank),unique(D2AShank));
    if ~isempty(NoAlignment)
        warning('GetHistologyOutput:ShankWithoutAlignment','%s %s %s: no histology alignment for shank(s) %s: their units get the areas of the nearest aligned shank (see clusinfo.AreaDist)', ...
            MiceOpt{midx},thisdate,thisprobe,mat2str(NoAlignment(:)'))
    end
    allDists = abs(clusinfo.AdjustedDepth - Depth2Area.Depth') + abs(clusinfo.Shank - D2AShank') * 1000;

    % Set distances to gray areas to Inf
    allDists(:, grayMask) = Inf;

    % Get closest non-gray region
    [mindist, DistIdx] = nanmin(allDists, [], 2);
    clusinfo.AreaDist = mindist; % um to the assigned histology point (+1000 per shank away; large: outside the aligned range or in white matter)

    % Assign
    clusinfo.Area = lower(Depth2Area.Area(DistIdx));
    clusinfo.Color = cell2mat(arrayfun(@(X) hex2rgb(Depth2Area.Color{X}), DistIdx, 'Uni', 0));
    clusinfo.Coordinates = Depth2Area.Coordinates(DistIdx, :);



    %% I only care about certain levels of histology - not layer specificity at this point
    % Group some areas
    Area= strrep(clusinfo.Area,'/','');

    [areaopt, id1] = unique(Area,'stable');
    atlastable = readtable('structure_tree_safe_2017.csv');
    if ~isfield(clusinfo,'AreaFN')
        clusinfo.AreaFN = cellfun(@(X) atlastable.name(find(ismember(lower(atlastable.acronym),X))),clusinfo.Area,'Uni',0)
    end
    areaoptFN = clusinfo.AreaFN(id1);
    uniquearean = length(areaopt);
    % Use Allen Brain Atlas to find all areas in the recordings
    % Add units to area specific cell
    % Area = Area(~cellfun(@isempty,Area));
    AutomaticAREASOfInterest = {};
    AutomaticAREASOfInterestFullName = {};
    ccfIdx = {};
    ColPerFullArea = nan(0,3);
    automaticareasofinterestid = nan(1,uniquearean); %Index for which area
    newareaname = cell(1,length(Area));
    newareaabrev = cell(1,length(Area));
    % Make names in table similar to those used here
    atlastable.acronym = lower(atlastable.acronym);
    atlastable.acronym = strrep(atlastable.acronym,'/','');
    for areaid=1:length(areaopt)
        %Find structure_id_path of area
        structure_id_path = atlastable.structure_id_path{find(ismember(lower(atlastable.acronym),areaopt{areaid}))};
        % Cut off last part
        parts = strsplit(structure_id_path,'/');
        if numel(parts)>3
            newstructure_id_path = fullfile(parts{2:end-2});
            newstructure_id_path = ['/' strrep(newstructure_id_path,'\','/') '/'];
        else
            newstructure_id_path = structure_id_path;
        end



        %     parts(cellfun(@isempty,parts))=[];
        if (length(parts)<=4 || sum(contains(atlastable.structure_id_path,newstructure_id_path)) > 17 || any(strfind(areaopt{areaid},'ca'))) && ~contains(areaoptFN{areaid},'layer')%if we cannot go further up or hippocampal area
            tmpareaname = areaopt{areaid}; %if we cannot go further up
            newstructure_id_path = structure_id_path;
        else
                % Find area with this new structure id
            tmpareaname = atlastable.acronym{find(ismember(atlastable.structure_id_path,newstructure_id_path))};
        end

        % Group auditory cortex regions together (audp/audd/audv/audpo/etc -> aud).
        [tmpareaname, forcedFullName] = local_groupAuditoryAcronym(tmpareaname);
        %Does this one already exist in Automatic AREAS Of Interest?
        if ~any(ismember(AutomaticAREASOfInterest,tmpareaname))
            %no, then create
            AutomaticAREASOfInterest = {AutomaticAREASOfInterest{:} tmpareaname};
            if ~isempty(forcedFullName)
                AutomaticAREASOfInterestFullName = {AutomaticAREASOfInterestFullName{:} forcedFullName};
            else
                AutomaticAREASOfInterestFullName = {AutomaticAREASOfInterestFullName{:} atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}};
            end
            % Use average color scheme
            ColPerFullArea = cat(1,ColPerFullArea,nanmean(clusinfo.Color(ismember(Area,areaopt{areaid}),:),1));
            ccfIdx = {ccfIdx{:} find(ismember(lower(atlastable.acronym),areaopt{areaid}))};

        else
            ccfIdx{find(ismember(AutomaticAREASOfInterest,tmpareaname))} = [ccfIdx{find(ismember(AutomaticAREASOfInterest,tmpareaname))} find(ismember(lower(atlastable.acronym),areaopt{areaid}))];
        end
        % Index correctly
        automaticareasofinterestid(areaid) = find(ismember(AutomaticAREASOfInterest,tmpareaname));
        newareaabrev(ismember(Area,areaopt{areaid})) = {tmpareaname};
        if ~isempty(forcedFullName)
            newareaname(ismember(Area,areaopt{areaid}))={forcedFullName}; %Save out per unit
        else
            newareaname(ismember(Area,areaopt{areaid}))={atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}}; %Save out per unit
        end
    end
    if any(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'})) % make these void
        AutomaticAREASOfInterestFullName(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'}))={'root'};
        AutomaticAREASOfInterest(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'}))={'root'};
    end
    [~,colidx] = ismember(newareaname,AutomaticAREASOfInterestFullName);
    [~,sortareaidx] = sortrows(ColPerFullArea);

    clusinfo.OriArea = Area;
    clusinfo.Area = newareaabrev';
    clusinfo.AreaFN = newareaname';

end

if exist('tmp') && isstruct(tmp) && isfield(tmp,'VRDat')
    clustershere = tmp.VRDat.clusid;
    for cid = 1:numel(clustershere)
        idx = find(clusinfo.cluster_id == tmp.VRDat.clusid(cid),1,'first');
        tmp.VRDat.Area(cid) = clusinfo.Area(idx);
        tmp.VRDat.Color(cid,:) = clusinfo.Color(idx,:);
        tmp.VRDat.Coordinates(cid,:) = clusinfo.Coordinates(idx,:);
    end

end

function [acrOut, forcedFullName] = local_groupAuditoryAcronym(acrIn)
% Group auditory cortical regions together.
% Allen acronyms typically include AUDp/AUDd/AUDv/AUDpo etc.
acrOut = acrIn;
forcedFullName = '';
if isempty(acrIn)
    return
end
if isstring(acrOut)
    acrOut = char(acrOut);
end
acrLow = lower(strtrim(acrOut));
if startsWith(acrLow, 'aud')
    acrOut = 'aud';
    forcedFullName = 'Auditory cortex';
end
end
