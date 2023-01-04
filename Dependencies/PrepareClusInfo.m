% Create saving directory
thisdate = params.thisdate;
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end

% Check for multiple subfolders?
subsesopt = dir(myKsDir);
subsesopt(~[subsesopt.isdir])=[];
subsesopt(ismember({subsesopt(:).name},{'..','.phy'}))=[];

if strcmp(RecordingType{midx},'Chronic')
    if ~RunPyKSChronic %MatchUnitsAcrossDays
        disp('Unit matching in Matlab')
        subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),{subsesopt(:).folder},'UniformOutput',0))) = []; %Use separate days and match units via matlab script
    else
        disp('Using chronic pyks option')
        subsesopt = subsesopt(cell2mat(cellfun(@(X) any(strfind(X,'Chronic')),{subsesopt(:).folder},'UniformOutput',0))); %Use chronic output from pyks
    end
end
%% Initialize everything
AllQMsPaths = cell(1,0);
AllRawPaths = cell(1,0);
AllUniqueTemplates = [];
recsesAll = [];
sp = cell(1,0);
channelmap=[];
channelpos = [];
cluster_id = [];
Label = [];
Good_ID = [];
depth = [];
channel = [];
Shank=[];
recses = [];
countid=1;
clear AllDecompPaths
% figure;
cols = jet(length(subsesopt));
for subsesid=1:length(subsesopt)
    if isempty(dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'*.npy')))
        continue
    end

    thissubses = str2num(subsesopt(subsesid).name)
    if isempty(thissubses)
        thissubses=1;
    end
    if isempty(thisdate)
        thisdatenow = strsplit(subsesopt(subsesid).folder,'\');
        thisdatenow = thisdatenow{end-1};
    else
        thisdatenow = thisdate;
    end
    %% Load Spike Data
    sp{countid} = loadKSdir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name),params); % Load Spikes with PCs
    [sp{countid}.spikeAmps, sp{countid}.spikeDepths, sp{countid}.templateDepths, sp{countid}.templateXpos, sp{countid}.tempAmps, sp{countid}.tempsUnW, sp{countid}.templateDuration, sp{countid}.waveforms] = templatePositionsAmplitudes(sp{countid}.temps, sp{countid}.winv, sp{countid}.ycoords, sp{countid}.xcoords, sp{countid}.spikeTemplates, sp{countid}.tempScalingAmps); %from the spikes toolbox

    %% Remove noise; spikes across all channels'
    figure
    scatter(sp{countid}.st,sp{countid}.spikeDepths,4,[0 0 0],'filled')
    hold on

    binsz = 0.001;
    edges = [min(sp{countid}.st)-binsz/2:binsz:max(sp{countid}.st)+binsz/2];
    timevec = [min(sp{countid}.st):binsz:max(sp{countid}.st)];

    depthstep = 25; %um
    depthedges = [min(sp{countid}.spikeDepths)-depthstep/2:depthstep:max(sp{countid}.spikeDepths)+depthstep/2];
    tmpact = arrayfun(@(X) histcounts(sp{countid}.st(sp{countid}.spikeDepths>depthedges(X)&sp{countid}.spikeDepths<depthedges(X+1)),edges),1:length(depthedges)-1,'UniformOutput',0);
    tmpact = cat(1,tmpact{:});
    tmpact = (tmpact-nanmean(tmpact,2))./nanstd(tmpact,[],2); %Z-score
    tpidx = find(sum(tmpact>3,1)>0.5*length(depthedges)-1);

    % Remove these timepoints
    rmidx = arrayfun(@(X) find(sp{countid}.st>timevec(X)-binsz/2&sp{countid}.st<timevec(X)+binsz/2),tpidx,'UniformOutput',0);
    rmidx = cat(1,rmidx{:});
    scatter(sp{countid}.st(rmidx),sp{countid}.spikeDepths(rmidx),4,[1 0 0],'filled')
    drawnow
    nori = length(sp{countid}.st);
    fields = fieldnames(sp{countid});
    for fid = 1:length(fields)
        eval(['tmp = sp{countid}. ' fields{fid} ';'])
        if any(size(tmp) == nori)
            tmp(rmidx,:,:)=[];
            eval(['sp{countid}.' fields{fid} '=tmp;'])
        end
    end

    %% Channel data
    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'channel_map.npy'));
    channelmaptmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));

    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'channel_positions.npy'));
    channelpostmp = readNPY(fullfile(myClusFile(1).folder,myClusFile(1).name));
    if length(channelmaptmp)<length(channelpostmp)
        channelmaptmp(end+1:length(channelpostmp))=length(channelmaptmp):length(channelpostmp)-1;
    end

    %% Is it correct channelpos though...?
    if strcmp(thisdatenow,'Chronic')
        sessionsIncluded = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'SessionsIncluded.mat'));
        sessionsIncluded = load(fullfile(sessionsIncluded.folder,sessionsIncluded.name));
        lfpD = arrayfun(@(X) dir(sessionsIncluded.ThesePaths{X}),1:length(sessionsIncluded.ThesePaths),'UniformOutput',0);
        lfpD = cat(2,lfpD{:});
    else
        myLFDir = fullfile(DataDir{DataDir2Use(midx)},MiceOpt{midx},thisdatenow,'ephys');
        lfpD = dir(fullfile([myLFDir],[thisdate '_' MiceOpt{midx} '*'], '**\*.ap.*bin')); % ap file from spikeGLX specifically
        if isempty(lfpD)
            lfpD = dir(fullfile([myLFDir],['*' MiceOpt{midx} '*'], '**\*.ap.*bin')); % ap file from spikeGLX specifically
        end
    end

    if isempty(lfpD)
        disp('No LFP data found')
    elseif length(lfpD)>=length(subsesopt)
        if ~strcmp(thisdatenow,'Chronic')
            probenr = cellfun(@(X) X(strfind(X,'imec')+4),{lfpD(:).name},'UniformOutput',0);
            if ~isempty(probenr)
                tmpprobe = strsplit(thisprobe,'Probe');
                disp('Take correct probe number data')
                lfpD = lfpD(( find(ismember(probenr,tmpprobe{2}))));
            else
                disp('Just take data from the last recording')
                lfpD = lfpD((thissubses));
            end

        end
    elseif length(lfpD)<length(subsesopt)
        disp('Should be a different amount of probes?')
        disp('Just take data from the last recording')
        lfpD = lfpD((thissubses));
    else
        lfpD = lfpD(thissubses);
    end

    if length(lfpD)>1
        if ~strcmp(thisdatenow,'Chronic')
            lfpD = lfpD(find(cellfun(@(X) any(strfind(X,['_' num2str(thissubses) '_'])),{lfpD(:).name})));
        end
    end
    channelpostmpconv = ChannelIMROConversion(lfpD(1).folder,1); % For conversion when not automatically done
    %% Load Cluster Info
    myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_info.tsv'));
    if isempty(myClusFile)
        disp('This data is not curated with phy!!')
        curratedflag=0;
        myClusFile = dir(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_group.tsv'));
        clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after phy)
        [clusinfo,sp{countid}] = ConvertTemplatesAfterPhy(clusinfo,sp{countid});

        clusidtmp = clusinfo.cluster_id;
        cluster_id = cat(1,cluster_id,clusinfo.cluster_id);
        tmpLabel = char(length(clusinfo.cluster_id));
        KSLabelfile = tdfread(fullfile(subsesopt(subsesid).folder,subsesopt(subsesid).name,'cluster_KSLabel.tsv'));
        tmpLabel(ismember(clusinfo.cluster_id,KSLabelfile.cluster_id)) = KSLabelfile.KSLabel(ismember(KSLabelfile.cluster_id,clusinfo.cluster_id));
        Label = [Label,tmpLabel];
        Good_IDtmp = ismember(tmpLabel,'g'); %Identify good clusters

        % Find depth and channel
        depthtmp = nan(length(clusinfo.cluster_id),1);
        xtmp = nan(length(clusinfo.cluster_id),1);
        channeltmp = nan(length(clusinfo.cluster_id),1);
        for clusid=1:length(depthtmp)
            % Depth according to PyKS2 output
            depthtmp(clusid)=round(sp{countid}.templateDepths(clusid));%round(nanmean(sp{countid}.spikeDepths(find(sp{countid}.clu==clusid-1))));
            xtmp(clusid)=sp{countid}.templateXpos(clusid);
            [~,minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1,channelpostmp(X,:),[xtmp(clusid),depthtmp(clusid)]),'euclidean'),1:size(channelpostmp,1),'UniformOutput',0)));
            try
                channeltmp(clusid) = channelmaptmp(minidx);
                depthtmp(clusid) = channelpostmpconv(minidx,2);
                xtmp(clusid) = channelpostmpconv(minidx,1);
            catch
                channeltmp(clusid)=channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmpconv(minidx-1,2);
                xtmp(clusid) = channelpostmpconv(minidx-1,1);
            end
            sp{countid}.spikeDepths(ismember(sp{countid}.clu,clusidtmp(clusid))) = depthtmp(clusid);

        end
        depth = cat(1,depth, depthtmp);
        channel = cat(1,channel,channeltmp);

    else
        CurationDone = 1;
        save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'CuratedResults.mat'),'CurationDone')
        clusinfo = tdfread(fullfile(myClusFile(1).folder,myClusFile(1).name));
        % Convert sp data to correct cluster according to phy (clu and
        % template are not necessarily the same after phy)
        [clusinfo,sp{countid}] = ConvertTemplatesAfterPhy(clusinfo,sp{countid});


        curratedflag=1;
        if isfield(clusinfo,'id')
            clusidtmp = clusinfo.id;
            cluster_id = [cluster_id,clusinfo.id];
        elseif isfield(clusinfo,'cluster_id')
            clusidtmp = clusinfo.cluster_id;

            cluster_id=[cluster_id,clusinfo.cluster_id];
        else
            keyboard
            disp('Someone thought it was nice to change the name again...')
        end
        KSLabel = clusinfo.KSLabel;
        Label = [Label,clusinfo.group]; % You want the group, not the KSLABEL!
        % Find depth and channel
        depthtmp = nan(length(clusidtmp),1);
        xtmp = nan(length(clusidtmp),1);
        channeltmp = nan(length(clusidtmp),1);
        for clusid=1:length(depthtmp)
            % Depth according to PyKS2 output
            depthtmp(clusid)=round(sp{countid}.templateDepths(clusid));%round(nanmean(sp{countid}.spikeDepths(find(sp{countid}.clu==clusid-1))));
            xtmp(clusid)=sp{countid}.templateXpos(clusid);
            [~,minidx] = min(cell2mat(arrayfun(@(X) pdist(cat(1,channelpostmp(X,:),[xtmp(clusid),depthtmp(clusid)]),'euclidean'),1:size(channelpostmp,1),'UniformOutput',0)));
            try
                channeltmp(clusid) = channelmaptmp(minidx);
                depthtmp(clusid) = channelpostmpconv(minidx,2);
                xtmp(clusid) = channelpostmpconv(minidx,1);
            catch
                channeltmp(clusid)=channelmaptmp(minidx-1);
                depthtmp(clusid) = channelpostmpconv(minidx-1,2);
                xtmp(clusid) = channelpostmpconv(minidx-1,1);
            end
            sp{countid}.spikeDepths(ismember(sp{countid}.clu,clusidtmp(clusid))) = depthtmp(clusid);

        end

        depth = [depth,depthtmp];

        channel = [channel,clusinfo.ch];
        Good_IDtmp = ismember(cellstr(clusinfo.group),'good')
        channeltmp = clusinfo.ch;
    end
    channelpostmp = channelpostmpconv;
    ypostmp = channelpostmp(:,2);
    xpostmp = channelpostmp(:,1);
    xposopt = (floor(xpostmp./250));% Assuming no new shank if not at least 100 micron apart
    Shanktmp = floor(xpostmp(channeltmp+1)./250);

    %     [~,minid] = arrayfun(@(X) (abs(floor(xpostmp(X)./250)-xposopt)),channeltmp+1,'UniformOutput',0);
    Shank = cat(1,Shank,Shanktmp);

    recses = cat(1, recses, repmat(countid,length(channeltmp),1));

    %     channelrecses = cat(2,channelrecses,repmat(countid,1,size(channelpostmp,1)));
    %     channelpos = cat(1,channelpos,channelpostmp);
    channelpos = channelpostmp; %only every need this here for mulitple sessions when chronic.
    channelmap = channelmaptmp;

    %     scatter(xtmp,depthtmp,10,cols(countid,:))
    %     hold on
    %     drawnow
    DecompressionFlag = 0;
    if RunQualityMatrix
        uniqueTemplates = [];
        unitTypeAcrossRec= [];
        %% Quality matrix - Bombcell (https://github.com/Julie-Fabre/bombcell)
        for id = 1:length(lfpD)
            ephysap_tmp = [];
            ephysap_path = fullfile(lfpD(id).folder,lfpD(id).name);
            if length(lfpD)>1
                savePath = fullfile(myClusFile(1).folder,num2str(id));
            else
                savePath =myClusFile(1).folder;
            end
            qMetricsExist = dir(fullfile(savePath, ['templates._jf_qMetrics.parquet']));
            idx = sp{countid}.SessionID==id;
            InspectionFlag = 0;
            if isempty(qMetricsExist) || RedoQM


                % First check if we want to use python for compressed data. If not, uncompress data first
                if any(strfind(lfpD(id).name,'cbin')) && DecompressLocal
                    if ~exist(fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','bin')))
                        disp('This is compressed data and we do not want to use Python integration... uncompress temporarily')
                        % Decompression
                        success = pyrunfile("MTSDecomp_From_Matlab.py","success",datapath = strrep(fullfile(lfpD(id).folder,lfpD(id).name),'\','/'),...
                            JsonPath =  strrep(fullfile(lfpD(id).folder,strrep(lfpD(id).name,'cbin','ch')),'\','/'), savepath = strrep(fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','bin')),'\','/'))
                        % Also copy metafile
                        copyfile(strrep(fullfile(lfpD(id).folder,lfpD(id).name),'cbin','meta'),strrep(fullfile(tmpdatafolder,lfpD(id).name),'cbin','meta'))
                    end
                    ephysap_tmp = fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','bin'));
                    DecompressionFlag = 1;
                    if InspectQualityMatrix
                        InspectionFlag=1;
                    end

                end
                clear param
                bc_qualityParamValues;
                % To make sure these parameters don't change (Github
                % pulls?)
                param.plotThis = 0;
                param.plotGlobal=1;
                param.verbose=1; % update user on progress
                param.reextractRaw=1; %Re extract raw waveforms
                param.rawWaveformMaxDef = 'firstSTD'
                %             idx = ismember(sp{countid}.spikeTemplates,clusidtmp(Good_IDtmp)); %Only include good units
                %careful; spikeSites zero indexed
                [qMetric, unitType] = bc_runAllQualityMetrics(param, round(sp{countid}.st(idx).*sp{countid}.sample_rate), sp{countid}.spikeTemplates(idx)+1, ...
                    sp{countid}.temps, sp{countid}.tempScalingAmps(idx),sp{countid}.pcFeat(idx,:,:),sp{countid}.pcFeatInd,channelpostmp,[], savePath);

            else
                ephysap_tmp = fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','bin'));

                load(fullfile(savePath, 'qMetric.mat'))
                bc_qualityParamValues;
                param.plotThis = 0;
                param.plotGlobal=1;
                bc_getQualityUnitType;
            end
            unitTypeAcrossRec{id} = unitType;
            uniqueTemplates{id} = unique(sp{countid}.spikeTemplates(idx));
            AllQMsPaths{countid}{id} = fullfile(savePath, 'qMetric.mat');
            AllRawPaths{countid}{id} = ephysap_path;
            AllDecompPaths{countid}{id} = ephysap_tmp;

            if InspectionFlag
                bc_getRawMemMap;

                %% view units + quality metrics in GUI
                % put ephys data into structure
                ephysData = struct;
                ephysData.spike_times = sp{countid}.st(idx).*round(sp{countid}.sample_rate);
                ephysData.ephys_sample_rate = sp{countid}.sample_rate;
                ephysData.spike_times_timeline = sp{countid}.st(idx);
                ephysData.spike_templates = sp{countid}.spikeTemplates(idx)+1; %0-indexed, make 1 indexed
                ephysData.templates = sp{countid}.temps;
                ephysData.template_amplitudes = sp{countid}.tempScalingAmps(idx);
                ephysData.channel_positions = channelpostmp;
                ephysData.waveform_t = 1e3*((0:size(sp{countid}.temps, 2) - 1) / ephysData.ephys_sample_rate);
                ephysParams = struct;
                plotRaw = 1;
                probeLocation=[];
                %% Inspect the results?
                unitQualityGuiHandle = bc_unitQualityGUI(memMapData,ephysData,qMetric, param, probeLocation, unitType, plotRaw);
                disp('Not continuing until GUI closed')
                while isvalid(unitQualityGuiHandle) && InspectQualityMatrix
                    pause(0.01)
                end
                disp('GUI closed')
            end
            %% Apply QM findings:
            disp('Only use units that passed the quality matrix parameters')
        end
        AllUniqueTemplates = cat(1,AllUniqueTemplates(:),cat(1,uniqueTemplates{:}));
        unitType = cat(1,unitTypeAcrossRec{:});
        Good_IDtmp(unitType ~= 1) = 0; % MUA
        Good_IDtmp(unitType == 1) =1; % Good
        Good_IDtmp(unitType == 0) =0;
        NoiseUnit = false(size(Good_IDtmp));
        NoiseUnit(unitType == 0)=1; % NOISE
    else
        AllUniqueTemplates = cat(1,AllUniqueTemplates,unique(sp{countid}.spikeTemplates));
        NoiseUnit = false(size(Good_IDtmp));
    end

    addthis = nanmax(recsesAll);
    if isempty(addthis)
        addthis=0;
    end
    if iscell(uniqueTemplates)
        recsesAlltmp = arrayfun(@(X) repmat(addthis+X,1,length(uniqueTemplates{X})),[1:length(uniqueTemplates)],'UniformOutput',0);
        recsesAll =cat(1,recsesAll(:), cat(2,recsesAlltmp{:})');
    else
        recsesAll = cat(1,recsesAll(:),repmat(addthis+1,1,length(uniqueTemplates))');
    end
    Good_ID = [Good_ID,Good_IDtmp]; %Identify good clusters



    sp{countid}.RecSes = sp{countid}.SessionID+addthis; %Keep track of recording session, as cluster IDs are not unique across sessions
    countid=countid+1;

    close all

end
sp = [sp{:}];
spnew = struct;

fields = fieldnames(sp(1));
for fieldid=1:length(fields)
    try
        eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
    catch
        eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])

    end
end
% Find correct dataset index
spikeTimes =cat(1,sp(:).st);
spikeRecSes = cat(1,sp(:).RecSes);
spikeSites = cat(1,sp(:).spikeTemplates);
spikeCluster = cat(1,sp(:).clu);
spikeAmps = cat(1,sp(:).spikeAmps);
spikeDepths = cat(1,sp(:).spikeDepths);
spikeShank = nan(length(spikeCluster),1);
ShankOpt = unique(Shank);
for shid = 1:length(ShankOpt)
    spikeShank(ismember(spikeCluster,cluster_id(Shank==ShankOpt(shid)))&ismember(spikeRecSes,recses(Shank==ShankOpt(shid)))) = ShankOpt(shid);
end
templateDepths = cat(1,sp(:).templateDepths);
tempAmps = cat(1,sp(:).tempAmps);
tempsUnW = cat(1,sp(:).tempsUnW);
templateDuration = cat(1,sp(:).templateDuration);
waveforms = cat(1,sp(:).waveforms);
templateWaveforms = cat(1,sp(:).temps);
templateAmplitudes = cat(1,sp(:).tempScalingAmps);
% pcFeatures = cat(1,sp(:).pcFeat);
% pcFeatureIdx = cat(1,sp(:).pcFeatInd);
channelPositionsX = cat(1,sp(:).xcoords);
channelPositionsY = cat(1,sp(:).ycoords);
channelPositions = cat(2,channelPositionsX,channelPositionsY);
%  [sp{countid}.spikeAmps, sp{countid}.spikeDepths, sp{countid}.templateDepths, sp{countid}.tempAmps, sp{countid}.tempsUnW, sp{countid}.templateDuration, sp{countid}.waveforms]
ShankOpt = unique(Shank);
ShankID = nan(size(Shank));
for shankid = 1:length(ShankOpt)
    ShankID(Shank==ShankOpt(shankid))=shankid;
end
clusinfo.Shank = Shank;
clusinfo.ShankID = ShankID;
clusinfo.RecSesID = recsesAll;
Good_IDx = find(Good_ID);
clusinfo.ch = channel;
clusinfo.depth = depth;
clusinfo.cluster_id = AllUniqueTemplates;
clusinfo.group = Label;
clusinfo.Good_ID = Good_ID;
clusinfo.Noise_ID = NoiseUnit;
AllQMsPaths = [AllQMsPaths{:}];
AllRawPaths = [AllRawPaths{:}];
if exist('AllDecompPaths') && ~isempty(AllDecompPaths)
    AllDecompPaths = [AllDecompPaths{:}];
else
    AllDecompPaths=[];
end
sp = spnew;
sp.sample_rate = sp.sample_rate(1);
clear spnew

%% Match different units?
 if MatchUnitsAcrossDays
     param.channelpos = channelpos;
     param.RunPyKSChronic = RunPyKSChronic;
     param.SaveDir = fullfile(SaveDir,MiceOpt{midx});
     [UniqueID, MatchTable] = UnitMatch(clusinfo,AllRawPaths,param,sp);
 end
 %% Remove temporary files
 if DecompressLocal && DecompressionFlag
     clear memMapData
     clear ap_data
     try
     delete(fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','bin')))
     delete(fullfile(tmpdatafolder,strrep(lfpD(id).name,'cbin','meta')))
     catch ME
         keyboard
     end
 end
