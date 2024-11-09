function [clusinfo, sp, Params, recsesidxIncluded] = LoadPreparedClusInfo(KiloSortPaths,Params)

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(KiloSortPaths));
addthis=0;
StoreRecSesID = cell(1,length(KiloSortPaths));
for subsesid=1:length(KiloSortPaths)
     if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy'))) || ~exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'))
        continue
    end

    disp(['Loading clusinfo for ' KiloSortPaths{subsesid}])
    tmp = load(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'clusinfo');
    clusinfo{subsesid} = tmp.clusinfo;

    % Replace recsesid with subsesid
    clusinfo{subsesid}.RecSesID = clusinfo{subsesid}.RecSesID+addthis;
    StoreRecSesID{subsesid} = unique(clusinfo{subsesid}.RecSesID);
    addthis=max(clusinfo{subsesid}.RecSesID);
end



% Add all cluster information in one 'cluster' struct - can be used for further analysis
clusinfo = [clusinfo{:}];
if isempty(clusinfo)
    clusinfo = [];
    sp = [];
    return
end
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

%% Here we're going to actually load in all the sessions requested - sp
sp = cell(1,length(KiloSortPaths));
countid=1;
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))  || ~exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'))
        continue
    end
    disp(['Loading spike data for ' KiloSortPaths{subsesid}])
    tmp = load(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'sp');
    sp{subsesid} = tmp.sp;
    % Replace recsesid with subsesid
    sp{subsesid}.RecSes = repmat(countid,size(sp{subsesid}.RecSes));
    countid=countid+1;
end
% Add all spikedata in one spikes struct - can be used for further analysis
sp = [sp{:}];
spnew = struct;

fields = fieldnames(sp(1));
for fieldid=1:length(fields)
    try
        eval(['spnew.' fields{fieldid} '= cat(1,sp(:).' fields{fieldid} ');'])
    catch ME
        if strcmp(ME.message,'Out of memory.')
            eval(['spnew.' fields{fieldid} ' = sp(1).' fields{fieldid} ';'])
            for tmpid = 2:length(sp)
                eval(['spnew.' fields{fieldid} ' = cat(1,spnew.' fields{fieldid} ', sp(tmpid).' fields{fieldid} ');'])
            end
        else
            eval(['spnew.' fields{fieldid} '= cat(2,sp(:).' fields{fieldid} ');'])
        end
    end
end
sp = spnew;
sp.sample_rate = sp.sample_rate(1);
clear spnew

% Save out sp unique cluster
recsesidxIncluded = false(1,length(KiloSortPaths));
nclus = length(clusinfo.cluster_id);
RawPathsUsed = {};
sp.UniqClu = sp.clu;
sp.SpikeOnProbe = nan(size(sp.clu));
clusinfo.UniqueID = (1:length(clusinfo.cluster_id))';
clusinfo.IMROID = repmat(0,length(clusinfo.cluster_id),1);

if Params.UnitMatch
    disp('Assigning correct unique ID')
    PartsPath = strsplit(Params.SaveDir,'\');
    UMOutputAll = dir(fullfile(PartsPath{1},PartsPath{2},PartsPath{3},'**','UnitMatch','UnitMatch.mat'));
    for imroid = 1:length(UMOutputAll)
        IMROId = strsplit(UMOutputAll(imroid).folder,'\');
        IMROId = strsplit(IMROId{end-1},'_');
        IMROId = str2num(IMROId{end});
        if Params.separateIMRO & isempty(IMROId)
            continue
        end

        UMOutput = load(fullfile(UMOutputAll(imroid).folder,UMOutputAll(imroid).name),'UMparam','UniqueIDConversion');
        UMparam = UMOutput.UMparam;
        recsesidx = find(ismember(UMparam.KSDir,KiloSortPaths)); % Find which recses id they should have in UM output
        recsesidx2 = find(ismember(KiloSortPaths,UMparam.KSDir));
        if ~any(recsesidx)
            continue
        end
        UniqueIDConversion = UMOutput.UniqueIDConversion;
        
        TheseClus = find(ismember(clusinfo.RecSesID,[StoreRecSesID{recsesidx2}]));
        if ~isempty(TheseClus)
            clusinfo.UniqueID(TheseClus) = UniqueIDConversion.UniqueID(ismember(UniqueIDConversion.recsesAll,recsesidx))'; %Assign correct UniqueID
            clusinfo.IMROID(TheseClus) = repmat(IMROId,sum(ismember(clusinfo.RecSesID,[StoreRecSesID{recsesidx2}])),1);
        end
        for clusid=1:length(TheseClus)
            sp.UniqClu(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.UniqueID(clusid);
            sp.SpikeOnProbe(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.ProbeID(clusid);
        end
        recsesidxIncluded(recsesidx2) = 1;
        for rrid = 1:numel(recsesidx)
            if isstruct(UMparam.RawDataPaths{recsesidx(rrid)})
                UMparam.RawDataPaths{recsesidx(rrid)} = fullfile(UMparam.RawDataPaths{recsesidx(rrid)}.folder,UMparam.RawDataPaths{recsesidx(rrid)}.name);
            end
        end

        RawPathsUsed = {RawPathsUsed{:} UMparam.RawDataPaths{recsesidx}};
    end
end
Params.RawDataPaths = RawPathsUsed;
Params.KSDir = KiloSortPaths;
%% add nans for missing UID values to not confuse them with other UIDs
if any(recsesidxIncluded==0)
    clusinfo.UniqueID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
    clusinfo.IMROID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
else % Organize units; only one entry per UID to avoid double counts of neurons
    [UOpt,ID1,ID2] = unique(cat(2,clusinfo.UniqueID,clusinfo.IMROID),'rows','stable'); % Look for unique ID/IMRO pairs

    thesefields = fieldnames(clusinfo);
    for fid = 1:numel(thesefields)
        clusinfo.(thesefields{fid}) =  clusinfo.(thesefields{fid})(ID1);
    end
end




return