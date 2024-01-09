function [clusinfo, sp, Params, recsesidxIncluded] = LoadPreparedClusInfo(KiloSortPaths,Params)

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(KiloSortPaths));
addthis=0;
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy'))) || ~exist(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'))
        continue
    end

    disp(['Loading clusinfo for ' KiloSortPaths{subsesid}])
    tmp = load(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'),'clusinfo');
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
sp.UniqClu = sp.clu;
if Params.UnitMatch
    disp('Assigning correct unique ID')
    PartsPath = strsplit(Params.SaveDir,'\');
    UMOutputAll = dir(fullfile(PartsPath{1},PartsPath{2},PartsPath{3},'**','UnitMatch','UnitMatch.mat'));
    for imroid = 1:length(UMOutputAll)
        IMROId = strsplit(UMOutputAll(imroid).folder,'\');
        IMROId = strsplit(IMROId{end-1},'_');
        IMROId = str2num(IMROId{end});

        UMOutput = load(fullfile(UMOutputAll(imroid).folder,UMOutputAll(imroid).name),'UMparam','UniqueIDConversion');
        UMparam = UMOutput.UMparam;
        recsesidx = find(ismember(UMparam.KSDir,KiloSortPaths)); % Find which recses id they should have in UM output
        recsesidx2 = find(ismember(KiloSortPaths,UMparam.KSDir));
        if ~any(recsesidx)
            continue
        end
        UniqueIDConversion = UMOutput.UniqueIDConversion;
        
        TheseClus = find(ismember(clusinfo.RecSesID,recsesidx2));
        clusinfo.UniqueID(TheseClus) = UniqueIDConversion.UniqueID(ismember(UniqueIDConversion.recsesAll,recsesidx))'; %Assign correct UniqueID
        clusinfo.IMROID(TheseClus) = repmat(IMROId,sum(ismember(clusinfo.RecSesID,recsesidx2)),1);

        for clusid=1:length(TheseClus)
            sp.UniqClu(sp.clu==clusinfo.cluster_id(TheseClus(clusid)) & sp.RecSes==clusinfo.RecSesID(TheseClus(clusid))) = clusinfo.UniqueID(clusid);
        end
        recsesidxIncluded(recsesidx2) = 1;
    end
else
    clusinfo.UniqueID = (1:length(clusinfo.cluster_id))';
end

%% add nans for missing UID values
if any(recsesidxIncluded==0)
    clusinfo.UniqueID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
    clusinfo.IMROID(ismember(clusinfo.RecSesID,find(recsesidxIncluded==0))) = nan;
end



return