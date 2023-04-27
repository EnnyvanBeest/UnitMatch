function [clusinfo, sp, Params] = LoadPreparedClusInfo(KiloSortPaths,Params)

%% Here we're going to actually load in all the sessions requested - only clusinfo to save memory for unitmatch
clusinfo = cell(1,length(KiloSortPaths));
addthis=0;
for subsesid=1:length(KiloSortPaths)
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
        continue
    end
   
    disp(['Loading clusinfo for ' KiloSortPaths{subsesid}])
    tmp = matfile(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'));
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
    if isempty(dir(fullfile(KiloSortPaths{subsesid},'*.npy')))
        continue
    end
    disp(['Loading spike data for ' KiloSortPaths{subsesid}])
    tmp = matfile(fullfile(KiloSortPaths{subsesid},'PreparedData.mat'));
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
nclus = length(clusinfo.cluster_id);
sp.UniqClu = sp.clu;
if Params.UnitMatch
    for clusid=1:nclus
        sp.UniqClu(sp.clu==clusinfo.cluster_id(clusid) & sp.RecSes==clusinfo.RecSesID(clusid)) = clusinfo.UniqueID(clusid);
    end
end


return