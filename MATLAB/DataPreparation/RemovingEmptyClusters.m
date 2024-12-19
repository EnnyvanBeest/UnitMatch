function   [clusinfo, sp, emptyclus] = RemovingEmptyClusters(clusinfo,sp)

%% Due to merging/ splitting etc. there are empty clusters. Remove these
% Convert sp data accordingly

% Initialize new vectors
nclus = length(clusinfo.cluster_id);
temps = nan(nclus,size(sp.temps,2),size(sp.temps,3));
pcFeatInd = nan(nclus,size(sp.pcFeatInd,2));
templateDepths = nan(nclus,1);
templateXpos = nan(nclus,1);
tempAmps = nan(nclus,1);
tempsUnW = nan(nclus,size(sp.tempsUnW,2),size(sp.tempsUnW,3));
templateDuration = nan(nclus,1);
waveforms = nan(nclus,size(sp.waveforms,2));

emptyclus = [];
for clusid=1:nclus
    oriclusid = unique(sp.clu(find(sp.clu == clusinfo.cluster_id(clusid))));
    if isempty(oriclusid)
        emptyclus = [emptyclus clusid];
    elseif numel(oriclusid)==1
        temps(clusid,:,:)=sp.temps(oriclusid+1,:,:);
        try
            pcFeatInd(clusid,:)=sp.pcFeatInd(oriclusid+1,:);
        catch ME

        end
        templateDepths(clusid) = sp.templateDepths(oriclusid+1); % 0-indexed!
        templateXpos(clusid) = sp.templateXpos(oriclusid+1);
        tempAmps(clusid) = sp.tempAmps(oriclusid+1);
        tempsUnW(clusid,:,:) = sp.tempsUnW(oriclusid+1,:,:);
        templateDuration(clusid) = sp.templateDuration(oriclusid+1);
        waveforms(clusid,:) = sp.waveforms(oriclusid+1,:);
    else % Average them together
        temps(clusid,:,:)=nanmean(sp.temps(oriclusid+1,:,:),1);
        try
            pcFeatInd(clusid,:)=nanmean(sp.pcFeatInd(oriclusid+1,:),1);
        catch ME

        end
        templateDepths(clusid) = nanmean(sp.templateDepths(oriclusid+1),1); % 0-indexed!
        templateXpos(clusid) = nanmean(sp.templateXpos(oriclusid+1),1);
        tempAmps(clusid) = nanmean(sp.tempAmps(oriclusid+1),1);
        tempsUnW(clusid,:,:) = nanmean(sp.tempsUnW(oriclusid+1,:,:),1);
        templateDuration(clusid) = nanmean(sp.templateDuration(oriclusid+1),1);
        waveforms(clusid,:) = nanmean(sp.waveforms(oriclusid+1,:),1);
    end
end
takeclus = 1:nclus;
takeclus(emptyclus)=[];

% Save out
if length(sp.cgs)~=length(takeclus)
    sp.cgs = sp.cgs(takeclus);
    sp.cids = sp.cids(takeclus);
end
sp.temps=temps(takeclus,:,:);
sp.pcFeatInd = pcFeatInd(takeclus,:);
sp.templateDepths = templateDepths(takeclus);
sp.templateXpos = templateXpos(takeclus);
sp.tempAmps = tempAmps(takeclus);
sp.tempsUnW = tempsUnW(takeclus,:,:);
sp.templateDuration = templateDuration(takeclus);
sp.waveforms = waveforms(takeclus,:);
if any(emptyclus)
    disp(['Found ' num2str(length(emptyclus)) ' empty clusters. Removing from clusinfo'])

  
    %% Remove empty clusters from sp
    removeidx = ismember(sp.clu,clusinfo.cluster_id(emptyclus));
    fields = fieldnames(sp);
    for id = 1:length(fields)
        if numel(sp.(fields{id})) == numel(removeidx)
            sp.(fields{id})(removeidx) = [];
        end
    end  
    
    %% Remove the empty clusters
    fields = fieldnames(clusinfo);
    for id = 1:length(fields)
        clusinfo.(fields{id})(emptyclus,:)=[];
    end

else
    disp(['Found no empty clusters.'])
end





















end
