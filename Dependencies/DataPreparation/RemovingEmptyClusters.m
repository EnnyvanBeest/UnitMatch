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
    oriclusid = unique(sp.spikeTemplates(find(sp.clu == clusid-1))); %0-indexed!
    if isempty(oriclusid)
        emptyclus = [emptyclus clusid];
    elseif length(oriclusid)==1
        temps(clusid,:,:)=sp.temps(oriclusid+1,:,:);
        try
            pcFeatInd(clusid,:)=sp.pcFeatInd(oriclusid+1,:);
        catch ME
            
        end
        templateDepths(clusid) = sp.templateDepths(oriclusid+1);
        templateXpos(clusid) = sp.templateXpos(oriclusid+1);
        tempAmps(clusid) = sp.tempAmps(oriclusid+1);
        tempsUnW(clusid,:,:) = sp.tempsUnW(oriclusid+1,:,:);
        templateDuration(clusid) = sp.templateDuration(oriclusid+1);
        waveforms(clusid,:) = sp.waveforms(oriclusid+1,:);
    else
        keyboard
    end
end
takeclus = 1:nclus;
takeclus(emptyclus)=[];

% Save out
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
end

%% Remove the empty clusters
fields = fieldnames(clusinfo);
for id = 1:length(fields)
    eval(['clusinfo.' fields{id} '(emptyclus)=[];'])
    %cluster_id(emptyclus)=[];
%     clusinfo.Amplitude(emptyclus)=[];
%     clusinfo.ContamPct(emptyclus)=[];
%     clusinfo.KSLabel(emptyclus,:)=[];
%     clusinfo.amp(emptyclus)=[];
%     clusinfo.ch(emptyclus)=[];
% 
%     clusinfo.depth(emptyclus)=[];
%     clusinfo.fr(emptyclus)=[];
%     clusinfo.group(emptyclus,:)=[];
%     clusinfo.n_spikes(emptyclus)=[];
%     clusinfo.sh(emptyclus)=[];

end



















end
