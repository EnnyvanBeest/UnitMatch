function spikeStruct = loadKSdir(ksDir, varargin)

if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = false;
end
if ~isfield(params, 'loadPCs')
    params.loadPCs = false;
end

% load spike data
spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

if exist(fullfile(ksDir, 'spike_times_corrected.npy'))
    ss = readNPY(fullfile(ksDir, 'spike_times_corrected.npy')); %For chronic sessions
    sessionid = ones(size(ss));
    idx = find([0; diff(ss)]<0);
    for id = 1:length(idx)
        sessionid(idx(id):end) = id+1;
    end
else
    ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
    sessionid = ones(size(ss));

end
st = double(ss)/spikeStruct.sample_rate;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed

if exist(fullfile(ksDir,'spike_datasets.npy'))
    datas = readNPY(fullfile(ksDir,'spike_datasets.npy'));
else
    datas = zeros(length(spikeTemplates),1);
    params.thisdate = [];
end

if ~isfield(params,'thisdate') || isempty(params.thisdate)
    datasetidx = unique(datas);
else 
    datasetidx = find(cell2mat(cellfun(@(X) any(strfind(X,params.thisdate)),strsplit(spikeStruct.dat_path,','),'UniformOutput',0)))-1; %Which dataset to take?
end

if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy')); %Changed by phy
else
    clu = spikeTemplates; % Original
end


tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv')) 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 


if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        datas = datas(~ismember(clu, noiseClusters));
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
        
        
    end
    
else
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
   
% Only take needed data
spikeTemplates = spikeTemplates(ismember(datas,datasetidx));
tempScalingAmps = tempScalingAmps(ismember(datas,datasetidx));
st = st(ismember(datas,datasetidx));
clu = clu(ismember(datas,datasetidx));
datas = datas(ismember(datas,datasetidx));
if params.loadPCs
    pcFeat = pcFeat(ismember(datas,datasetidx), :,:);
    %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
end
        

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

spikeStruct.st = st;
spikeStruct.SessionID = sessionid;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;
spikeStruct.dataset = datas;