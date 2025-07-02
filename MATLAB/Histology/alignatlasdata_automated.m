% Fully Automated alignatlasdata Function
% Uses Bayesian probability estimation with Gaussian smoothing
% Incorporates functional clustering, flexible histology alignment, and hierarchical structure
% Supports multiple shanks with linked transformations

function [Depth2Area, ProbMatrix] = alignatlasdata_automated(histinfo, AllenCCFPath, sp, clusinfo, removenoise, trackcoordinates)
% Load Allen CCF Structure Tree
structureTree = readtable(fullfile(AllenCCFPath, 'structure_tree_safe_2017.csv'));
% correct firing rates where necessary
acronyms = lower(structureTree.acronym);
color_hex = structureTree.color_hex_triplet;

% Ensure valid histinfo input
if ~iscell(histinfo)
    histinfo = {histinfo};
end
if nargin > 5 && ~iscell(trackcoordinates)
    trackcoordinates = {trackcoordinates};
end

% Extract spike and cluster data
spikeTimes = sp.st;
spikeDepths = sp.spikeDepths;
spikeRecSes = sp.RecSes;
amplitudes = sp.spikeAmps;
if ~all(isnan(clusinfo.UniqueID))
    spikeCluster = sp.UniqClu;
    cluster_id = clusinfo.UniqueID;
else
    spikeCluster = sp.clu;
    cluster_id = clusinfo.cluster_id;
end

% we may have multiple shanks, but the current recording may cover just a
% selection of shanks:
ShanksUsed = unique(clusinfo.Shank);
histinfo = histinfo(ShanksUsed+1); % ID + 1 (0-indexed)
numShanks = numel(histinfo);

% identify max spike depth distibution
activityClusters = cell(numShanks, 1);
MUACorr = cell(numShanks, 1);
depthEdges = cell(numShanks, 1);
probeLengths = cell(numShanks,1);
totalProbeLength = cell(numShanks,1);
for shank = 1:numShanks

    % Structural
    % need to flip the probe?
    if trackcoordinates{shank}(end,2)>trackcoordinates{shank}(1,2)
        trackcoordinates{shank} = flipud(trackcoordinates{shank});
        histinfo{shank} = flipud(histinfo{shank});
    end
    probeLengths{shank} = sqrt(sum(trackcoordinates{shank}-trackcoordinates{shank}(end,:),2).^2); % Distance per segment
    totalProbeLength{shank} = max(probeLengths{shank}); % Total depth from the brain surface

    % Functional
    if removenoise
        spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID) & clusinfo.Shank'==ShanksUsed(shank)));
    else
        spikeID = ismember(spikeCluster, cluster_id(clusinfo.Shank'==ShanksUsed(shank)));
    end

    if ~any(spikeID)
        continue
    end


    % Autocorrelation of spikes
    % [autocorr{shank}, depthCenters{shank}] = computeAutocorrByDepth(spikeTimes(spikeID), spikeDepths(spikeID), spikeRecSes(spikeID), 1, 50, height(histinfo));
    % Activity per bin
    [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}, FreqIntens{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});
    % Identify valid depths:
    minDepth = min(spikeDepths)-20; %always the tip
    firingRates{shank} = smoothdata(firingRates{shank},'gaussian',floor(500./nanmedian(unique(diff(depthEdges{shank})))));
    % firstactivity = find(firingRates{shank}>0.05,1,'first');
    % maxDepthtmp1 =  depthEdges{shank}(find(firingRates{shank}(firstactivity:end)>0.05,1,'last')+firstactivity-1);
    % % CrossCorrVar = CrossCorrVar./nanmax(CrossCorrVar);
    % maxDepthtmp2 = depthEdges{shank}(find(CrossCorrVar>0.15,1,'last'));
    % Good unit based
    nGood = smoothdata(histcounts(clusinfo.depth(logical(clusinfo.Good_ID) & clusinfo.Shank'==shank-1),depthEdges{shank}),'gaussian',floor(500./nanmedian(unique(diff(depthEdges{shank})))));

    if all(nGood==0)
        nGood = nGood+1;
    end
    % Ampltidues
    binID = discretize(spikeDepths(spikeID), depthEdges{shank});
    amps = double(amplitudes(spikeID));

    spikeAmps{shank} = accumarray(binID(~isnan(binID)), amps(~isnan(binID)),[numel(depthEdges{shank})-1,1], @mean, NaN);                         
    spikeAmps{shank} = smoothdata(spikeAmps{shank},1,'gaussian',floor(500./nanmedian(unique(diff(depthEdges{shank})))));
    spikeAmps{shank}(isnan(spikeAmps{shank})) = 0;
    spikeAmps{shank} = abs(spikeAmps{shank}-nanmedian( spikeAmps{shank}));


    % Structural
    StructuralPossibility = probeLengths{shank}>=minDepth & probeLengths{shank}<=max(depthEdges{shank}) & depthEdges{shank}(1:end-1)'<totalProbeLength{shank};
    depthscore = (StructuralPossibility'.*(nGood./max(nGood))+0.1).*(firingRates{shank}'./nanmax(firingRates{shank}) + (1-(spikeAmps{shank}'./nanmax(spikeAmps{shank}))) + nanmean(FreqIntens{shank},2)'./nanmax(nanmean(FreqIntens{shank},2)));
    depthscore = depthscore ./ 4;
    
    figure;
    subplot(5,1,1)
    plot(depthEdges{shank}(1:end-1),firingRates{shank}./nanmax(firingRates{shank}));
    title('Firing Rate')

    subplot(5,1,2)
    plot(depthEdges{shank}(1:end-1),nGood./max(nGood));
    title('nGoodClus')

    subplot(5,1,3)   
    plot(depthEdges{shank}(1:end-1),1-(spikeAmps{shank}./nanmax(spikeAmps{shank})))
    title(' Average Amplitude')

    subplot(5,1,4)
    plot(depthEdges{shank}(1:end-1),FreqIntens{shank}./nanmax(FreqIntens{shank},[],1))
    title('Power at Frequencies')

    subplot(5,1,5)
    plot(depthEdges{shank}(1:end-1),depthscore);
    title('Depth Score')
    hold on
    line(get(gca,'xlim'),[0.05 0.05],'color',[1 0 0])

    maxDepthtmp = depthEdges{shank}(find(depthscore>0.05,1,'last')+1);

    maxDepth{shank} = maxDepthtmp;
       
end
if ~exist('maxDepth')
    Depth2Area = [];
    return
end
maxDepth = nanmax(cat(1,maxDepth{:}));
DepthRange = maxDepth - minDepth;

spikeTimes(spikeDepths<minDepth | spikeDepths>maxDepth) = [];
spikeCluster(spikeDepths<minDepth | spikeDepths>maxDepth) = [];
spikeDepths(spikeDepths<minDepth | spikeDepths>maxDepth) = [];


% Clean up histinfo a little bit
for shank = 1:numShanks
    if totalProbeLength{shank}>DepthRange
        removeId = probeLengths{shank}>maxDepth | probeLengths{shank}<minDepth;  %  Remove these entries
        histinfo{shank}(removeId,:) = [];%-
        trackcoordinates{shank}(removeId,:) = [];
        probeLengths{shank}(removeId) = [];
        depthEdges{shank}(removeId) = [];
        spikeAmps{shank}(removeId) = [];

        % Need to rerun
        if removenoise
            spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID) & clusinfo.Shank'==ShanksUsed(shank)));
        else
            spikeID = ismember(spikeCluster, cluster_id(clusinfo.Shank'==ShanksUsed(shank)));
        end
        [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}, FreqIntens{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});

    end

    [histinfo{shank} removeid] = cleanHistinfo(histinfo{shank});
    if ~isempty(removeid)
        trackcoordinates{shank}(removeid,:) = [];
        probeLengths{shank}(removeid) = [];
        depthEdges{shank}(removeid) = [];
        spikeAmps{shank}(removeid) = [];

        % Need to rerun
        if removenoise
            spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID) & clusinfo.Shank'==ShanksUsed(shank)));
        else
            spikeID = ismember(spikeCluster, cluster_id(clusinfo.Shank'==ShanksUsed(shank)));
        end
        [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}, FreqIntens{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});

    end
    if any(depthEdges{shank}>maxDepth | depthEdges{shank}<minDepth) % We still need to cut our depthEdges down

        % Need to rerun
        if removenoise
            spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID) & clusinfo.Shank'==ShanksUsed(shank)));
        else
            spikeID = ismember(spikeCluster, cluster_id(clusinfo.Shank'==ShanksUsed(shank)));
        end
        [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}, FreqIntens{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});
    end
end



% Align histology data with functional clusters using the same depth bins
alignedHistology = cell(numShanks, 1);
for shank = 1:numShanks
    histDepths = depthEdges{shank}(1:end-1) - diff(depthEdges{shank})/2;
    [alignedHistology{shank}, Features{shank}] = stretchHistology(histinfo{shank}, activityClusters{shank}, MUACorr{shank}, FreqIntens{shank}, spikeAmps{shank}, histDepths, firingRates{shank}, trackcoordinates{shank}, structureTree);
end


% Generate deterministic depth-to-area mapping
Depth2Area = cell(numShanks, 1);
for shank = 1:numShanks
    Depth2Area{shank} = assignAreas(alignedHistology{shank}, acronyms, shank, color_hex);
end

% Visualize results
for shank = 1:numShanks
    visualizeResults(Depth2Area{shank}, Features{shank});
end
Depth2Area=cat(1,Depth2Area{:});

end

function [histinfo, removeidx] = cleanHistinfo(histinfo)
removeidx = find(ismember(histinfo.RegionAcronym,'Not found in brain'));
histinfo(removeidx,:) = [];
[uniqueAreas,id1,id2] = unique(histinfo.RegionAcronym,'stable');

% Check variable names
if ~any(ismember(histinfo.Properties.VariableNames,'Position')) & any(ismember(histinfo.Properties.VariableNames,'Index'))
    histinfo.Position = histinfo.Index;
end
for areaid = 1:numel(uniqueAreas)
    idx = find(id2==areaid);
    if numel(idx)==1
        id2(idx) = nan;
    end
    difvec = [1; diff(idx)];
    weirdidx = find(difvec>1);
    if idx(end) == height(histinfo)
        difvec = [difvec; 2];
    else
        difvec = [difvec; nan];
    end
    for idd = weirdidx'
        if (difvec(idd)>1 & difvec(idd-1)>1) || (difvec(idd)>1 & difvec(idd+1)>1)
            id2(idx(idd)) = nan;
        end
    end
end
id2 = round(fillmissing(id2,'linear'));
histinfo.RegionAcronym = uniqueAreas(id2);
end

function [activityClusters, MUACorr, firingRates, depthEdges, FreqIntens] = computeFunctionalClusters(spikeTimes, spikeDepths, histinfo)
nDepthBins = height(histinfo);
depthEdges = linspace(min(spikeDepths), max(spikeDepths), nDepthBins + 1);
timeEdges = min(spikeTimes):1:max(spikeTimes);
nSmooth = ceil(40./nanmedian(diff(depthEdges))); % smooth for MUA corr with 50 micron
spikeHist = histcounts2(spikeDepths, spikeTimes, depthEdges, timeEdges);

% Detect noise - and correct
noiselevel = quantile(spikeHist(spikeHist>0),0.99);
spikeHist(:,sum(spikeHist>noiselevel,1)>5) = [];
noiselevel = quantile(spikeHist(spikeHist>0),0.99);
spikeHist(spikeHist>noiselevel) = nan;
spikeHist = smoothdata(spikeHist,1,'gaussian',nSmooth);
MUACorr = corr(spikeHist');

% MUACorr = nan(numel(depthEdges)-1,numel(depthEdges)-1,5);
% for did1 = 1:numel(depthEdges)-1
%     for did2 = 1:numel(depthEdges)-1
%         MUACorr(did1,did2,:) = xcorr(spikeHist(did1,:)',spikeHist(did2,:)',2);
%     end
% end

nAreas = numel(unique(histinfo.RegionAcronym));
spikeHist(isnan(spikeHist)) = 0;
[W, ~] = nnmf(spikeHist, nAreas);
activityClusters = W;
activityClusters = activityClusters./nanmax(activityClusters,[],1);
activityClusters(:,isnan(nanmax(activityClusters,[],1))) = [];
firingRates = mean(spikeHist, 2);
activityClusters = cat(2,activityClusters,firingRates/ max(mean(spikeHist, 2)));

% Spectogram
Fs         = 500;           % desired “sampling rate” (Hz) for your rate trace
dt         = 1/Fs;
tMax       = max(spikeTimes);
edges      = 0:dt:tMax;      % bin‐edges for histogram
newDepthEdges = min(spikeDepths):50:max(spikeDepths);

%--- 2) Bin spikes into a spike‐train
spikeTrain = histcounts2(spikeDepths, spikeTimes, newDepthEdges, edges);

%--- 3) (Optional) Smooth the binned counts to get an instantaneous rate
Bands = [0.1 0.5; 0.5 4; 4 8; 8 12; 13 30; 30 60; 60 100; 100 200];
FreqIntens = nan(numel(newDepthEdges)-1,size(Bands,1));
disp('Computing frequency band strength, can take a few minutes...')
for bid = 1:size(Bands,1)
    if size(spikeTrain,2)>Fs*7*60 % Take just 5 mins of data
        sampleidx = Fs*2*60:Fs*7*60;
    else
        sampleidx = 1:size(spikeTrain,2);
    end
    FreqIntens(:,bid) = bandpower(spikeTrain(:,sampleidx)', Fs, Bands(bid,:));
end

% 1) compute bin‐centers
oldZ = (newDepthEdges(1:end-1) + newDepthEdges(2:end))/2;    % 1×M
newZ = (depthEdges(1:end-1) + depthEdges(2:end))/2; % 1×K

% 2) interpolate
%    FreqIntens is M×N, so interp1 will return K×N
FreqIntens = interp1(oldZ, FreqIntens, newZ, 'linear', 'extrap' );

end

function [autocorrMatrix, depthCenters] = computeAutocorrByDepth(spikeTimes, spikeDepths, spikeRecSes, binSizeMs, maxLagMs, depthBinSize)
% Compute and plot spike autocorrelation across probe depth bins
%
% Inputs:
%   spikeTimes     - vector of spike times (in seconds)
%   spikeDepths    - vector of spike depths (in microns)
%   binSizeMs      - temporal bin size for autocorrelation (ms)
%   maxLagMs       - max lag for autocorrelation (ms)
%   depthBinSize   - depth bin size (microns)

% Set up depth bins
minDepth = floor(min(spikeDepths));
maxDepth = ceil(max(spikeDepths));
depthEdges = minDepth:depthBinSize:maxDepth;
depthCenters = depthEdges(1:end-1) + depthBinSize/2;

% Preallocate
maxLagBins = round(maxLagMs / binSizeMs);
lagRange = -maxLagBins:maxLagBins;
nLags = length(lagRange);
nBins = length(depthCenters);
RecSesOpt = unique(spikeRecSes);
autocorrMatrix = nan(nBins, nLags);

% Loop through depth bins
for recid = 1:numel(RecSesOpt)
    for i = 1:nBins
        depthRange = [depthEdges(i), depthEdges(i+1)];

        % Filter spike times
        inRange = spikeDepths >= depthRange(1) & spikeDepths < depthRange(2) & spikeRecSes == RecSesOpt(recid);
        times = spikeTimes(inRange);

        if numel(times) < 50
            continue; % Skip bins with too few spikes
        end

        % Bin spikes
        binSizeSec = binSizeMs / 1000;
        minTime = min(times);
        maxTime = max(times);
        edges = minTime:binSizeSec:maxTime;
        binnedSpikes = histcounts(times, edges);

        % Compute autocorrelation
        [ac, ~] = xcorr(binnedSpikes, maxLagBins, 'coeff');

        % Remove zero-lag if desired
        ac(maxLagBins+1) = 0;

        % Store
        autocorrMatrix(i, :) = nanmean(cat(1,autocorrMatrix(i, :),ac),1);
    end
end

% Plot
figure;
imagesc(lagRange * binSizeMs, depthCenters, autocorrMatrix);
xlabel('Lag (ms)');
ylabel('Depth (µm)');
title('Spike Autocorrelation Across Depth');
colormap turbo;
colorbar;
set(gca, 'YDir', 'normal'); % depth upwards
end

function [alignedHistology, F] = stretchHistology(histinfo, activityClusters, MUACorr, FreqIntens, spikeAmps, histDepths, firingRates, trackcoordinates, structureTree)

% Normalize all scores:
nDepth = size(MUACorr,1);
FreqIntens = (FreqIntens-nanmin(FreqIntens,[],1))./(nanmax(FreqIntens,[],1)-nanmin(FreqIntens,[],1));
spikeAmps = (spikeAmps-nanmin(spikeAmps))./(nanmax(spikeAmps)-nanmin(spikeAmps));
firingRates =  (firingRates-nanmin(firingRates))./(nanmax(firingRates)-nanmin(firingRates));
% pull out the PCA scores
scores = nan(nDepth,size(MUACorr,3));
for tl = 1:size(MUACorr,3)
    C0 = squeeze(MUACorr(:,:,tl));
    C0=fillmissing(C0,'constant',0,1);
    C0=fillmissing(C0,'constant',0,2);

    % each row i is the full‐depth correlation profile of depth i
    [coeff, scorestmp, ~, ~, explained] = pca(C0, 'NumComponents', 1);

    scores(:,tl) = scorestmp;
 
end
scores = (scores-nanmin(scores,[],1))./(nanmax(scores,[],1)-nanmin(scores,[],1));

% concatenate
F = [FreqIntens, spikeAmps, firingRates, activityClusters, scores];
F(:,sum(isnan(F),1)==size(F,1)) = [];

% PCA → first PC
pcnum = 0;
explained = 0;
while sum(explained(1:pcnum))<65
    pcnum = pcnum + 1;

    [~,PC1,~,~,explained] = pca(double(F),'NumComponents',pcnum);
    fprintf('PCA: %.0f PC explain %.1f%% of variance\n', pcnum, sum(explained(1:pcnum)));
end
PC1 = sum(PC1,2);
PC1 = detrend(PC1);
PC1 = (PC1-nanmin(PC1))./(nanmax(PC1)-nanmin(PC1));

% Unique areas and gray-like identification
[uniqueAreas, ~, table2areaidx] = unique(histinfo.RegionAcronym, 'stable');

% Identify a switch:
SwitchIdx = [1; find(abs(diff(table2areaidx))>0)+1];
uniqueAreas = uniqueAreas(table2areaidx(SwitchIdx));

% colours and expected firing rates
grayscaleThreshold = 0.1; % Threshold for gray-like areas
grayLike = zeros(1, numel(uniqueAreas));
EFR = zeros(1, numel(uniqueAreas));
for areaid = 1:numel(uniqueAreas)
    areaColors = structureTree.color_hex_triplet(ismember(lower(structureTree.acronym), lower(uniqueAreas(areaid))));
    EFR(areaid) = structureTree.expected_firing_rate(ismember(lower(structureTree.acronym), lower(uniqueAreas(areaid))));
    grayLike(areaid) = cellfun(@(c) nanmean(abs(diff(hex2rgb(c)))) < grayscaleThreshold, areaColors);
end

EFR(logical(grayLike)) = 0; % we still want gray like areas to be squeezed

% Get parent area
parentArea = cell(1,numel(uniqueAreas));
for areaid = 1:numel(uniqueAreas)
    parentArea{areaid} = getParentArea(uniqueAreas{areaid}, structureTree);
end
parentArea(logical(grayLike)) = {'WhiteMatter'}

% Identify expected depth per parent area
[uniqueParentArea, ~, idToParent] = unique(parentArea,'stable');
% Identify a switch:
SwitchIdx = [1; find(abs(diff(idToParent))>0)+1];
uniqueParentArea = uniqueParentArea(idToParent(SwitchIdx));

numPerParent = nan(1,numel(uniqueParentArea));
EFRPerParent = nan(1,numel(uniqueParentArea));
alreadyused = false(1,numel(parentArea));
startidx = 1;
for parid = 1:numel(uniqueParentArea)
    if parid < numel(uniqueParentArea)
        NextAreaidx = find(ismember(parentArea,uniqueParentArea(parid+1)) & ~alreadyused,1,'first');%identify first occurance of this area
        stopidx = find(ismember(histinfo.RegionAcronym(startidx:end),uniqueAreas(ismember(parentArea,parentArea(NextAreaidx)))),1,'first')+startidx-1;
    else
        NextAreaidx = numel(parentArea);%identify first occurance of this area
        stopidx = height(histinfo);
    end
    ThisAreaidx = find(ismember(parentArea,uniqueParentArea(parid)) & ~alreadyused,1,'first');%identify first occurance of this area

    if ~grayLike(ThisAreaidx)
        numPerParent(parid) = numel(histinfo.Position(find(ismember(histinfo.RegionAcronym(startidx:stopidx),uniqueAreas(ismember(parentArea,parentArea(ThisAreaidx)))))+startidx-1));
        EFRPerParent(parid) = nanmean(EFR(ismember(parentArea,parentArea(ThisAreaidx))));
        alreadyused(1:NextAreaidx-1) = true;
    else
        numPerParent(parid) = 0; % Artifically set gray areas to 0 coverage
        EFRPerParent(parid) = nanmean(EFR(ismember(parentArea,parentArea(ThisAreaidx))));
    end
  

    startidx = stopidx;
end

% if a parent area occured more than once, divide numPerParent accordingly
[~,id1,id2] = unique(uniqueParentArea);
for id = 1:numel(id1)
    id3 = id2 == id2(id);
    numPerParent(id3) = (nansum(numPerParent(id3)))./sum(id3);
end
% Convert this into distance
CoveragePerParent = numPerParent.*(max(histDepths)-min(histDepths))/height(histinfo);

% 1) compute absolute first difference
d = (diff(smoothdata(PC1,'gaussian',50)));       % (N-1)×1

% 2) sort descending, pick top nChanges
topIdx = [];
qval = 1;
while numel(topIdx)<numel(CoveragePerParent) & qval>0.75
    qval = qval-0.05;
    [pks,topIdx] = findpeaks(d,'MinPeakDistance',5,'MinPeakHeight',quantile(d,qval),'NPeaks',numel(CoveragePerParent));
end
changePts = sort(topIdx,'ascend') + 1;

% first, identify all the distances between two 'low firing depths'
starterid = 1;
stopid = 1;

CoverageBetweenLowFRs = [];
SaveStarterID = starterid;
SaveStopID = [];
while stopid<numel(histDepths)
    stopid = changePts(find(changePts > starterid, 1, 'first'));
    if isempty(stopid)
        break
    end
    SaveStopID = [SaveStopID stopid];
    if isempty(stopid), stopid = numel(histDepths); end
    CoverageBetweenLowFRs = [CoverageBetweenLowFRs histDepths(stopid) - histDepths(starterid)];
    starterid = stopid+1;
    if starterid==numel(histDepths)
        break
    end
    SaveStarterID = [SaveStarterID starterid];
end
SaveStopID = [SaveStopID numel(histDepths)];



CoverageBetweenLowFRs = CoverageBetweenLowFRs./(sum(CoverageBetweenLowFRs)./sum(CoveragePerParent));
% Compute cumulative sums for optimal matching
cumulativeHistology = cumsum((CoveragePerParent),'omitnan');
cumulativeFunctional = cumsum((CoverageBetweenLowFRs));


% Find best overal alignment using least-squares approach
[~, assignment] = min(abs(cumulativeHistology' - cumulativeFunctional), [], 2);
assignment(EFRPerParent==0) = nan;
nulasignment = find(~ismember(1:numel(cumulativeFunctional),assignment));
counter = 0;
for id = 1:numel(nulasignment)
    assignment(assignment>nulasignment(id)-counter) = assignment(assignment>nulasignment(id)-counter)-1; % Correct index in assignment
    SaveStarterID(nulasignment(id)-counter+1) = SaveStarterID(nulasignment(id)-counter); % broaden the start and stop
    SaveStarterID(nulasignment(id)-counter) = []; % Remove the starter id for this assignment
    SaveStopID(nulasignment(id)-counter) = []; % Remove the stop id for this assignment
    counter = counter + 1;
end
SaveStarterID = SaveStarterID(1:max(assignment));

SaveStopID = SaveStarterID(2:end)-1;
SaveStopID = [SaveStopID numel(histDepths)];
% Iterative refinement of depth realignment
newDepths = nan(1, numel(histDepths));
newTrackcoordinates = nan(size(trackcoordinates));
for id = 1:max(assignment)
    starterid = SaveStarterID(id);
    stopid = SaveStopID(id);
    histdepthtmp = histDepths(starterid:stopid);
    trackcoordinatestmp = trackcoordinates(starterid:stopid, :);

    % Areas to consider in this assignment round
    TheseAreas = unique(uniqueParentArea(assignment == id & EFRPerParent'~=0),'stable');
    IncludeThis = false(1,numel(TheseAreas));
    % Are all these areas represented?
    for areaid = 1:numel(TheseAreas)
        if ismember(uniqueAreas(ismember(parentArea,TheseAreas{areaid})),histinfo.RegionAcronym(starterid:stopid))
            IncludeThis(areaid) = true;
        end
    end
    TheseAreas = TheseAreas(IncludeThis);

    nAreas = numel(TheseAreas);

    % Where to draw the line between these areas?
    tmpVec = starterid:stopid;
    if nAreas>1
        d = diff(PC1(starterid:stopid));
        topIdx = [];
        qval = 1;
        while numel(topIdx)<nAreas-1
            qval = qval-0.05;
            if qval<0.5
                % Use ratio in histology instead
                [~,id1] = ismember(histinfo.RegionAcronym(starterid:stopid),uniqueAreas(ismember(parentArea,TheseAreas)));
                nPerRegion = arrayfun(@(X) sum(id1==X),1:nAreas);
                nPerRegion = nPerRegion./sum(nPerRegion);
                ratioPerRegion = cumsum(nPerRegion);
                transitionPoints = [1 round(ratioPerRegion.*numel(tmpVec))];
                transitionPoints(transitionPoints==0) = [];
                break
            end
            [pks,topIdx] = findpeaks(d,'MinPeakDistance',2,'MinPeakHeight',quantile(d,qval),'NPeaks',nAreas-1);
            changePts = sort(topIdx,'ascend') + 1;
            transitionPoints = [1 changePts' numel(tmpVec)];
        end
     

    else
        transitionPoints = [1 numel(tmpVec)];
    end


    % Interpolation
    for aid2 = 1:nAreas
        tblidx = find(ismember(histinfo.RegionAcronym, uniqueAreas(ismember(parentArea,TheseAreas{aid2}))));
        tblidx(tblidx<starterid) = [];
        tblidx(tblidx>stopid) = [];
        newDepths(tblidx) = linspace(histdepthtmp(transitionPoints(aid2)), histdepthtmp(transitionPoints(aid2+1)), numel(tblidx));
        if any(newDepths(tblidx)' - newDepths(1:starterid-1)<0,'all')
            keyboard
        end
        for dimid = 1:3
            newTrackcoordinates(tblidx, dimid) = linspace(trackcoordinatestmp(transitionPoints(aid2), dimid), trackcoordinatestmp(transitionPoints(aid2+1), dimid), numel(tblidx));
        end
    end
end
PosOrNeg = sign(nanmean(diff(histDepths)));

if isnan(newDepths(1))
    if PosOrNeg>0
        newDepths(1) = min(histDepths);
    else
        newDepths(1) = max(histDepths);
    end
end
if isnan(newDepths(end))
    if PosOrNeg>0
        newDepths(end) = max(histDepths);
    else
        newDepths(end) = min(histDepths);
    end
end
newDepths = fillmissing(newDepths,'linear','EndValues','none');
newDepths(newDepths<0) = 0; % Remove extremes
newDepths(newDepths>max(histDepths)) = max(histDepths);

for dimid = 1:size(newTrackcoordinates,2)
    PosOrNeg = sign(nanmean(diff(trackcoordinates(:,dimid))));
    if isnan( newTrackcoordinates(1,dimid))
        if PosOrNeg>0
            newTrackcoordinates(1,dimid) = min(trackcoordinates(:,dimid));
        else
            newTrackcoordinates(1,dimid) = max(trackcoordinates(:,dimid));
        end
    end
    if isnan( newTrackcoordinates(end,dimid))
        if PosOrNeg>0
            newTrackcoordinates(end,dimid) = max(trackcoordinates(:,dimid));
        else
            newTrackcoordinates(end,dimid) = min(trackcoordinates(:,dimid));
        end
    end
    newTrackcoordinates(:,dimid) = fillmissing(newTrackcoordinates(:,dimid),'linear','EndValues','none');
end
% [newDepths,sortidx] = sort(newDepths,'ascend');
histinfo.Depth = newDepths';
histinfo.trackcoordinates = newTrackcoordinates;
alignedHistology = histinfo;
end




function parentArea = getParentArea(area, structureTree)
idx = find(strcmpi(structureTree.acronym, area));
if isempty(idx)
    parentArea = area; % If area not found, return itself
else
    structurePath = structureTree.structure_id_path{idx};
    pathParts = strsplit(structurePath, '/');
    if numel(pathParts) > 3
        parentID = str2double(pathParts{end-2});
        parentIdx = find(structureTree.id == parentID);
        parentArea = structureTree.acronym{parentIdx};
    else
        parentArea = area; % If no parent, return itself
    end
end
end

function Depth2Area = assignAreas(alignedHistology, acronyms, shank, color_hex)
% Assigns brain areas to depths based on aligned histology
depths = alignedHistology.Depth; % Compute bin centers
colPerDepth = cell(1,numel(alignedHistology.RegionAcronym));
for aid = 1:numel(alignedHistology.RegionAcronym)

    colPerDepth(aid) = color_hex(find(ismember(acronyms,lower(alignedHistology.RegionAcronym(aid)))));

end

% Ensure correct mapping between depths and areas
if numel(depths) ~= numel(alignedHistology.RegionAcronym)
    error('Mismatch between depth bins and area labels. Check alignment.');
end

% Create table mapping each depth to an area
Depth2Area = table(depths, repmat(shank, numel(depths),1), alignedHistology.RegionAcronym, colPerDepth', alignedHistology.trackcoordinates, 'VariableNames', {'Depth','Shank','Area','Color','Coordinates'});

end


function visualizeResults(Depth2Area, Features)
% Order areas by depth
[sortedDepths, orderIdx] = sort(Depth2Area.Depth,'descend');
sortedAreas = Depth2Area.Area(orderIdx);
sortedCols = Depth2Area.Color(orderIdx);

% Convert area names to colors
[uniqueAreas, idx] = unique(sortedAreas, 'stable');
colorMap = containers.Map(uniqueAreas, sortedCols(idx));
areaColors = cellfun(@(x) colorMap(x), sortedAreas, 'UniformOutput', false);

% Convert hex colors to RGB
areaColorsRGB = cellfun(@(c) sscanf(c, '%2x%2x%2x')' / 255, areaColors, 'UniformOutput', false);
areaColorsRGB = vertcat(areaColorsRGB{:});

figure('Position',[22 14 560 900]);
h(1) = subplot(1,3,1);
scatter(categorical(sortedAreas,unique(sortedAreas, 'stable')), sortedDepths, 25, areaColorsRGB, 'filled');
set(gca,'YDir','normal')
title('Ordered Depth-to-Area Mapping');
ylabel('Depth (µm)');
xlabel('Brain Area');
ylim([min(sortedDepths) max(sortedDepths)])
makepretty

h(2) = subplot(1,3,[2,3]);
imagesc([],sortedDepths,Features(orderIdx, :));
set(gca,'ydir','normal')
colormap(hot)
colorbar;
title('Features');
xlabel('Feature');
ylabel('Depth');
ylim([min(sortedDepths) max(sortedDepths)])
makepretty

linkaxes(h,'y')
end
