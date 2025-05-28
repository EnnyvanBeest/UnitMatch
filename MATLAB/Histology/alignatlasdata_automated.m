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
spikeCluster = sp.clu;
spikeTimes = sp.st;
spikeDepths = sp.spikeDepths;
cluster_id = clusinfo.cluster_id;
numShanks = numel(histinfo);

% Optional: Remove noise clusters
if removenoise
    spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID)));
else
    spikeID = true(size(spikeCluster));
end
% identify max spike depth distibution
activityClusters = cell(numShanks, 1);
MUACorr = cell(numShanks, 1);
depthEdges = cell(numShanks, 1);
for shank = 1:numShanks

    % Autocorrelation of spikes
    [autocorr{shank}, depthCenters{shank}] = computeAutocorrByDepth(spikeTimes, spikeDepths, 1, 50, 50);
    % Activity per bin
    [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});
    % Identify valid depths:
    minDepth = min(spikeDepths); %always the tip
    firingRates{shank} = smoothdata(firingRates{shank},'gaussian',floor(500./nanmedian(unique(diff(depthEdges{shank})))));
    firstactivity = find(firingRates{shank}>0.1,1,'first');
    maxDepthtmp1 =  depthEdges{shank}(find(firingRates{shank}(firstactivity:end)>1,1,'last')+firstactivity-1);
    autoCorrVar = smoothdata(nanvar(autocorr{shank},[],2),'gaussian',floor(500./nanmedian(unique(diff(depthCenters{shank})))));
    autoCorrVar = autoCorrVar./nanmax(autoCorrVar);
    maxDepthtmp2 = depthCenters{shank}(find(autoCorrVar>0.1,1,'last'));
    maxDepthtmp = min([maxDepthtmp1,maxDepthtmp2]);
    if isempty(maxDepthtmp)
        maxDepthtmp = depthEdges{shank}(end)./2;
    end
    maxDepth{shank} = maxDepthtmp;
end
maxDepth = nanmean(cat(1,maxDepth{:}));
DepthRange = maxDepth - minDepth;

% Compute functional clustering (NMF, firing rates, correlation) for each shank
% Adjust histinfo
probeLengths = cell(numShanks,1);
totalProbeLength = cell(numShanks,1);
for shank = 1:numShanks
    % need to flip the probe?
    if trackcoordinates{shank}(end,2)>trackcoordinates{shank}(1,2)
        trackcoordinates{shank} = flipud(trackcoordinates{shank});
        histinfo{shank} = flipud(histinfo{shank});
    end
    probeLengths{shank} = sqrt(sum(trackcoordinates{shank}-trackcoordinates{shank}(end,:),2).^2); % Distance per segment
    totalProbeLength{shank} = max(probeLengths{shank}); % Total depth from the brain surface

    if totalProbeLength{shank}>DepthRange
        removeId = probeLengths{shank}>maxDepth | probeLengths{shank}<minDepth;  %  Remove these entries 
        histinfo{shank}(removeId,:) = [];%- 
        trackcoordinates{shank}(removeId,:) = [];
        probeLengths{shank}(removeId) = [];
    end

end
spikeTimes(spikeDepths<minDepth | spikeDepths>maxDepth) = [];
spikeCluster(spikeDepths<minDepth | spikeDepths>maxDepth) = [];
spikeDepths(spikeDepths<minDepth | spikeDepths>maxDepth) = [];


% Optional: Remove noise clusters
if removenoise
    spikeID = ismember(spikeCluster, cluster_id(logical(clusinfo.Good_ID)));
else
    spikeID = true(size(spikeCluster));
end

% Clean up histinfo a little bit
for shank = 1:numShanks
    histinfo{shank} = cleanHistinfo(histinfo{shank});
    trackcoordinates{shank} = trackcoordinates{shank}(1:size(histinfo{shank},1),:);

end

% Redo
activityClusters = cell(numShanks, 1);
MUACorr = cell(numShanks, 1);
depthEdges = cell(numShanks, 1);
for shank = 1:numShanks
    [activityClusters{shank}, MUACorr{shank}, firingRates{shank}, depthEdges{shank}] = computeFunctionalClusters(spikeTimes(spikeID), spikeDepths(spikeID), histinfo{shank});
end



% Align histology data with functional clusters using the same depth bins
alignedHistology = cell(numShanks, 1);
for shank = 1:numShanks
    histDepths = depthEdges{shank}(1:end-1) - diff(depthEdges{shank})/2;
    alignedHistology{shank} = stretchHistology(histinfo{shank}, activityClusters{shank}, MUACorr{shank}, histDepths, firingRates{shank}, trackcoordinates{shank}, structureTree)
end


% Generate deterministic depth-to-area mapping
Depth2Area = cell(numShanks, 1);
for shank = 1:numShanks
    Depth2Area{shank} = assignAreas(alignedHistology{shank}, acronyms, shank, color_hex);
end

% Visualize results
for shank = 1:numShanks
    visualizeResults(Depth2Area{shank}, MUACorr{shank}, activityClusters{shank});
end
Depth2Area=cat(1,Depth2Area{:});

end

function histinfo = cleanHistinfo(histinfo)
histinfo(ismember(histinfo.RegionAcronym,'Not found in brain'),:) = [];
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

function [activityClusters, MUACorr, firingRates, depthEdges] = computeFunctionalClusters(spikeTimes, spikeDepths, histinfo)
nDepthBins = height(histinfo);
depthEdges = linspace(min(spikeDepths), max(spikeDepths), nDepthBins + 1);
timeEdges = min(spikeTimes):1:max(spikeTimes);
nSmooth = ceil(40./nanmedian(diff(depthEdges))); % smooth for MUA corr with 50 micron
spikeHist = histcounts2(spikeDepths, spikeTimes, depthEdges, timeEdges);
spikeHist = smoothdata(spikeHist,1,'gaussian',nSmooth);
MUACorr = corr(spikeHist');

nAreas = numel(unique(histinfo.RegionAcronym));
[W, ~] = nnmf(spikeHist, nAreas);
activityClusters = W;
activityClusters = activityClusters./nanmax(activityClusters,[],1);
activityClusters(:,isnan(nanmax(activityClusters,[],1))) = [];
firingRates = mean(spikeHist, 2);
activityClusters = cat(2,activityClusters,firingRates/ max(mean(spikeHist, 2)));
end

function [autocorrMatrix, depthCenters] = computeAutocorrByDepth(spikeTimes, spikeDepths, binSizeMs, maxLagMs, depthBinSize)
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
autocorrMatrix = nan(nBins, nLags);

% Loop through depth bins
for i = 1:nBins
    depthRange = [depthEdges(i), depthEdges(i+1)];
    
    % Filter spike times
    inRange = spikeDepths >= depthRange(1) & spikeDepths < depthRange(2);
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
    autocorrMatrix(i, :) = ac;
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

function alignedHistology = stretchHistology(histinfo, activityClusters, MUACorr, histDepths, firingRates, trackcoordinates, structureTree)
% Compute functional area transitions using both NMF and MUA correlations
gradientActivity = cat(1, zeros(1, size(activityClusters, 2)), abs(diff(activityClusters, 1, 1))); % Gradient across NMF components
% transitionScoreNMF = var(gradientActivity, 0, 2); % Variance across NMF components
transitionScoreNMF = [0; sum(diff(gradientActivity,1).^2,2)];
MUACorr(isnan(MUACorr)) = 0;

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
% Convert this into distance
CoveragePerParent = numPerParent.*(max(histDepths)-min(histDepths))/height(histinfo);

% Identify low firing depths
% smoothFiringRates = smoothdata(firingRates, 'gaussian', 4);
lowFiring = firingRates < 1.5;

% Adaptive transition clustering
numClusters = numel(uniqueParentArea);
transitionScoreMUA = [0; abs(diff(kmeans(MUACorr, numClusters, 'Replicates', 5)))];
transitionScoreMUA = double(imbinarize(transitionScoreMUA));
transitionScore = smoothdata(transitionScoreNMF ./ nanmax(transitionScoreNMF) + transitionScoreMUA ./ nanmax(transitionScoreMUA), 'gaussian', 1);
transitionScore(isnan(transitionScore)) = 0;
transitionScore(~isfinite(transitionScore)) = max(transitionScore(isfinite(transitionScore)), [], 'omitnan');

PeakHeight =  quantile(transitionScore, 1 - numel(uniqueParentArea) ./ numel(transitionScore));
[~,  parenttransitionpoints] = findpeaks(transitionScore, 'MinPeakHeight',PeakHeight, 'SortStr', 'descend');

Borders = transitionScore>PeakHeight | lowFiring; % Borders are likely to suddenly have low firing or transition the functional dynamics

% first, identify all the distances between two 'low firing depths'
starterid = 1;
stopid = 1;

CoverageBetweenLowFRs = [];
SaveStarterID = starterid;
SaveStopID = [];
PeakVal = quantile(Borders,1-sum(isnan(CoveragePerParent))./numel(Borders));
while stopid<numel(histDepths)
    stopid = find(Borders(starterid+2:end) >= PeakVal, 1, 'first') + starterid + 1;
    if isempty(stopid)
        break
    end
    SaveStopID = [SaveStopID stopid];
    if isempty(stopid), stopid = numel(histDepths); end
    CoverageBetweenLowFRs = [CoverageBetweenLowFRs histDepths(stopid) - histDepths(starterid)];
    starterid = find(Borders(stopid:end) < PeakVal, 1, 'first') + stopid - 2;
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
    TheseAreas = uniqueParentArea(assignment == id & EFRPerParent'~=0);
    nAreas = numel(TheseAreas);

    transitionPoints = parenttransitionpoints(parenttransitionpoints>starterid&parenttransitionpoints<stopid);
    if ~isempty(transitionPoints) & numel(transitionPoints)>=nAreas-1
        transitionPoints = [starterid,sort(transitionPoints(1:nAreas-1)','ascend'),stopid];
    else
        transitionPoints = round(linspace(starterid,stopid,nAreas+1));
    end
    transitionPoints = transitionPoints-starterid+1;

    % Interpolation
    for aid2 = 1:nAreas
        tblidx = find(ismember(histinfo.RegionAcronym, uniqueAreas(ismember(parentArea,TheseAreas{aid2}))));
        tblidx(tblidx<starterid) = [];
        tblidx(tblidx>stopid) = [];
        newDepths(tblidx) = linspace(histdepthtmp(transitionPoints(aid2)), histdepthtmp(transitionPoints(aid2+1)), numel(tblidx));
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


function visualizeResults(Depth2Area, MUACorr, activityClusters)
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

figure('Position',[896 376 1507 848]);
subplot(1,3,1);
scatter(categorical(sortedAreas,unique(sortedAreas, 'stable')), sortedDepths, 25, areaColorsRGB, 'filled');
set(gca,'YDir','normal')
title('Ordered Depth-to-Area Mapping');
ylabel('Depth (µm)');
xlabel('Brain Area');
ylim([min(sortedDepths) max(sortedDepths)])
makepretty

subplot(1,3,2);
imagesc([],sortedDepths,activityClusters(orderIdx, :));
set(gca,'ydir','normal')
colormap(hot)
colorbar;
title('Activity Clusters');
xlabel('NMF');
ylabel('Depth');
ylim([min(sortedDepths) max(sortedDepths)])
makepretty


subplot(1,3,3);
imagesc(sortedDepths,sortedDepths,MUACorr(orderIdx, orderIdx));
set(gca,'ydir','normal')
colormap(hot)
colorbar;
title('MUA correlation');
xlabel('Depth');
ylabel('Depth');
ylim([min(sortedDepths) max(sortedDepths)])

makepretty
end
