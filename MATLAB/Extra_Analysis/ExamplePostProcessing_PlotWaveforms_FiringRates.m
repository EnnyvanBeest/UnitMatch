%% Example code to plot waveforms after UnitMatch:
GoodIdx = logical(UniqueIDConversion.GoodID); % Good units only

[UID,id1,id2] = unique(UniqueIDConversion.UniqueID(GoodIdx)); % Find unique IDs

TrackedID = find(arrayfun(@(X) sum(id2==X)>1,unique(id2))); % number of neurons tracked
if numel(TrackedID)>25 % Don't overload the plotting, just take 25 examples
    TrackedId = datasample(TrackedID,25,'replace',false);
end

figure('name','ExampleWaveforms')
for unitid = 1:numel(TrackedId)
    subplot(ceil(sqrt(numel(TrackedId))),round(sqrt(numel(TrackedId))),unitid) % subplot
    hold on
    % all neurons with this UID
    tmpwaveform = nanmean(WaveformInfo.ProjectedWaveform(:,id2==TrackedID(unitid),:),3); %Average over two halves of recording
    plot(tmpwaveform)
    xlabel('Time')
    ylabel('mV')
    makepretty
    hold off

end

%% Firing rate of a neuron across recordings
RecSes = UniqueIDConversion.recsesAll(GoodIdx); % Recording session for good units
OriClusID = UniqueIDConversion.OriginalClusID(GoodIdx); % Original ClusterID
timeBinSize = 0.1% % per 100ms
timeEdges = -timeBinSize/2:timeBinSize:max(sp.st)+timeBinSize/2;% Define timebins

figure('name','FR')
for unitid = 1:numel(TrackedId)
    subplot(ceil(sqrt(numel(TrackedId))),round(sqrt(numel(TrackedId))),unitid) % subplot
    hold on

    
    % all neurons with this UID
    Idx = find(UniqueIDConversion.UniqueID(GoodIdx) == UID(TrackedId(unitid)));

    spikesForNeuron = nan(1,numel(Idx));
    % Their recording:
    for iidx = 1:numel(Idx)
        sp = loadKSdir(UMparam.KSDir{RecSes(Idx(iidx))});
      
        spikesForNeuron(iidx) = nanmean(histcounts(sp.st(sp.clu==OriClusID(Idx(iidx))),timeEdges))./timeBinSize;
    end
    scatter(RecSes(Idx),spikesForNeuron,20,[0 0 0],'filled')
    xlabel('Recording')
    ylabel('spikes/sec')
    makepretty
    hold off

end