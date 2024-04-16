%% for Enny + Celian
mouse = 'JF067'%% Get mouse behavior 
RedoUM = 1;
TakingExampleCells = [7,8,11]; % were used for the paper
% smoothdata(PSTH, 'movmean', [10-70])
% or 0-70
nShuffle = 1000;


%% Behavior data
bhvData = cl_task_performance({mouse});

%% Recording info
% using iMouse = 1 for now, but JF078, 84 (2 probes) and 82 (3 probes) are also in the same task + all the recs are processed - so they could be integrated. 
savedirs = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\Learning_Striatum_new'
UMDays = [1:3];%[2,3,4,5,6,7];
RelevantDays = [1,2,3];%[2,4,5];%[2,7,8];%7:9;
PickChar = 'KS4';%'KS4'; %  PyKS, KS2

nRec = numel(RelevantDays);
ColOpt = flipud(gray(nRec+1));
ColOpt = ColOpt(2:end,:);
ColOpt = cat(3,ColOpt,ColOpt.*0.5);

%% Get match data

SaveDir = dir(fullfile(savedirs,mouse,'*','*','UnitMatch.mat'));
%[UniqueID, MatchTable] = AssignUniqueID_POSTUM(SaveDir);
PickID = find(cellfun(@(X) any(strfind(X,PickChar)),{SaveDir(:).folder}));
load(fullfile(SaveDir(PickID).folder,SaveDir(PickID).name))
%

% Rerun UnitMatch for relevant days only
% Get clusinfo
if RedoUM
    UMparam.KSDir = UMparam.KSDir(UMDays);
    UMparam.RawDataPaths = UMparam.RawDataPaths(UMDays);
    UMparam.AllChannelPos = UMparam.AllChannelPos(UMDays);
    % 
    % for ksid = 1:numel(UMparam.KSDir)
    %     myClusFile = dir(fullfile(UMparam.KSDir{ksid}, 'channel_map.npy'));
    %     channelmaptmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));
    % 
    %     myClusFile = dir(fullfile(UMparam.KSDir{ksid}, 'channel_positions.npy'));
    %     channelpostmp = readNPY(fullfile(myClusFile(1).folder, myClusFile(1).name));
    % 
    %     UMparam.AllChannelPos{ksid} = channelpostmp;
    % end
    UMparam.AllProbeSN = UMparam.AllProbeSN(UMDays);
    UMparam.AllDecompPaths = UMparam.AllDecompPaths(UMDays);
    UMparam.Coordinates = UMparam.Coordinates(UMDays);
    UMparam.AllRawPaths = UMparam.AllRawPaths(UMDays);
    UMparam.SaveDir = fullfile(UMparam.SaveDir,'UMDays');

    clusinfo = getClusinfo(UMparam.KSDir);
    UMparam.spikeWidth = 61; % width of spikes in samples (typically assuming 30KhZ sampling)
    UMparam.NewPeakLoc =  22; % floor(UMparam.spikeWidth./2);
    UMparam.waveidx =  UMparam.NewPeakLoc-7:UMparam.NewPeakLoc+10;
    [UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
    [UniqueIDConversion, MatchTable] = AssignUniqueID(UMparam.SaveDir);
    RelevantDays = find(ismember(UMDays,RelevantDays));

end
%% 

% % Relevant for relevant days?
idx = ismember(UniqueIDConversion.recsesAll, RelevantDays) & UniqueIDConversion.GoodID' == 1;
idx2 = ismember(UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID)),RelevantDays);
UniqueIDConversion.OriginalClusID = UniqueIDConversion.OriginalClusID(idx);
UniqueIDConversion.recsesAll = UniqueIDConversion.recsesAll(idx);
UniqueIDConversion.GoodID = UniqueIDConversion.GoodID(idx);
UniqueIDConversion.Path4UnitNPY = UniqueIDConversion.Path4UnitNPY(idx2);
UniqueIDConversion.UniqueIDConservative = UniqueIDConversion.UniqueIDConservative(idx);
UniqueIDConversion.UniqueID = UniqueIDConversion.UniqueID(idx);

WaveformInfo.MaxChannel = WaveformInfo.MaxChannel(idx2,:);
WaveformInfo.ProjectedLocation = WaveformInfo.ProjectedLocation(:,idx2,:);
WaveformInfo.ProjectedWaveform = WaveformInfo.ProjectedWaveform(:,idx2,:);
WaveformInfo.ProjectedLocationPerTP = WaveformInfo.ProjectedLocationPerTP(:,idx2,:,:);

% % Redo AssignUID for relevant subset of days
% [MatchTable, UniqueIDConversion] = AssignUniqueIDAlgorithm(MatchTable, UniqueIDConversion, UMparam);
UMparam.KSDir = UMparam.KSDir(RelevantDays);


UniqueIDConversion.UniqueIDConservative = UniqueIDConversion.UniqueID; % We take intermediate
%UniqueIDConversion = TmpFile.UniqueIDConversion;
GoodId = logical(UniqueIDConversion.GoodID);
UniqueID = UniqueIDConversion.UniqueIDConservative(GoodId);
recses = UniqueIDConversion.recsesAll(GoodId);
[UniqueIDOpt, idx1, idx2] = unique(UniqueID); %UID options
RecSesOpt = unique(recses); %Recording sessions options
RecSesPerUID = arrayfun(@(X) ismember(RecSesOpt, recses(idx2 == X)), 1:numel(UniqueIDOpt), 'Uni', 0); % Extract which recording sessions a unite appears in
RecSesPerUID = cat(2, RecSesPerUID{:});


%% Continue extracting behavior
DatesKS = cellfun(@(X) strsplit(X,'\'),UMparam.KSDir,'Uni',0);
DatesKS = cellfun(@(X) X{find(cellfun(@(X) length(X)==10,DatesKS{1}))},DatesKS,'Uni',0);
BehaviorIdx = cellfun(@(X) find(ismember(bhvData.dates,X)),DatesKS,'Uni',0); % For some reason bhdata is +1 indexed
BehaviorIdx2Take = cell2mat(BehaviorIdx(~cellfun(@isempty,BehaviorIdx)));
BehaviorIdx = find(~cellfun(@isempty,BehaviorIdx));
Performance = nan(1,nRec);
Performance(BehaviorIdx) = bhvData.goLeft(BehaviorIdx2Take, 1) ./ bhvData.nTrials(BehaviorIdx2Take, 1);

figure('name','Behavior');
scatter(1:nRec,Performance,35,ColOpt(1:nRec,:,1),'filled')
ylim([0.5 1])
set(gca,'XTick',[1:nRec],'XTickLabel',DatesKS)
makepretty
offsetAxes

%%
% units to get PSTH for 
Or_UniqueID = arrayfun(@(x) {UniqueIDConversion.Path4UnitNPY{x}(end -17:end - 14)}, 1:size(UniqueIDConversion.Path4UnitNPY, 2));
UUIDs = UniqueID(idx1);
RecOpt = unique(recses);
allRecordings = recses(idx1);
oriID_good = UniqueIDConversion.OriginalClusID(GoodId);
theseUnits = find(sum(RecSesPerUID) >= 1); % get all cells 
UniqueIDConversion.Path4UnitNPY_noGoodID = cell(size(UniqueIDConversion.UniqueIDConservative, 2), 1);
UniqueIDConversion.Path4UnitNPY_noGoodID(GoodId) = UniqueIDConversion.Path4UnitNPY;

% experiment info 
animal = mouse;
protocol = 'choiceworld'; % (this is the name of the Signals protocol)
experiments = cl_find_experiments(animal, protocol, true);
experiments = experiments([experiments.ephys]);

% Find experiments that we want to keep
ExpIndx = cellfun(@(X) find(ismember({experiments.thisDate},X)),DatesKS);
experiments = experiments(ExpIndx);
for iRecording = 1:length(experiments)
    experiments(iRecording).ephys_ks_paths = UMparam.KSDir{iRecording};
end
loadClusters = 0; % whether to load phy output or not 

%% Get PSTHs for matched cells 

% PSTH and ACG parameters
raster_window = [-0.5, 1];
psth_bin_size = 0.001;
ACGbinSize = 0.001;
ACGduration = 1;
maxnTrials = 2500;
ISIbins = [0 5*10.^(-4:0.2:0)];

% PSTH time 
psth_time = raster_window(1)+psth_bin_size/2:psth_bin_size:raster_window(2)-psth_bin_size/2;
save(fullfile(savedirs,mouse,'psth_time.mat'), 'psth_time')

% pre-allocate array space 
acg_long_track = nan(size(theseUnits, 2), nRec, 500);
vis_long_track_pass = nan(size(theseUnits, 2), nRec, 3, 1500);
dPrimePerUnit = nan(size(theseUnits, 2), nRec);
VisRespPerUnit = nan(size(theseUnits,2), nRec, 3, maxnTrials);
ISIsPerUnit = nan(size(theseUnits,2), nRec, length(ISIbins)-1);


UnitIsMerged = false(size(theseUnits, 2),1);

for iRecording = 1:nRec
    thisRecording = iRecording;
    site = 1;
    recording = [];
    thisDate = experiments(thisRecording).thisDate;
    n_trials = zeros(size(experiments(thisRecording).experiment, 2), 1);
    if n_trials>maxnTrials
        keyboard
    end
    keepMe = zeros(size(experiments(thisRecording).experiment, 2), 1);

    for iExperiment = 1:size(experiments(thisRecording).experiment, 2)
        exp = experiments(thisRecording).experiment(iExperiment);
        % cl_cortexlab_filename(animal, thisDate, experiment, file, site, recording, shank, KSVersion)
        [block_filename, block_exists] = cl_cortexlab_filename(animal, thisDate, exp, 'block');
        try % QQ hacky 
            load(block_filename)
            keepMe(iExperiment) = contains(block.expDef, 'choiceworld');
            if isfield(block.events, 'stim_idValues')
                n_trials(iExperiment) = length(block.events.stim_idValues);
            elseif isfield(block.events, 'stimulusOnTimes')
                n_trials(iExperiment) = length(block.events.stimulusOnTimes);
            end
        catch
            n_trials(iExperiment) = NaN;
            keepMe(iExperiment) = false;
        end
    end

    if sum(keepMe) == 0
        continue; 
    end

    experiment = experiments(thisRecording).experiment(1);
    if length(experiment) > 1
        experiment = 2;
    end

    % filename = cl_cortexlab_filename(animal, thisDate, '', 'ephys', 1, '', '');
    filename = UMparam.KSDir{thisRecording};

    try % QQ hacky 
        cl_load_experiment;
    catch ME
        warning('error')
        continue;
    end

    spike_templates_0dx_unique = unique(spike_templates_0idx);

    for iUnit = 1:size(theseUnits, 2)
      
        thisUID = UUIDs(theseUnits(iUnit));
        thisMatchTableIdx = GoodId & UniqueIDConversion.UniqueIDConservative == thisUID;
        recordings_unique = unique([UniqueIDConversion.recsesAll(thisMatchTableIdx)]);
        if ~ismember(RecOpt(iRecording), recordings_unique)
            continue;
        end

        thisUnit_0idx = UniqueIDConversion.OriginalClusID(GoodId' & ...
            UniqueIDConversion.UniqueIDConservative' == thisUID & ...
            UniqueIDConversion.recsesAll == RecOpt(thisRecording));

        if numel(thisUnit_0idx)>1
            UnitIsMerged(iUnit)=1;
        end
        theseSpikeTimes = spike_times_timeline(ismember(uint32(spike_templates_0idx), uint32(thisUnit_0idx)));
        if isempty(theseSpikeTimes)
            keyboard
        end
        % get ACG
        [acg, ~] = CCGBz([double(theseSpikeTimes); double(theseSpikeTimes)], [ones(size(theseSpikeTimes, 1), 1); ...
            ones(size(theseSpikeTimes, 1), 1) * 2], 'binSize', ACGbinSize, 'duration', ACGduration, 'norm', 'rate'); %function
        ACG = acg(:, 1, 1);
        acg_long_track(iUnit, thisRecording, :) = ACG(501:1000);

        % StimOffTimes = stimOn_times+1;
        % exclst = arrayfun(@(X,Y) theseSpikeTimes<X|theseSpikeTimes>Y,stimOn_times,StimOffTimes,'Uni',0);
        % exclst = any(cat(2,exclst{:}),2);
        % theseSpikeTimes(exclst) = nan;
        ISI2 = histcounts(diff(theseSpikeTimes),ISIbins, 'Normalization','probability');
        ISIsPerUnit(iUnit,thisRecording,:) = ISI2./sum(ISI2);

        % get visual PSTH
        [align_group_a, align_group_b] = ismember(trial_conditions(:, 2), unique(trial_conditions(:, 2)));
        [curr_psth, curr_raster, t, raster_x, raster_y] = cl_raster_psth(double(spike_templates_0idx), double(spike_times_timeline), ...
            double(thisUnit_0idx), raster_window, psth_bin_size, stimOn_times, align_group_b(1:size(stimOn_times, 1)));
        vis_long_track_pass(iUnit, thisRecording, 1:size(curr_psth, 1), :) = curr_psth;

        tmpl = nansum(curr_raster(align_group_b==1,psth_time>0&psth_time<0.5),2).*2; %sp/sec
        tmpc = nansum(curr_raster(align_group_b==2,psth_time>0&psth_time<0.5),2).*2;% sp/sec
        dPrimePerUnit(iUnit,thisRecording) = (nanmean(tmpc,1)-nanmean(tmpl,1))./...
            (0.5*sqrt(nanvar(tmpl,[],1)+nanvar(tmpc,[],1)));

        VisRespPerUnit(iUnit, thisRecording, 1, 1:numel(tmpl)) = tmpl;
        VisRespPerUnit(iUnit, thisRecording, 2, 1:numel(tmpc)) = tmpc;
    end
end
vis_long_track_pass(:,:,3,:) = []; %remove 90* stimulus, did not always show it. 
%% save data
% unit numbers 
%% All cells
% vis_long_track_pass = vis_long_track_pass.*1000; % spike count per bin
% vis_long_track_pass = vis_long_track_pass(:,RelevantDays,:,:);

unit_UniqueIDs = UUIDs(theseUnits);
save(fullfile(savedirs,mouse,'all_units_uniqueID.mat'), 'unit_UniqueIDs')

% PSTHs : unit_number x recording_number x stimulus_type x time (samples). 
%   stimulus_type = 1 means a controlateral stimulus
%   stimulus type = 2 means a central stimulus. 
save(fullfile(savedirs,mouse,'all_units_psth.mat'), 'vis_long_track_pass')


%%

nClus = size(vis_long_track_pass,1);
nStim = 2;
stepsz = 0.05; % in ms
binEdges = [round(min(psth_time)./stepsz).*stepsz:stepsz:round(max(psth_time)./stepsz).*stepsz];
timeline = round(min(psth_time)./stepsz).*stepsz+0.5.*stepsz:stepsz:round(max(psth_time)./stepsz).*stepsz-stepsz./2;
[bincount,idx1,idx2] = histcounts(psth_time,binEdges); % idx2 ident,ifies which bin each sample belongs
% smoothdata(PSTH, 'movmean', [10-70])

PSTH = smoothdata(vis_long_track_pass,4,'movmean',[10 70]);
% PSTH = arrayfun(@(X) nansum(PSTH(:,:,:,idx2==X),4)./stepsz,unique(idx2),'Uni',0);
% PSTH = cat(4,PSTH{:});
timeline = psth_time;
PSTH_Z = (PSTH - nanmean(PSTH,4))./(nanstd(PSTH,[],4));

%% Neurons to include
[nPres,sortidx] = sort(nansum(PSTH(:,:,1,5)~=0 & ~isnan(PSTH(:,:,1,5)),2),'descend');
sortidx = sortidx(nPres>=nRec);
sortidx(UnitIsMerged(sortidx)) = [];
% sortidx = sort(sortidx);
% sortidx = [89,97,102,103,105,160,184]; % For JF067
figure;
for recid = 1:nRec
    subplot(2,nRec,recid)
    h=imagesc(timeline,[],squeeze(PSTH_Z(sortidx,recid,1,:)),[-5 5]);
    set(h,'AlphaData',~isnan(squeeze(PSTH_Z(sortidx,recid,1,:))))
    title(['Day ' num2str(RelevantDays(recid))])
    if recid>1
        set(gca,'YTickLabel',[])
    end
    set(gca,'XTickLabel',[])
    makepretty
    offsetAxes

    subplot(2,nRec,recid+nRec)
    h=imagesc(timeline,[],squeeze(PSTH_Z(sortidx,recid,2,:)),[-5 5]);
    set(h,'AlphaData',~isnan(squeeze(PSTH_Z(sortidx,recid,1,:))))
    if recid>1
        set(gca,'YTickLabel',[])
    end
    
    makepretty
    offsetAxes
end
colormap(redblue)

%% Match probability vs baselinen firing
baselinefr = nanmean(nanmean(PSTH(sortidx,:,:,timeline<0),4),3);
baselinefrdiff = nan(numel(sortidx),nRec);
MatchProb = nan(numel(sortidx), nRec);
countid = 1;

for recid = 1:nRec-1
    for recid2 = recid+1:nRec
        for uid = 1:numel(sortidx)
            thisU = unit_UniqueIDs(sortidx(uid));
            tblidx1 = MatchTable.UID1==thisU & MatchTable.UID2 == thisU & MatchTable.RecSes1 == recid & MatchTable.RecSes2 == recid2;
            tblidx2 = MatchTable.UID1==thisU & MatchTable.UID2 == thisU & MatchTable.RecSes1 == recid2 & MatchTable.RecSes2 == recid;
            MatchProb(uid,countid) = nanmean([MatchTable.MatchProb(tblidx1); MatchTable.MatchProb(tblidx2)]);

            disp([num2str(sortidx(uid)) ' day ' num2str(recid) ' versus' num2str(recid2) ' p=' num2str(MatchProb(uid,countid))])
        end
        baselinefrdiff(:,countid) = baselinefr(:,recid2)-baselinefr(:,recid);

        countid = countid + 1;

    end
end

%% Baseline firing rates

T = array2table(baselinefr);
withinDesign = table([1:nRec]','VariableNames',{'Day'});
withinDesign.Day = categorical(withinDesign.Day);
% Repeated measures model
rm = fitrm(T,['baselinefr1-baselinefr' num2str(size(T,2)) ' ~ 1'],'WithinDesign',withinDesign);

AT = ranova(rm,'WithinModel','Day')
multcompare(rm,'Day')

MatchFig = figure('name','MatchProb vs Functional');
subplot(4,1,1)
hold on
for countid = 1:size(MatchProb,2)
    scatter(MatchProb(:,countid),baselinefrdiff(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(MatchProb(TakingExampleCells,countid),baselinefrdiff(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))

end
hold on
line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
ylabel('Baseline firing')
% set(gca,'xscale','log')
makepretty
offsetAxes

[r,pval] = corr(MatchProb(:),baselinefrdiff(:));
pShuf = nan(1,nShuffle);
for shufid = 1:nShuffle
    [pShuf(shufid)] = corr(MatchProb(:),datasample(baselinefrdiff(:),numel(baselinefrdiff(:))));
end
adjPval = invprctile(pShuf,r)/100;
disp(['Adjusted p=' num2str(adjPval) ' for firing rate'])
title(['p= ' num2str(adjPval)])

%% population
dRR = (PSTH-nanmean(PSTH(:,:,:,timeline<0),4))./nanmean(PSTH(:,:,:,timeline<0),4);
clear h
figure('name','Population');
ylims = [-1 5];
%   stimulus_type = 1 means a lateral stimulus
%   stimulus type = 2 means a central stimulus. 
for recid = 1:nRec
    subplot(1,nRec,recid)
    h(1) = plot(timeline,squeeze(nanmean(dRR(sortidx,recid,2,:),1)),'color',ColOpt(recid,:,2),'LineWidth',3);
    hold on
    h(2) = plot(timeline,squeeze(nanmean(dRR(sortidx,recid,1,:),1)),'color',ColOpt(recid,:,1),'LineWidth',1);

    line([0 0],ylims,'color',[0.5 0.5 0.5])
    title(['Day ' num2str(RelevantDays(recid))])
    if recid==1
        ylabel('dRR')
    else
        axis off
    end
    xlabel('time (s)')
    ylim(ylims)

    offsetAxes
    makepretty
    legend([h(:)],{'Central','Lateral'})


end



%% Example cells
clear h
% sortidx = [63, 97, 102, 103, 105, 111];
nExample = 4;
nExampleOri = nExample;
nBatch = ceil(numel(sortidx)/nExample);
nCols = 3 + nRec;    

ISIstoSave = nan(numel(sortidx),length(ISIbins)-1,nRec);
RespToSave = nan(numel(sortidx),length(timeline),nRec,2);
% ExampleID = [478,215,279]
%   stimulus_type = 1 means a lateral stimulus
%   stimulus type = 2 means a central stimulus.
ExclSingU = false(1,numel(sortidx));
for batchid = 1:nBatch
    figure('name',['Example Cells, batch ' num2str(batchid) '/' num2str(nBatch)]);
    takeindx = (batchid-1)*nExample+1:batchid*nExample;
    takeindx(takeindx>numel(sortidx)) = [];
    ExampleID = sortidx(takeindx);%datasample(sortidx,nExample,'replace',false);
    % Waveforms of example cells
    tmprecses = UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID));
    OriID = UniqueIDConversion.OriginalClusID(logical(UniqueIDConversion.GoodID));
    UID = UniqueIDConversion.UniqueIDConservative(logical(UniqueIDConversion.GoodID));
    UID2take = unit_UniqueIDs(ExampleID);
    nExample = numel(ExampleID);
    for exid = 1:nExample
        % idx = find(ismember(UniqueIDConversion.UniqueIDConservative(logical(UniqueIDConversion.GoodID)),UID2take(exid)));
        % idx(~ismember(tmprecses(idx),RelevantDays)) = [];
        % tmpday = tmprecses(idx);

        thisUID = UUIDs(sortidx(exid));
        
        subplot(nExample,nCols,(exid-1)*nCols+1)
        hold on
        for did = 1:nRec    
            thiswavIdx = find(UniqueIDConversion.UniqueIDConservative(find(GoodId)) == thisUID & UniqueIDConversion.recsesAll(find(GoodId))'==RecOpt(did));

            plot(squeeze(nanmean(nanmean(WaveformInfo.ProjectedWaveform(:,thiswavIdx,:),3),2)),'color',ColOpt(did,:,1));
            % plot(squeeze(waveforms_raw_long_track_enny(ExampleID(exid),did,:)),'color',ColOpt(did,:,1));
        end
        if exid<nExample
            axis off
        end
        ylim([-70 30])
        title(num2str(ExampleID(exid)))

        makepretty
        offsetAxes

        subplot(nExample,nCols,(exid-1)*nCols+2)
        hold on
        thiswavIdx = find(UniqueIDConversion.UniqueIDConservative(find(GoodId)) == thisUID);
        waveidx = find(sum(~isnan(squeeze(nanmean(WaveformInfo.ProjectedLocationPerTP(2,thiswavIdx,:,:),4))),1)==numel(thiswavIdx));

        for did = 1:nRec
            thiswavIdx = find(UniqueIDConversion.UniqueIDConservative(find(GoodId)) == thisUID & UniqueIDConversion.recsesAll(find(GoodId))'==RecOpt(did));
            plot(squeeze(nanmean(nanmean(WaveformInfo.ProjectedLocationPerTP(2,thiswavIdx,waveidx,:),2),4)),squeeze(nanmean(nanmean(WaveformInfo.ProjectedLocationPerTP(3,thiswavIdx,waveidx,:),2),4)),'color',ColOpt(did,:,1))
        end
        if exid<nExample
            axis off
        end
        makepretty
        offsetAxes
    end
    %

    UID2take = unit_UniqueIDs(ExampleID);
    for exid = 1:nExample
        subplot(nExample,nCols,(exid-1)*nCols+3)
        hold on
        for did = 1:nRec
            stairs(ISIbins(1:end-1)*1000,smooth(squeeze(ISIsPerUnit(ExampleID(exid),did,:)),5),'color',ColOpt(did,:,1), 'LineWidth', 2.0);
            ISIstoSave((batchid-1)*nExampleOri+exid,:,did) =squeeze(ISIsPerUnit(ExampleID(exid),did,:));
        end

        tmpcorr = corr(squeeze(ISIsPerUnit(ExampleID(exid),:,:))',squeeze(ISIsPerUnit(ExampleID(exid),:,:))');
        if any(tmpcorr(:))<0.5
            title('Low ISI corr)')
            ExclSingU((batchid-1)*nExampleOri+exid) = true;
        end
        xticks([5 50 500])
        yticks([0 0.1])
        xlabel('Time (ms)')
        ylabel('Firing rate (sp/s)')
        ylim([0 0.2])
        set(gca,'XScale','log')
        if exid<nExample
            axis off
        end
        makepretty
    end

    for exid = 1:numel(ExampleID)
        ylims = [floor(nanmin(nanmin(nanmin(dRR(ExampleID(exid),:,:,:))))-0.2) ceil(nanmax(nanmax(nanmax(dRR(ExampleID(exid),:,:,:))))+0.2)];

        for recid = 1:nRec

            subplot(nExample,nCols,(exid-1)*nCols+3+recid)

            hold on
            if ~any(dRR(ExampleID(exid),recid,1,:)~=0)
                continue
            end
            h(2) = plot(timeline,(squeeze(dRR(ExampleID(exid),recid,1,:))),'color',ColOpt(recid,:,1));
            h(1) = plot(timeline,(squeeze(dRR(ExampleID(exid),recid,2,:))),'color',ColOpt(recid,:,2));

            RespToSave((batchid-1)*nExampleOri+exid,:,recid,1) = (squeeze(dRR(ExampleID(exid),recid,1,:))); % Lateral
            RespToSave((batchid-1)*nExampleOri+exid,:,recid,2) = (squeeze(dRR(ExampleID(exid),recid,2,:))); %Central

            line([0 0],ylims,'color',[0.7 0.7 0.7])
            if exid==1
                title(['Day ' num2str(RelevantDays(recid))])
            end
            if recid==1
                ylabel('dRR')
            end
            if exid == nExample
                xlabel('time (s)')
            end
            if ~(recid==1)
                axis off
            end
            offsetAxes
            makepretty
            ylim(ylims)
        end

    end
end



%% AUC values
redCol = cat(2,linspace(0.3,1,3)',repmat(0,3,1),repmat(0,3,1));
blueCol = cat(2,repmat(0,3,1),repmat(0,3,1),linspace(0.3,1,3)');

AUCISI = nan(nRec,nRec);
AUCCentral = nan(nRec,nRec);
AUClateral = nan(nRec,nRec);
figure('name','AUC values')
countid = 1;
for did1 = 1:nRec
    for did2 = 1:nRec
        if did1>=did2
            continue
        end
        ISICorr = corr(squeeze(ISIstoSave(:,:,did1))',squeeze(ISIstoSave(:,:,did2))');
        scores = [diag(ISICorr);ISICorr(logical(triu(ones(size(ISICorr)),1)))];
        labels = [ones(size(ISICorr,1),1); zeros(sum(sum(triu(ones(size(ISICorr)),1))),1)];
        [~,~,~,AUCISI(did1,did2)] = perfcurve(labels,scores,1);

        LatCorr = corr(squeeze(RespToSave(:,:,did1,1))',squeeze(RespToSave(:,:,did2,1))');
        scores = [diag(LatCorr);LatCorr(logical(triu(ones(size(LatCorr)),1)))];
        labels = [ones(size(LatCorr,1),1); zeros(sum(sum(triu(ones(size(LatCorr)),1))),1)];
        [~,~,~,AUClateral(did1,did2)] = perfcurve(labels,scores,1);

        CentCorr = corr(squeeze(RespToSave(:,:,did1,2))',squeeze(RespToSave(:,:,did2,2))');
        scores = [diag(CentCorr);CentCorr(logical(triu(ones(size(CentCorr)),1)))];
        labels = [ones(size(CentCorr,1),1); zeros(sum(sum(triu(ones(size(CentCorr)),1))),1)];
        [~,~,~,AUCCentral(did1,did2)] = perfcurve(labels,scores,1);


        h(1) = scatter(AUCISI(did1,did2),AUCCentral(did1,did2),35,redCol(countid,:),'filled');
        hold on;
        h(2) = scatter(AUCISI(did1,did2),AUClateral(did1,did2),35,blueCol(countid,:),'filled');
        countid = countid+1;

    end
end


xlabel('ISI')
ylabel('Visual response')
line([0 1],[0 1],'color',[0.2 0.2 0.2])

xlim([0 1])
ylim([0 1])
legend({'Central','Lateral'})
makepretty
offsetAxes

%%
RespToSave = (RespToSave-nanmin(reshape(RespToSave,numel(sortidx),[]),[],2))./(nanmax(reshape(RespToSave,numel(sortidx),[]),[],2)-nanmin(reshape(RespToSave,numel(sortidx),[]),[],2));
ISIstoSave = (ISIstoSave-nanmin(reshape(ISIstoSave,numel(sortidx),[]),[],2))./(nanmax(reshape(ISIstoSave,numel(sortidx),[]),[],2)-nanmin(reshape(ISIstoSave,numel(sortidx),[]),[],2));

%%
CorrVals = nan(numel(sortidx),nRec);
LatVals = nan(numel(sortidx),nRec);
CentVals = nan(numel(sortidx),nRec);
countid = 1;
for did1 = 1:nRec-1
    for did2 = did1+1:nRec;
        % ISICorr = nanmean((squeeze(ISIstoSave(:,:,did1))'-squeeze(ISIstoSave(:,:,did2))').^2);
        % CorrVals(:,countid) = (ISICorr);

        ISICorr = corr(squeeze(ISIstoSave(:,:,did1))',squeeze(ISIstoSave(:,:,did2))');
        CorrVals(:,countid) = diag(ISICorr);

        LatCorr = nanmean((squeeze(RespToSave(:,:,did1,1))'-squeeze(RespToSave(:,:,did2,1))').^2);
        LatVals(:,countid) = (LatCorr);

        % LatCorr = corr(squeeze(RespToSave(:,:,did1,1))',squeeze(RespToSave(:,:,did2,1))');
        % LatVals(:,countid) = diag(LatCorr);

        CentCorr = nanmean((squeeze(RespToSave(:,:,did1,2))'-squeeze(RespToSave(:,:,did2,2))').^2);
        CentVals(:,countid) = (CentCorr);

        %  CentCorr = corr(squeeze(RespToSave(:,:,did1,2))',squeeze(RespToSave(:,:,did2,2))');
        % CentVals(:,countid) = diag(CentCorr);
        countid = countid+1;
    end
end


[r,pval] = corr(MatchProb(:),CorrVals(:));
pShuf = nan(1,nShuffle);
for shufid = 1:nShuffle
    [pShuf(shufid)] = corr(MatchProb(:),datasample(CorrVals(:),numel(CorrVals(:))));
end
    adjPval = invprctile(pShuf,r)/100;
    disp(['Adjusted p=' num2str(adjPval) ' for ISI'])


figure(MatchFig)
subplot(4,1,2)
hold on
for countid = 1:size(MatchProb,2)
    scatter(MatchProb(:,countid),CorrVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(MatchProb(TakingExampleCells,countid),CorrVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
title(['p= ' num2str(adjPval)])

line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
ylabel('ISI')
makepretty
offsetAxes

[r,pval] = corr(MatchProb(:),LatVals(:));
pShuf = nan(1,nShuffle);
for shufid = 1:nShuffle
    [pShuf(shufid)] = corr(MatchProb(:),datasample(LatVals(:),numel(LatVals(:))));
end
    adjPval = invprctile(pShuf,r)/100;
    disp(['Adjusted p=' num2str(adjPval) ' for Lateral visual response'])

subplot(4,1,3)
hold on
for countid = 1:size(MatchProb,2)
    scatter(MatchProb(:,countid),LatVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(MatchProb(TakingExampleCells,countid),LatVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
% ylim([0 0.1])
title(['p= ' num2str(adjPval)])

ylabel('Lateral visual response')
makepretty
offsetAxes

[r,pval] = corr(MatchProb(:),CentVals(:));
pShuf = nan(1,nShuffle);
for shufid = 1:nShuffle
    [pShuf(shufid)] = corr(MatchProb(:),datasample(CentVals(:),numel(CentVals(:))));
end
    adjPval = invprctile(pShuf,r)/100;
    disp(['Adjusted p=' num2str(adjPval) ' for Central visual response'])


subplot(4,1,4)
hold on
for countid = 1:size(MatchProb,2)
    scatter(MatchProb(:,countid),CentVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(MatchProb(TakingExampleCells,countid),CentVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
title(['p= ' num2str(adjPval)])

line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
% ylim([0 1.5])
xlabel('p(match)')
ylabel('Central visual response')
makepretty
offsetAxes
%%
figure('name','ISI vsLateral')
subplot(3,1,1)
hold on
for countid = 1:size(MatchProb,2)
    scatter(CorrVals(:,countid),CentVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(CorrVals(TakingExampleCells,countid),CentVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
% ylim([0 1.5])
xlabel('ISI')
ylabel('Central')
makepretty
offsetAxes

subplot(3,1,2)
hold on
for countid = 1:size(MatchProb,2)
    scatter(CorrVals(:,countid),LatVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(CorrVals(TakingExampleCells,countid),LatVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
% ylim([0 1.5])
xlabel('ISI')
ylabel('Lateral')
makepretty
offsetAxes

subplot(3,1,3)
hold on
for countid = 1:size(MatchProb,2)
    scatter(CentVals(:,countid),LatVals(:,countid),35,ColOpt(countid,:,1),'filled')
    scatter(CentVals(TakingExampleCells,countid),LatVals(TakingExampleCells,countid),37,hsv(length(TakingExampleCells)))
end
line(get(gca,'xlim'),[0 0],'LineStyle',':','color',[0.1 0.1 0.1])
% ylim([0 1.5])
xlabel('Central')
ylabel('Lateral')
makepretty
offsetAxes
%% Conclusion about population?
% Ratio Lateral-Central
%   stimulus_type = 1 means a lateral stimulus
%   stimulus type = 2 means a central stimulus. 
ExampleID = sortidx;
ColOpt = ColOpt(1:nRec,:,1); % only need one

NeuronCols = distinguishable_colors(numel(sortidx));
ExIds = find(ismember(sortidx,ExampleID));

figure('name','dprime'); 
subplot(2,1,1)
hold on
tmpl = squeeze(nanmean(dRR(sortidx,:,1,timeline>0&timeline<0.5),4)); % tmp l
tmpc = squeeze(nanmean(dRR(sortidx,:,2,timeline>0&timeline<0.5),4)); % tmp c
dprime = (nanmean(tmpc,1)-nanmean(tmpl,1))./...
    (0.5*sqrt(nanvar(tmpl,[],1)+nanvar(tmpc,[],1)));
scatter(RelevantDays,dprime,35,ColOpt,'filled')
line([RelevantDays(1) RelevantDays(end)],[0 0],'Color',[0.2 0.2 0.2])
xlim([RelevantDays(1)-0.5 RelevantDays(end)+0.5])
xlabel('Day')
ylabel('dPrime')
title('Population')
makepretty
offsetAxes
%
% dPrimePerUnit(iUnit,thisRecording) = (nanmean(tmpc,1)-nanmean(tmpl,1))./...
%           (0.5*sqrt(nanvar(tmpl,[],1)+nanvar(tmpc,[],1)));
%
%       VisRespPerUnit(iUnit, thisRecording, 1, 1:numel(tmpl)) = tmpl;
%       VisRespPerUnit(iUnit, thisRecording, 2, 1:numel(tmpc)) = tmpc;


Sig = zeros(1,numel(sortidx));
for uid = 1:numel(sortidx)
    tmpl = squeeze(VisRespPerUnit(sortidx(uid),:,1,:));
    tmpc = squeeze(VisRespPerUnit(sortidx(uid),:,2,:));

    % Difference
    tmpdiff = tmpc - tmpl;
    tmpl(:,sum(isnan(tmpl),1)==nRec) = [];
    tmpc(:,sum(isnan(tmpc),1)==nRec) = [];
    tmpdiff(:,sum(isnan(tmpdiff),1)==nRec) = [];

    tmpdiff = cat(1,tmpdiff(1,:),nanmean(tmpdiff([2:3],:),1));

    p = anova1(tmpdiff',[],'off');

    if p<0.05
        if sum(sign(diff(nanmean(tmpdiff,2))))>0
            Sig(uid) = 1;
            dPrimePerUnit(sortidx(uid),:)

        elseif sum(sign(diff(nanmean(tmpdiff,2))))<0
            Sig(uid) = -1;
        end
    end


end
subplot(2,1,2)
hold on
for uid = 1:numel(sortidx)
    if Sig(uid)==-1
        plot(1:nRec, dPrimePerUnit(sortidx(uid),:),'.-','color',[0 0.2 1],'MarkerSize',30)
    elseif Sig(uid)==1
        plot(1:nRec, dPrimePerUnit(sortidx(uid),:),'.-','color',[1 0.2 0],'MarkerSize',30)
    else
        plot(1:nRec, dPrimePerUnit(sortidx(uid),:),'.-','color',[0 0 0],'MarkerSize',30)
    end
end
xlabel('Day')
ylabel('dPrime')
title('Individual neurons')
makepretty
offsetAxes

figure('name','Pie')
piechart(categorical(Sig))



%% Correlation of functional scores
tmpl = squeeze(dRR(sortidx,:,1,:)); % tmp l
tmpc = squeeze(dRR(sortidx,:,2,:)); % tmp c


tmplcorr = nan(numel(sortidx),nRec,nRec);
tmpccorr = nan(numel(sortidx),nRec,nRec);
ISIcorr = nan(numel(sortidx),nRec,nRec);
counter = 0;

figure('name','Functional Correlations')
for did1 = 1:nRec
    for did2 = 1:nRec
        if did1>=did2
            continue
        end

        for uid = 1:numel(sortidx)
            tmplcorr(uid,did1,did2) = corr(squeeze(tmpl(uid,did1,:)),squeeze(tmpl(uid,did2,:)));
            tmpccorr(uid,did1,did2) = corr(squeeze(tmpc(uid,did1,:)),squeeze(tmpc(uid,did2,:)));
            ISIcorr(uid,did1,did2) = corr(squeeze(ISIstoSave(uid,:,did1))',squeeze(ISIstoSave(uid,:,did2))');
        end
        counter = counter+1;
        subplot(3,1,1)
        hold on
        scatter(repmat(counter,1,numel(sortidx)),tmplcorr(:,did1,did2),30,nanmean(ColOpt([did1,did2],:),1),'filled')
        title('Lateral Stimulus response')



        subplot(3,1,2)
        hold on
        scatter(repmat(counter,1,numel(sortidx)),tmpccorr(:,did1,did2),30,nanmean(ColOpt([did1,did2],:),1),'filled')
        title('Central Stimulus response')
   


        subplot(3,1,3)
        hold on
        scatter(repmat(counter,1,numel(sortidx)),ISIcorr(:,did1,did2),30,nanmean(ColOpt([did1,did2],:),1),'filled')
        title('ISIh')
   

    end
end
for subid = 1:3
    subplot(3,1,subid)
    ylabel('Correlation')
         makepretty
        offsetAxes
end
linkaxes


% Normalize for MSE
tmpl = reshape(tmpl,numel(sortidx),[]);
tmpl = (tmpl-nanmin(tmpl,[],2))./(nanmax(tmpl,[],2) - nanmin(tmpl,[],2));
tmpl = reshape(tmpl,numel(sortidx),nRec,[]);

tmpc = reshape(tmpc,numel(sortidx),[]);
tmpc = (tmpc-nanmin(tmpc,[],2))./(nanmax(tmpc,[],2) - nanmin(tmpc,[],2));
tmpc = reshape(tmpc,numel(sortidx),nRec,[]);


ISIstoSave = reshape(ISIstoSave,numel(sortidx),[]);
ISIstoSave = (ISIstoSave-nanmin(ISIstoSave,[],2))./(nanmax(ISIstoSave,[],2) - nanmin(ISIstoSave,[],2));
ISIstoSave = reshape(ISIstoSave,numel(sortidx),[],nRec);

%% Shuffle
tmpl = squeeze(dRR(:,:,1,:)); % tmp l
tmpc = squeeze(dRR(:,:,2,:)); % tmp c

tmplcorr = nan(numel(sortidx),nRec-1);
tmpccorr = nan(numel(sortidx),nRec-1);
ISIcorr = nan(numel(sortidx),nRec-1);

tmplcorrShuf = nan(nShuffle,nRec-1);
tmpccorrShuf = nan(nShuffle,nRec-1);
ISIcorrShuf = nan(nShuffle,nRec-1);
counter = 0;

for did1 = 1:nRec-1
    for did2 = did1+1
        if did1>=did2
            continue
        end
        counter = counter+1;

        for uid = 1:numel(sortidx)
            % tmplcorr(uid,counter) = nanmean((squeeze(tmpl(uid,did1,:))-squeeze(tmpl(uid,did2,:))).^2);
            % tmpccorr(uid,counter) = nanmean((squeeze(tmpc(uid,did1,:))-squeeze(tmpc(uid,did2,:))).^2);
            % ISIcorr(uid,counter) = nanmean((squeeze(ISIstoSave(uid,:,did1))-squeeze(ISIstoSave(uid,:,did2))).^2);
            % %             % 
            tmplcorr(uid,counter) = corr(squeeze(tmpl(sortidx(uid),did1,:)),squeeze(tmpl(sortidx(uid),did2,:)),'Type','pearson');
            tmpccorr(uid,counter) = corr(squeeze(tmpc(sortidx(uid),did1,:)),squeeze(tmpc(sortidx(uid),did2,:)),'Type','pearson');
            ISIcorr(uid,counter) = corr(squeeze(ISIsPerUnit(sortidx(uid),(did1),:)),squeeze(ISIsPerUnit(sortidx(uid),(did2),:)),'Type','pearson');
        end
        for shufid = 1:nShuffle
            uid1 = datasample(find(~isnan(nanmean(tmpl(:,did1,:),3))),1);
            uid2 = datasample(find(~isnan(nanmean(tmpl(:,did2,:),3))),1);
            % tmplcorrShuf(shufid,counter) = nanmean((squeeze(tmpl(uid1,did1,:))-squeeze(tmpl(uid(2),did2,:))).^2);
            % tmpccorrShuf(shufid,counter) = nanmean((squeeze(tmpc(uid1,did1,:))-squeeze(tmpc(uid(2),did2,:))).^2);
            % ISIcorrShuf(shufid,counter) = nanmean((squeeze(ISIstoSave(uid1,:,did1))-squeeze(ISIstoSave(uid(2),:,did2))).^2);
            % %             % 
    
            tmplcorrShuf(shufid,counter) = corr(squeeze(tmpl(uid1,did1,:)),squeeze(tmpl(uid2,did2,:)),'Type','pearson');
            tmpccorrShuf(shufid,counter) = corr(squeeze(tmpc(uid1,did1,:)),squeeze(tmpc(uid2,did2,:)),'Type','pearson');
            ISIcorrShuf(shufid,counter) = corr(squeeze(ISIsPerUnit(uid1,(did1),:)),squeeze(ISIsPerUnit(uid2,(did2),:)),'Type','pearson');
        end
  

    end
end
figure('name','Functional Scores')

subplot(3,1,1)
hold on; 
% h=violinplot(tmplcorrShuf); 
plot(1:counter,quantile(tmplcorrShuf,0.05,1),'-','color',[0.5 0.5 0.5])
plot(1:counter,quantile(tmplcorrShuf,0.95,1),'-','color',[0.5 0.5 0.5])
scatter(1:counter,nanmean(tmplcorr,1),50,[0 0 0],'filled')
title('Lateral stimulus')
makepretty
offsetAxes
ylabel('Corr')

subplot(3,1,2)
hold on; 
% h=violinplot(tmplcorrShuf); 
plot(1:counter,quantile(tmpccorrShuf,0.05,1),'-','color',[0.5 0.5 0.5])
plot(1:counter,quantile(tmpccorrShuf,0.95,1),'-','color',[0.5 0.5 0.5])
scatter(1:counter,nanmean(tmpccorr,1),50,[0 0 0],'filled')
title('Central stimulus')

makepretty
offsetAxes
ylabel('Corr')

subplot(3,1,3)
hold on
plot(1:counter,quantile(ISIcorrShuf,0.05,1),'-','color',[0.5 0.5 0.5])
plot(1:counter,quantile(ISIcorrShuf,0.95,1),'-','color',[0.5 0.5 0.5])
scatter(1:counter,nanmean(ISIcorr,1),50,[0 0 0],'filled')
title('ISI')
makepretty
offsetAxes
ylabel('Corr')

%%
figure('name','Functional Scores')

subplot(3,1,1)
hold on; 
h=violinplot(tmplcorr); 
title('Lateral stimulus')
makepretty
offsetAxes
ylabel('Corr')

subplot(3,1,2)
hold on; 
h=violinplot(tmpccorr); 
title('Central stimulus')
makepretty
offsetAxes
ylabel('Corr')


subplot(3,1,3)
hold on; 
h=violinplot(ISIcorr); 
title('ISI')
makepretty
offsetAxes
ylabel('Corr')
linkaxes

T = cat(2,tmplcorr,tmpccorr,ISIcorr);
% repeated measures design
T = array2table(T);
withinDesign = table([repmat(0,1,nRec-1) repmat(1,1,nRec-1) repmat(2,1,nRec-1)]',[1:nRec-1 1:nRec-1 1:nRec-1]','VariableNames',{'Stimulus','Day'});
withinDesign.Day = categorical(withinDesign.Day);
withinDesign.Stimulus = categorical(withinDesign.Stimulus);
% Repeated measures model
rm = fitrm(T,['T1-T' num2str(size(T,2)) ' ~ 1'],'WithinDesign',withinDesign);

AT = ranova(rm,'WithinModel','Day*Stimulus');
multcompare(rm,'Stimulus')
multcompare(rm,'Day')


%% Ratio
tmpl = squeeze(nanmean(dRR(:,:,1,timeline>0&timeline<0.5),4)); % tmp l
tmpc = squeeze(nanmean(dRR(:,:,2,timeline>0&timeline<0.5),4)); % tmp c

RatioCL = squeeze(tmpc-tmpl)./(tmpc+tmpl);
scatter(RelevantDays,nanmean(RatioCL(sortidx,:),1),35,ColOpt,'filled')
   


RatioChange = diff(RatioCL(sortidx,:)');
RatioChange(abs(RatioChange)<0.2) = 0;
figure
subplot(2,2,[1,2])
hold on
for exid = 1:nExample
    plot(RelevantDays,RatioCL(ExampleID(exid),:),'.-','color',NeuronCols(ExIds(exid),:))
end
line([RelevantDays(1) RelevantDays(end)],[0 0],'Color',[0.2 0.2 0.2])
xlim([RelevantDays(1)-0.5 RelevantDays(end)+0.5])
xlabel('Day')
ylabel('ratio')
makepretty
offsetAxes


subplot(2,2,3)
piechart(categorical((sign(RatioChange(1,:)))))

subplot(2,2,4)
piechart(categorical((sign(RatioChange(2,:)))))


%%
figure
hold on
for exid = 1:nExample
    plot(RelevantDays,RatioCL(ExampleID(exid),:),'.-','color',NeuronCols(ExIds(exid),:))
end
line([RelevantDays(1) RelevantDays(end)],[0 0],'Color',[0.2 0.2 0.2])
xlim([RelevantDays(1)-0.5 RelevantDays(end)+0.5])
xlabel('Day')
ylabel('ratio')
makepretty
offsetAxes

%%
figure('Name','repeated vs non-repeated');
T = cat(2,tmpc(sortidx,:),tmpl(sortidx,:));
Stim = cat(2,zeros(numel(sortidx),size(tmpc,2)),ones(numel(sortidx),size(tmpl,2)));
Day = cat(2,repmat(1:nRec,numel(sortidx),1),repmat(1:nRec,numel(sortidx),1));
Factors = table(Stim(:),Day(:),'VariableNames',{'Stimulus','Day'});
aov = anova(Factors,T(:),'Model','Interactions');
display(['non-repeated, interaction p=' num2str(round(aov.stats.pValue(3)*100)./100)])
display(['non-repeated, day p=' num2str(round(aov.stats.pValue(2)*100)./100)])
display(['non-repeated, Stimulus p=' num2str(round(aov.stats.pValue(1)*100)./100)])

subplot(2,2,1)

violinplot(tmpc(sortidx,:));
xlabel('Day')
ylabel('dRR')
title('Central')
makepretty
offsetAxes


subplot(2,2,2)
violinplot(tmpl(sortidx,:));
% errorbar(1,nanmean(nanmean(RatioCL(sortidx,3:5),2)),nanstd(nanmean(RatioCL(sortidx,3:5),2))./sqrt(numel(DiffPhases)-1))
% hold on
% errorbar(2,nanmean(nanmean(RatioCL(sortidx,6:8),2)),nanstd(nanmean(RatioCL(sortidx,6:8),2))./sqrt(numel(DiffPhases)-1))
xlabel('Day')
ylabel('dRR')
title('Lateral')
makepretty
offsetAxes
linkaxes

% repeated measures design
T = array2table(T);
withinDesign = table([repmat(0,1,nRec) repmat(1,1,nRec)]',[1:nRec 1:nRec]','VariableNames',{'Stimulus','Day'});
withinDesign.Day = categorical(withinDesign.Day);
withinDesign.Stimulus = categorical(withinDesign.Stimulus);
% Repeated measures model
rm = fitrm(T,['T1-T' num2str(size(T,2)) ' ~ 1'],'WithinDesign',withinDesign);

AT = ranova(rm,'WithinModel','Day*Stimulus');
multcompare(rm,'Stimulus')
subplot(2,2,[3])
hold on
for uid = 1:numel(sortidx)
    plot(1:nRec,tmpc(sortidx(uid),:),'.-','color',NeuronCols(uid,:))
end
xlabel('Day')
ylabel('dRR')
makepretty
offsetAxes
title('Central')

clear h
subplot(2,2,[4])
hold on
for uid = 1:numel(sortidx)
     h(uid) = plot(1:nRec,tmpl(sortidx(uid),:),'.-','color',NeuronCols(uid,:));
end
xlabel('Day')
ylabel('dRR')
legend(h(:))

title('Lateral')
makepretty
offsetAxes
linkaxes
display(['repeated, interaction p=' num2str(round(AT.pValue(7)*100)./100)])
display(['repeated, day p=' num2str(round(AT.pValue(5)*100)./100)])
display(['repeated, Stimulus p=' num2str(round(AT.pValue(3)*100)./100)])


% 
% figure('name','VennDiagrams')
% 
% venn()
% % scatter(Ydata(RelevantDays),squeeze(nanmean(RatioCL(sortidx,:,timeline>0&timeline<0.5),3)),30,'filled')




