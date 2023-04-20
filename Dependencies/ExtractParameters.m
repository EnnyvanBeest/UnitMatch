function AllWVBParameters = ExtractParameters(Path4UnitNPY,clusinfo,param)

%% Extract relevant information
nclus = length(Path4UnitNPY);
spikeWidth = param.spikeWidth;
Allchannelpos = param.channelpos;
waveidx = param.waveidx;
NewPeakLoc = param.NewPeakLoc;

recsesAll = clusinfo.RecSesID;
Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
recsesGood = recsesAll(Good_Idx);

%% Initialize
ProjectedLocation = nan(2,nclus,2);
ProjectedLocationPerTP = nan(2,nclus,spikeWidth,2);
ProjectedWaveform = nan(spikeWidth,nclus,2); % Just take waveform on maximal channel
PeakTime = nan(nclus,2); % Peak time first versus second half
MaxChannel = nan(nclus,2); % Max channel first versus second half
waveformduration = nan(nclus,2); % Waveformduration first versus second half
Amplitude = nan(nclus,2); % Maximum (weighted) amplitude, first versus second half
spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second half
WaveIdx = false(nclus,spikeWidth,2);

%% Take geographically close channels (within 50 microns!), not just index!
timercounter = tic;
fprintf(1,'Extracting waveform information. Progress: %3d%%',0)
for uid = 1:nclus
    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    % load data
    spikeMap = readNPY(Path4UnitNPY{uid});

    % Detrending
    spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
    spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
    spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order

    try
        channelpos = Allchannelpos{recsesGood(uid)};
    catch ME
        % assume they all have the same configuration
        channelpos = Allchannelpos{recsesGood(uid)-1};
    end

    % Extract channel positions that are relevant and extract mean location
    [~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);

    % Extract unit parameters -
    % Cross-validate: first versus second half of session
    for cv = 1:2
        % Find maximum channels:
        [~,MaxChannel(uid,cv)] = nanmax(nanmax(abs(spikeMap(35:70,ChanIdx,cv)),[],1)); %Only over relevant channels, in case there's other spikes happening elsewhere simultaneously
        MaxChannel(uid,cv) = ChanIdx(MaxChannel(uid,cv));

        % Mean location:
        mu = sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1).*Locs',2)./sum(repmat(nanmax(abs(nanmean(spikeMap(:,ChanIdx,cv),3)),[],1),size(Locs,2),1),2);
        ProjectedLocation(:,uid,cv) = mu;
        % Use this waveform - weighted average across channels:
        Distance2MaxProj = sqrt(nansum(abs(Locs-ProjectedLocation(:,uid,cv)').^2,2));
        weight = (param.TakeChannelRadius-Distance2MaxProj)./param.TakeChannelRadius;
        ProjectedWaveform(:,uid,cv) = nansum(spikeMap(:,ChanIdx,cv).*repmat(weight,1,size(spikeMap,1))',2)./sum(weight);
        % Find significant timepoints
        wvdurtmp = find(abs(ProjectedWaveform(:,uid,cv))>abs(nanmean(ProjectedWaveform(1:20,uid,cv)))+2.5*nanstd(ProjectedWaveform(1:20,uid,cv))); % More than 2. std from baseline
        if isempty(wvdurtmp)
            wvdurtmp = waveidx;
        end
        wvdurtmp(~ismember(wvdurtmp,waveidx)) = []; %okay over achiever, gonna cut you off there
        if isempty(wvdurtmp)
            % May again be empty
            wvdurtmp = waveidx;
        end

        % Peak Time - to be safe take a bit of smoothing
        [~,PeakTime(uid,cv)] = nanmax(abs(smooth(ProjectedWaveform(wvdurtmp(1):wvdurtmp(end),uid,cv),2)));
        PeakTime(uid,cv) = PeakTime(uid,cv)+wvdurtmp(1)-1;
    end
    % Give each unit the best opportunity to correlate the waveform; cross
    % correlate to find the difference in peak
    [tmpcor, lags] = xcorr(ProjectedWaveform(:,uid,1),ProjectedWaveform(:,uid,2));
    [~,maxid] = max(tmpcor);
 
    % Shift accordingly
    ProjectedWaveform(:,uid,2) = circshift(ProjectedWaveform(:,uid,2),lags(maxid));
    spikeMap(:,:,2) = circshift(spikeMap(:,:,2),lags(maxid),1);
    if lags(maxid)>0
        ProjectedWaveform(1:lags(maxid),uid,2) = nan;
        spikeMap(1:lags(maxid),:,2) = nan;
    elseif lags(maxid)<0
        ProjectedWaveform(spikeWidth+lags(maxid):spikeWidth,uid,2) = nan;
        spikeMap(spikeWidth+lags(maxid):spikeWidth,:,2) = nan;
    end
    for cv = 1:2
        % Shift data so that peak is at timepoint x
        if PeakTime(uid,1)~=NewPeakLoc % Yes, take the 1st CV on purpose!
            ProjectedWaveform(:,uid,cv) = circshift(ProjectedWaveform(:,uid,cv),-(PeakTime(uid,1)-NewPeakLoc));
            spikeMap(:,:,cv) = circshift(spikeMap(:,:,cv),-(PeakTime(uid,1)-NewPeakLoc),1);
            if PeakTime(uid,1)-NewPeakLoc<0
                ProjectedWaveform(1:-(PeakTime(uid,1)-NewPeakLoc),uid,cv) = nan;
                spikeMap(1:-(PeakTime(uid,1)-NewPeakLoc),:,cv) = nan;
            else
                ProjectedWaveform(spikeWidth-(PeakTime(uid,cv)-NewPeakLoc):spikeWidth,uid,cv) = nan;
                spikeMap(spikeWidth-(PeakTime(uid,1)-NewPeakLoc):spikeWidth,:,cv) = nan;
            end
        end

        %     % Mean waveform - first extract the 'weight' for each channel, based on
        %     % how close they are to the projected location (closer = better)
        Distance2MaxChan = sqrt(nansum(abs(Locs-channelpos(MaxChannel(uid,cv),:)).^2,2));
        % Difference in amplitude from maximum amplitude
        spdctmp = (nanmax(abs(spikeMap(:,MaxChannel(uid,cv),cv)),[],1)-nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1))./nanmax(abs(spikeMap(:,MaxChannel(uid,cv),cv)),[],1);
        % Spatial decay (average oer micron)
        spatialdecay(uid,cv) = nanmean(spdctmp./Distance2MaxChan');
        Peakval = ProjectedWaveform(PeakTime(uid,cv),uid,cv);
        Amplitude(uid,cv) = Peakval;

        % Full width half maximum
        wvdurtmp = find(abs(sign(Peakval)*ProjectedWaveform(waveidx,uid,cv))>0.25*sign(Peakval)*Peakval);
        if ~isempty(wvdurtmp)
            wvdurtmp = [wvdurtmp(1):wvdurtmp(end)]+waveidx(1)-1;
            waveformduration(uid,cv) = length(wvdurtmp);
        else
            waveformduration(uid,cv) = nan;
        end

        % Mean Location per individual time point:
        ProjectedLocationPerTP(:,uid,wvdurtmp,cv) = cell2mat(arrayfun(@(tp) sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1).*Locs',2)./sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1),2),wvdurtmp,'Uni',0));
        WaveIdx(uid,wvdurtmp,cv) = 1;
        % Save spikes for these channels
        %         MultiDimMatrix(wvdurtmp,1:length(ChanIdx),uid,cv) = nanmean(spikeMap(wvdurtmp,ChanIdx,wavidx),3);

    end
    if  norm(channelpos(MaxChannel(uid,1),:)-channelpos(MaxChannel(uid,2),:))>param.TakeChannelRadius
        % keyboard
    end
end
fprintf('\n')
disp(['Extracting raw waveforms and parameters took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])

%% Put in struct
AllWVBParameters.ProjectedLocation = ProjectedLocation;
AllWVBParameters.ProjectedLocationPerTP = ProjectedLocationPerTP;
AllWVBParameters.ProjectedWaveform = ProjectedWaveform;
AllWVBParameters.PeakTime = PeakTime;
AllWVBParameters.MaxChannel = MaxChannel;
AllWVBParameters.waveformduration = waveformduration;
AllWVBParameters.Amplitude = Amplitude;
AllWVBParameters.spatialdecay = spatialdecay;
AllWVBParameters.WaveIdx = WaveIdx;
