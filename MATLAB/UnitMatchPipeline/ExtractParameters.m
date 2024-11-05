function [AllWVBParameters,param] = ExtractParameters(Path4UnitNPY,clusinfo,param)
% Prepare fitting
Interpolate = 0;

opts = optimset('Display','off');
%% Extract relevant information
nclus = length(Path4UnitNPY);
spikeWidth = param.spikeWidth;
Allchannelpos = param.AllChannelPos;
RecSes = clusinfo.RecSesID;
if isfield(param,'ImposeDrift') % To test against artifical drift
    ImposeDrift = param.ImposeDrift;
else
    ImposeDrift = 0;
end
if isfield(clusinfo,'Coordinates') && param.UseHistology% Allow for real coordinates
    for recid = 1:max(RecSes)
        % Cluster
        Coordinates = clusinfo.Coordinates(RecSes == recid);
        Depth = clusinfo.depth(RecSes == recid);
        Shank = clusinfo.Shank(RecSes == recid);

        % Probe
        tmpchan = Allchannelpos{recid};
        ypostmp = tmpchan(:,2);
        xpostmp = floor(tmpchan(:,1)./250); % Assuming no new shank if not at least 100 micron apart
        ShankOpt = unique(xpostmp);
        DepthOpt = unique(ypostmp);

        newchan = nan(size(tmpchan,1),3); % in 3D

        for shid = 1:length(ShankOpt)
            for did = 1:length(DepthOpt)
                tmpcoord = cell2mat(Coordinates(find(Depth == DepthOpt(did) & Shank == ShankOpt(shid),1,'first')));
                if isempty(tmpcoord)
                    continue
                end
                tmpchanidx = find(ypostmp == DepthOpt(did) & xpostmp == ShankOpt(shid));
                newchan(tmpchanidx,:) = repmat(tmpcoord,length(tmpchanidx),1);
            end
            tmpchanidx = find(xpostmp == ShankOpt(shid));
            newchan2 = newchan(tmpchanidx,:);

            % Interpolate in case of missing channels
            newchan2 = fillmissing(newchan2,'linear',1);
            % Fit line
            polymodel = polyfitn(cat(2,newchan2(:,3),newchan2(:,2)),newchan2(:,1),1);
            yi = polyvaln(polymodel,cat(2,newchan2(:,3),newchan2(:,2)));
%             figure; scatter3(newchan2(:,3),newchan2(:,2),newchan2(:,1),5,[0 0 0],'filled');
%             hold on
%             plot3(newchan2(:,3),newchan2(:,2),yi,'r-')
%             xlabel('X')
%             ylabel('Y')
%             zlabel('Z')
%             makepretty

            % Find parallel lines and individual channel points
            chdist = unique(abs(diff(tmpchan(tmpchanidx,1))));
            XPosOpt = unique(tmpchan(tmpchanidx,1));
            newchan3 = newchan2;
            for xid = 1:length(XPosOpt)
                polymodel2 = polymodel;
                polymodel2.Coefficients(3) = polymodel2.Coefficients(3) - chdist/2 + (chdist * (xid-1));
                yp = polyvaln(polymodel2,cat(2,newchan2(tmpchan(tmpchanidx,1)==XPosOpt(xid),3),newchan2(tmpchan(tmpchanidx,1)==XPosOpt(xid),2)));
%                 plot3(newchan2(tmpchan(tmpchanidx,1)==XPosOpt(xid),3),newchan2(tmpchan(tmpchanidx,1)==XPosOpt(xid),2),yp,'b-')
                newchan3(tmpchan(tmpchanidx,1)==XPosOpt(xid),1) = yp;
            end
%             yp = polyvaln(polymodel2,cat(2,newchan2(:,3),newchan2(:,2)));
%             vecnorm(newchan2(1,:)-newchan3(1,:))
%             scatter3(newchan3(:,3),newchan3(:,2),newchan3(:,1),5,[0 1 0],'filled');

            % Save shank
            newchan(tmpchanidx,:) = newchan3;
        
        end
        % Save all shanks
        Allchannelpos{recid} = newchan;


    end

else
    % Add 3rd dimension of probe-ids
    ProbeOpt = unique(clusinfo.ProbeID);
    RecOpt = unique(RecSes);
    for recid = 1:length(RecOpt)
        probeid = unique(clusinfo.ProbeID(clusinfo.RecSesID==RecOpt(recid)));
        if length(Allchannelpos) >= recid
            Allchannelpos{recid} = cat(2,repmat(find(ismember(ProbeOpt,probeid)),size(Allchannelpos{recid},1),1),Allchannelpos{recid});
        end
    end
end
param.Coordinates = Allchannelpos;

waveidx = param.waveidx;
NewPeakLoc = param.NewPeakLoc;

recsesAll = clusinfo.RecSesID;
if param.GoodUnitsOnly
    Good_Idx = find(clusinfo.Good_ID); %Only care about good units at this point
else
    Good_Idx = 1:length(clusinfo.Good_ID);
    disp('Use all units including MUA and noise')
end

recsesGood = recsesAll(Good_Idx);

%% Initialize
if Interpolate
    upsampling = 0.1;
else
    upsampling = 1;
end
spikeWidth_up = ((spikeWidth-1)/upsampling)+1;
NewPeakLoc_up = NewPeakLoc/upsampling;

ProjectedLocation = nan(3,nclus,2);
ProjectedLocationPerTP = nan(3,nclus,spikeWidth_up,2);
ProjectedWaveform = nan(spikeWidth_up,nclus,2); % Just take waveform on maximal channel
PeakTime = nan(nclus,2); % Peak time first versus second half
MaxChannel = nan(nclus,2); % Max channel first versus second half
waveformduration = nan(nclus,2); % Waveformduration first versus second half
Amplitude = nan(nclus,2); % Maximum (weighted) amplitude, first versus second half
spatialdecay = nan(nclus,2); % how fast does the unit decay across space, first versus second halfs
spatialdecayfit = nan(nclus,2); % Same but now exponential fit
WaveIdx = false(nclus,spikeWidth_up,2);
A0Distance = nan(nclus,2); % Distance at which amplitudes are 0
expFun = @(p,d) p(1)*exp(-p(2)*d);%+p(3); % For spatial decay
GoodnessofFit = nan(nclus,2); % To test the expFun
% expFun = @(p,d) p(1).^2.*p(2)./d.^2;%+p(3); % For spatial decay


%% Take geographically close channels (within 50 microns!), not just index!
timercounter = tic;
fprintf(1,'Extracting waveform information. Progress: %3d%%',0)
for uid = 1:nclus
    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    % load data
    spikeMap = readNPY(Path4UnitNPY{uid});

    % Detrending
    spikeMap = detrend(spikeMap,1); 

    % Interpolate
    if Interpolate
        s = size(spikeMap);
        spikeMap_up = nan(numel(1:upsampling:s(1)), s(2),s(3));
        time_up = 1:upsampling:s(1);
        for ii = 1:s(2)
            for jj = 1:s(3)
                spikeMap_up(:,ii,jj) = interp1(1:s(1),spikeMap(:,ii,jj),1:upsampling:s(1),'spline');
            end
        end
        spikeMap = spikeMap_up;
        waveidx_up = find(time_up > waveidx(1) & time_up < waveidx(end));
    else
        waveidx_up = waveidx;
    end

    tmp1 = spikeMap(:,:,1);
    tmp2 = spikeMap(:,:,2);
    if all(isnan(tmp1(:))) || all(isnan(tmp2(:)))
        continue
    end

    try
        channelpos = Allchannelpos{recsesGood(uid)};
    catch ME
        % assume they all have the same configuration
        channelpos = Allchannelpos{1};
    end
    % Figure out drift (by number of channels)
    if ImposeDrift~=0
        [val,id1,id2] = unique(diff(channelpos(:,3)));
        id1(val==0)=[];
        val(val==0) = [];
        [~,id3] = max(arrayfun(@(X) (sum(id2==X)),1:length(val)));
        Channelspacing = val(id3);

        % Apply this drift to data
        ShiftDataBy = ImposeDrift./Channelspacing;
        [UniqueXPos, id1, id2] = unique(channelpos(:,2));
        for xid = 1:length(UniqueXPos)
            tmp3 = tmp2(:,id2==xid); % All channels on this xpos
            % figure; subplot(2,2,1); imagesc(1:size(tmp3,1),channelpos(id2==xid,3),tmp3'); xlabel('time'); ylabel('Depth'); title('Original')
            % Upsample
            tmp4 = nan(size(tmp3,1),size(tmp3,2)*10);
            for tp = 1:size(tmp3,1)
                tmp4(tp,:) = interp(tmp3(tp,:),10);
            end         
            % subplot(2,2,3); imagesc(1:size(tmp4,1),channelpos(id2==xid,3),tmp4'); xlabel('time'); ylabel('Depth'); title('Upsampled')
         
            tmp4 = circshift(tmp4,round(ShiftDataBy*10),2);
            % subplot(2,2,4); imagesc(1:size(tmp4,1),channelpos(id2==xid,3),tmp4'); xlabel('time'); ylabel('Depth'); title([num2str(ImposeDrift) ' Drift added'])

            for tp = 1:size(tmp3,1)
                tmp3(tp,:) = downsample(tmp4(tp,:),10);
            end
                % subplot(2,2,2); imagesc(1:size(tmp3,1),channelpos(id2==xid,3),tmp3'); xlabel('time'); ylabel('Depth'); title([num2str(ImposeDrift) ' Drift added'])
            tmp2(:,id2==xid) = tmp3;
        end
        spikeMap(:,:,2) = tmp2;
    end

    % Extract channel positions that are relevant and extract mean location
    [~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(waveidx_up,:,:),3)),[],1));
    try
    OriChanIdx = find(cell2mat(arrayfun(@(Y) vecnorm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    catch
        keyboard;
    end
    if isempty(OriChanIdx)
        disp('Warning, could not find proper channels.. try another channel')
        OriChanIdx = find(cell2mat(arrayfun(@(Y) vecnorm(channelpos(MaxChanneltmp+1,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    end
    OriLocs = channelpos(OriChanIdx,:);

    % Extract unit parameters -
    % Cross-validate: first versus second half of session
    for cv = 1:2
        ChanIdx = OriChanIdx;
        Locs = OriLocs;
        % Find maximum channels:
        [~,MaxChannel(uid,cv)] = nanmax(nanmax(abs(spikeMap(waveidx_up,ChanIdx,cv)),[],1)); %Only over relevant channels, in case there's other spikes happening elsewhere simultaneously
        MaxChannel(uid,cv) = ChanIdx(MaxChannel(uid,cv));


        %     % Mean waveform - first extract the 'weight' for each channel, based on
        %     % how close they are to the projected location (closer = better)
        Distance2MaxChan = vecnorm(Locs-channelpos(MaxChannel(uid,cv),:),2,2);

        % Difference in amplitude from maximum amplitude
        spdctmp = abs(spikeMap(NewPeakLoc_up,ChanIdx,cv)); %(abs(spikeMap(NewPeakLoc_up,MaxChannel(uid,cv),cv))-abs(spikeMap(NewPeakLoc_up,ChanIdx,cv)))./abs(spikeMap(NewPeakLoc_up,MaxChannel(uid,cv),cv));
        % Remove zero
        spdctmp(Distance2MaxChan==0) = [];
        Distance2MaxChan(Distance2MaxChan==0) = [];

        if all(isnan(spdctmp(:)))
            continue
        end
        try
            p = lsqcurvefit(expFun,[max(spdctmp) 0.05],double(Distance2MaxChan'),(spdctmp),[],[],opts);
            GoodnessofFit(uid,cv) = nansum((spdctmp-expFun(p,double(Distance2MaxChan')))./max(spdctmp).^2);

            %             p = lsqcurvefit(expFun,[max(spdctmp) 0.05 min(spdctmp)],Distance2MaxChan',spdctmp,[],[],opts);
        catch ME
            p = [nan nan];
        end

        spatialdecayfit(uid,cv) = p(2); % The fit?
        spatialdecay(uid,cv) = nanmean(spdctmp./Distance2MaxChan'); % Or just the average?

        % Determine distance at which it's just noise
        tmpmin = (log(10)/p(2)); %use spatial decay back to 10% of origina lvalue
        if tmpmin>param.TakeChannelRadius || tmpmin<0
            tmpmin = param.TakeChannelRadius;
        end

        A0Distance(uid,cv) = tmpmin;
        ChanIdx = find(cell2mat(arrayfun(@(Y) vecnorm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))< A0Distance(uid,cv)); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);
        % Mean location:
        mu = sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1).*Locs',2)./sum(repmat(nanmax(abs(spikeMap(:,ChanIdx,cv)),[],1),size(Locs,2),1),2);
        ProjectedLocation(:,uid,cv) = mu;
        % Use this waveform - weighted average across channels:
        Distance2MaxProj = vecnorm(Locs-ProjectedLocation(:,uid,cv)',2,2);
        weight = (A0Distance(uid,cv)-Distance2MaxProj)./A0Distance(uid,cv);
        ProjectedWaveform(:,uid,cv) = nansum(spikeMap(:,ChanIdx,cv).*repmat(weight,1,size(spikeMap,1))',2)./sum(weight);
        % Find significant timepoints
        wvdurtmp = find(abs(ProjectedWaveform(:,uid,cv) - nanmean(ProjectedWaveform(1:20/upsampling,uid,cv)))>2.5*nanstd(ProjectedWaveform(1:20/upsampling,uid,cv))); % More than 2. std from baseline
        if isempty(wvdurtmp)
            wvdurtmp = waveidx_up;
        end
        wvdurtmp(~ismember(wvdurtmp,waveidx_up)) = []; %okay over achiever, gonna cut you off there
        if isempty(wvdurtmp)
            % May again be empty
            wvdurtmp = waveidx_up;
        end

        % Peak Time - to be safe take a bit of smoothing
        [~,PeakTime(uid,cv)] = nanmax(abs(ProjectedWaveform(wvdurtmp(1):wvdurtmp(end),uid,cv)));
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
        ProjectedWaveform(spikeWidth_up+lags(maxid):spikeWidth_up,uid,2) = nan;
        spikeMap(spikeWidth_up+lags(maxid):spikeWidth_up,:,2) = nan;
    end
  
    for cv = 1:2
        ChanIdx = OriChanIdx;
        Locs = channelpos(ChanIdx,:);
        % Shift data so that peak is at timepoint x
        if PeakTime(uid,1)~=NewPeakLoc_up % Yes, take the 1st CV on purpose!
            ProjectedWaveform(:,uid,cv) = circshift(ProjectedWaveform(:,uid,cv),-(PeakTime(uid,1)-NewPeakLoc_up));
            spikeMap(:,:,cv) = circshift(spikeMap(:,:,cv),-(PeakTime(uid,1)-NewPeakLoc_up),1);
            if PeakTime(uid,1)-NewPeakLoc_up<0
                ProjectedWaveform(1:-(PeakTime(uid,1)-NewPeakLoc_up),uid,cv) = nan;
                spikeMap(1:-(PeakTime(uid,1)-NewPeakLoc_up),:,cv) = nan;
            else
                ProjectedWaveform(spikeWidth_up-(PeakTime(uid,1)-NewPeakLoc_up):spikeWidth_up,uid,cv) = nan;
                spikeMap(spikeWidth_up-(PeakTime(uid,1)-NewPeakLoc_up):spikeWidth_up,:,cv) = nan;
            end
        end

        Peakval = ProjectedWaveform(NewPeakLoc_up,uid,cv);
        Amplitude(uid,cv) = Peakval;

        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))< A0Distance(uid,cv)); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);
        % Full width half maximum
        wvdurtmp = find(abs(sign(Peakval)*ProjectedWaveform(waveidx_up,uid,cv))>0.25*sign(Peakval)*Peakval);
        if ~isempty(wvdurtmp)
            wvdurtmp = [wvdurtmp(1):wvdurtmp(end)]+waveidx_up(1)-1;
            waveformduration(uid,cv) = length(wvdurtmp);
        else
            waveformduration(uid,cv) = nan;
        end

        % Mean Location per individual time point:
        tmp = cell2mat(arrayfun(@(tp) sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1).*Locs',2)./sum(repmat(abs(spikeMap(tp,ChanIdx,cv)),size(Locs,2),1),2),wvdurtmp,'Uni',0));
        % Smooth
        for dim = 1:size(tmp,1)
            tmp(dim,:) = smoothdata(tmp(dim,:),'gaussian',5);
        end
        ProjectedLocationPerTP(:,uid,wvdurtmp,cv) = tmp;
        WaveIdx(uid,wvdurtmp,cv) = 1;
        % Save spikes for these channels
        %         MultiDimMatrix(wvdurtmp,1:length(ChanIdx),uid,cv) = nanmean(spikeMap(wvdurtmp,ChanIdx,wavidx),3);

    end
    if  norm(channelpos(MaxChannel(uid,1),:)-channelpos(MaxChannel(uid,2),:))>param.TakeChannelRadius
        % keyboard
    end
end

% Downsample
PeakTime = PeakTime * upsampling;
waveformduration = waveformduration * upsampling;
ProjectedWaveform = ProjectedWaveform(1:1/upsampling:end,:,:);
ProjectedLocationPerTP = ProjectedLocationPerTP(:,:,1:1/upsampling:end,:);
WaveIdx = WaveIdx(:,1:1/upsampling:end,:);

fprintf('\n')
disp(['Extracting raw waveforms and parameters took ' num2str(toc(timercounter)) ' seconds for ' num2str(nclus) ' units'])
if nanmedian(A0Distance(:))>0.75*param.TakeChannelRadius
    disp('Warning, consider larger channel radius')
end
%% Put in struct
AllWVBParameters.ProjectedLocation = ProjectedLocation;
AllWVBParameters.ProjectedLocationPerTP = ProjectedLocationPerTP;
AllWVBParameters.ProjectedWaveform = ProjectedWaveform;
AllWVBParameters.PeakTime = PeakTime;
AllWVBParameters.MaxChannel = MaxChannel;
AllWVBParameters.waveformduration = waveformduration;
AllWVBParameters.Amplitude = Amplitude;
AllWVBParameters.spatialdecay = spatialdecay;
AllWVBParameters.spatialdecayfit = spatialdecayfit;
AllWVBParameters.WaveIdx = WaveIdx;

return
% Images for example neuron

figure; histogram(GoodnessofFit)
if 0


    figure; histogram(A0Distance(:),'FaceColor',[0 0 0],'EdgeColor',[0 0 0])
    hold on;
    line([nanmedian(A0Distance(:)) nanmedian(A0Distance(:))],get(gca,'ylim'),'color',[1 0 0])
    line([quantile(A0Distance(:),0.025) quantile(A0Distance(:),0.025)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    line([quantile(A0Distance(:),0.975) quantile(A0Distance(:),0.975)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    xlabel('d10 distance (\mum)')
    ylabel('Number units')
    makepretty


    uid = [10] % Example (AL032, take 10)
    fprintf(1,'\b\b\b\b%3.0f%%',uid/nclus*100)
    % load data
    spikeMap = readNPY(Path4UnitNPY{uid});

    % Detrending
    spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
    spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
    spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order

    %%% NOT UPSAMPLED FOR PLOTTING

    try
        channelpos = Allchannelpos{recsesGood(uid)};
    catch ME
        % assume they all have the same configuration
        channelpos = Allchannelpos{recsesGood(uid)-1};
    end

    % Extract channel positions that are relevant and extract mean location
    [~,MaxChanneltmp] = nanmax(nanmax(abs(nanmean(spikeMap(35:70,:,:),3)),[],1));
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius.*0.8); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:,:);

    % Plot
    tmp = nanmean(spikeMap(:,ChanIdx(channelpos(ChanIdx,2)==0),:),3);
    timevec = [-(param.NewPeakLoc-(1:size(spikeMap,1)))].*(1/30000)*1000; % In MS
    lims = [-150 150];
    figure;
    subplot(2,3,1)
    imagesc(timevec,channelpos(ChanIdx(channelpos(ChanIdx,2)==0),3),tmp',lims)
    colormap redblue
    makepretty
    xlabel('Time (ms)')
    ylabel('depth (\mum)')
    title('sites at 0um')
    set(gca,'ydir','normal')
    ylims = get(gca,'ylim');

    subplot(2,3,2) % Peak profile
    tmp = nanmean(spikeMap(param.NewPeakLoc,ChanIdx(channelpos(ChanIdx,2)==0),:),3);
    plot(tmp,channelpos(ChanIdx(channelpos(ChanIdx,2)==0),3),'k-')
    xlabel('\muV at peak')
    ylabel('depth (\mum)')
    ylim(ylims)
    makepretty

    subplot(2,3,3) % average waveform
    tmp = nanmean(spikeMap(:,MaxChannel(uid,1),:),3);
    plot(timevec,tmp,'k-')
    xlabel('Time (ms)')
    ylabel('\muV at peak channel')
    makepretty


    tmp = nanmean(spikeMap(:,ChanIdx(channelpos(ChanIdx,2)==32),:),3);

    subplot(2,3,4)
    imagesc(timevec,channelpos(ChanIdx(channelpos(ChanIdx,2)==32),3),tmp',lims)
    xlabel('Time (ms)')
    ylabel('Depth (\mum)')
    title('Sites at 32um')
    set(gca,'ydir','normal')
    ylim(ylims)

    colormap redblue

    makepretty
    subplot(2,3,5) % Peak profile
    tmp = nanmean(spikeMap(param.NewPeakLoc,ChanIdx(channelpos(ChanIdx,2)==32),:),3);
    plot(tmp,channelpos(ChanIdx(channelpos(ChanIdx,2)==0),3),'k-')
    xlabel('\muV at peak')
    ylabel('Depth (\mum)')
        ylim(ylims)

    makepretty

    subplot(2,3,6) % average waveform
    tmp = nanmean(spikeMap(:,MaxChannel(uid,1),:),3);
    plot(timevec,tmp,'k-')
    xlabel('Time (ms)')
    ylabel('\muV at peak channel')
    makepretty



    %         makepretty
    % Spatial decay plot
    figure('name','Spatial Decay Plot')
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(MaxChanneltmp,:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);
    Distance2MaxChan = vecnorm(Locs-channelpos(MaxChannel(uid,cv),:),2,2);
    % Difference in amplitude from maximum amplitude
    spdctmp = abs(spikeMap(NewPeakLoc,ChanIdx,cv)); %(abs(spikeMap(NewPeakLoc,MaxChannel(uid,cv),cv))-abs(spikeMap(NewPeakLoc,ChanIdx,cv)))./abs(spikeMap(NewPeakLoc,MaxChannel(uid,cv),cv));
    % Remove zero
    % spdctmp(Distance2MaxChan==0) = [];
    % Distance2MaxChan(Distance2MaxChan==0) = [];

    % Spatial decay (average oer micron)
    scatter(Distance2MaxChan,spdctmp,20,[0 0 0],'filled')
    xlabel('Distance to max channel (\mum)')
    ylabel('Amplitude (\muV)')

    hold on
    p = lsqcurvefit(expFun,[max(spdctmp) 0.05],double(Distance2MaxChan'),spdctmp,[],[],opts);
    plot(sort(Distance2MaxChan)',expFun(p,sort(Distance2MaxChan)'))
    % text(min(Distance2MaxChan),max(expFun(p,sort(Distance2MaxChan)))*0.99,['A = ' num2str(round(p(1)*100)/100) '.^2 * ' num2str(p(2)) './d^2'])
    text(min(Distance2MaxChan),max(expFun(p,sort(Distance2MaxChan)))*0.99,['A = ' num2str(round(p(1)*100)/100) '* exp(-' num2str(round(p(2)*100)/100) 'd)'])
    line([log(10)/p(2) log(10)/p(2)],get(gca,'ylim'),'color',[0.5 0.5 0.5],'LineStyle','--')

    makepretty
    %
    %     for chid = 1:length(spdctmp./Distance2MaxChan')
    %         plot([0 Distance2MaxChan(chid)],[0 spdctmp(chid)],'k--')
    %     end
    %     nanmean(spdctmp./Distance2MaxChan')

    %     subplot(1,4,4)
    %     hold on
    %     for chid = 1:length(spdctmp./Distance2MaxChan')
    %     scatter(1,spdctmp(chid)./Distance2MaxChan(chid),20,[ 0 0 0
    %         ],'filled')
    %     end
    %     hold on
    %     scatter(1,nanmean(spdctmp./Distance2MaxChan'),50,[0 0 1],'filled')
    %





    %     p = lsqcurvefit(expFun,[1 1],Distance2MaxChan',spdctmp,[],[],opts)
    %     hold on
    %     plot(sort(Distance2MaxChan)',expFun(p,sort(Distance2MaxChan)'))
    %     Error1 = nansum((expFun(p,sort(Distance2MaxChan)') - spdctmp).^2);
    %     p = polyfit(Distance2MaxChan',spdctmp,1)
    %     hold on
    %     plot(sort(Distance2MaxChan)',polyval(p,sort(Distance2MaxChan)'))
    %     Error2 = nansum(polyval(p,sort(Distance2MaxChan)' - spdctmp).^2);

end
