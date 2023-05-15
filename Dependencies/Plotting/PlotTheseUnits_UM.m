function PlotTheseUnits_UM(Pairs,MatchTable,UniqueIDConversion,WaveformInfo,AllSessionCorrelations,param,VisibleSetting)
% Plot UM results for groups of units
if nargin<4
    VisibleSetting = 'off';
end

if ~isfield(param,'TakeChannelRadius')
        param.TakeChannelRadius = 200;
        param.waveidx = 41-7:41+15;
end
%% Extracting all relevant data/parameters
if ~isdir(fullfile(param.SaveDir,'MatchFigures'))
    mkdir(fullfile(param.SaveDir,'MatchFigures'))
end

% Load sp
disp('Loading spike information...')
ndays = length(param.KSDir);
sp = cell(1,ndays);
for did = 1:ndays
    tmp = matfile(fullfile(param.KSDir{did},'PreparedData.mat'));
    sptmp = tmp.sp;
    clear tmp

    % Only keep parameters used
    sp{did}.st = sptmp.st;
    sp{did}.spikeTemplates = sptmp.spikeTemplates;
    sp{did}.spikeAmps = sptmp.spikeAmps;
     % Replace recsesid with subsesid
    sp{did}.RecSes = repmat(did,size(sp{did}.st));
end
clear sptmp

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
clear spnew


% Calculations for ISI!
tmpst=sp.st;
maxtime = 0;
for did=1:ndays
    tmpst(sp.RecSes==did)= tmpst(sp.RecSes==did)+maxtime;
    maxtime = max(tmpst(sp.RecSes==did));
end

% Find session switch
recsesGood = UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID)); %Rec session of these units
OriClusID = UniqueIDConversion.OriginalClusID(logical(UniqueIDConversion.GoodID));
nclus = length(OriClusID);
SessionSwitch = arrayfun(@(X) find(recsesGood==X,1,'first'),1:ndays,'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch) nclus+1];

ncellsperrecording = diff(SessionSwitch);

Allchannelpos = param.channelpos;
if ~iscell(Allchannelpos)
    Allchannelpos = {Allchannelpos};
end

% Extract matchtable scores
MatchProbability = reshape(MatchTable.MatchProb,nclus,nclus);
FingerprintR = reshape(MatchTable.FingerprintCor,nclus,nclus);
RankScoreAll = reshape(MatchTable.RankScore,nclus,nclus);
for scid = 1:length(param.Scores2Include)
    eval([param.Scores2Include{scid} ' = reshape(MatchTable.' param.Scores2Include{scid} ',nclus,nclus);'])
end
%% Plot figures
timercounter = tic;
disp('Plotting pairs...')
cv=1;
for pairid=1:length(Pairs)
    tmpfig = figure('visible',VisibleSetting);

    cols =  jet(length(Pairs{pairid}));
    clear hleg
    addforamplitude=0;
    for uidx = 1:length(Pairs{pairid})
        if cv==2 %Alternate between CVs
            cv=1;
        else
            cv=2;
        end
        uid = Pairs{pairid}(uidx);
        channelpos = Allchannelpos{recsesGood(uid)};
        % Load raw data
        try
            spikeMap = readNPY(fullfile(param.KSDir{recsesGood(uid)},'RawWaveforms',['Unit' num2str(OriClusID(uid)+1) '_RawSpikes.npy'])); %0-indexed to 1-indexed
        catch
            keyboard
        end
        % Detrending
        spikeMap = permute(spikeMap,[2,1,3]); %detrend works over columns
        spikeMap = detrend(spikeMap,1); % Detrend (linearly) to be on the safe side. OVER TIME!
        spikeMap = permute(spikeMap,[2,1,3]);  % Put back in order
        %Load channels
        ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(WaveformInfo.MaxChannel(uid,cv),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
        Locs = channelpos(ChanIdx,:);
    
        subplot(3,3,[1,4])
        hold on
         scatter(Locs(:,1)*10,Locs(:,2)*10,20,[0.5 0.5 0.5],'filled') % Indicate sites
        for id = 1:length(Locs)
            plot(Locs(id,1)*10+[1:size(spikeMap,1)],Locs(id,2)*10+spikeMap(:,ChanIdx(id),cv),'-','color',cols(uidx,:),'LineWidth',1)
        end
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv)*10,WaveformInfo.ProjectedLocation(2,uid,cv)*10,20,[0 0 0],'filled') %Indicate Centroid
        hleg(uidx) = plot(WaveformInfo.ProjectedLocation(1,uid,cv)*10+[1:size(spikeMap,1)],WaveformInfo.ProjectedLocation(2,uid,cv)*10+WaveformInfo.ProjectedWaveform(:,uid,cv),'-','color',cols(uidx,:),'LineWidth',2);


        subplot(3,3,[2])
        hold on
        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled')
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv),WaveformInfo.ProjectedLocation(2,uid,cv),20,[0 0 0],'filled')

        takesamples = param.waveidx;
        takesamples = unique(takesamples(~isnan(takesamples)));
        h(1) = plot(squeeze(WaveformInfo.ProjectedLocationPerTP(1,uid,takesamples,cv)),squeeze(WaveformInfo.ProjectedLocationPerTP(2,uid,takesamples,cv)),'-','color',cols(uidx,:));
        scatter(squeeze(WaveformInfo.ProjectedLocationPerTP(1,uid,takesamples,cv)),squeeze(WaveformInfo.ProjectedLocationPerTP(2,uid,takesamples,cv)),30,takesamples,'filled')
        colormap(hot)


        subplot(3,3,5)
        if uidx==1
            plot(channelpos(:,1),channelpos(:,2),'k.')
            hold on
        end
        h(1)=plot(channelpos(WaveformInfo.MaxChannel(uid,cv),1),channelpos(WaveformInfo.MaxChannel(uid,cv),2),'.','color',cols(uidx,:),'MarkerSize',15);


        subplot(3,3,3)
        hold on
        h(1)=plot(spikeMap(:,WaveformInfo.MaxChannel(uid,cv),cv),'-','color',cols(uidx,:));

        % Scatter spikes of each unit
        subplot(3,3,6)
        hold on
        idx1=find(sp.spikeTemplates == OriClusID(uid)-1 & sp.RecSes == recsesGood(uid));
        if length(unique(Pairs{pairid}))==1
            if cv==1
                idx1 = idx1(1:floor(length(idx1)/2));
            else
                idx1 = idx1(ceil(length(idx1)/2):end);
            end
        end
        scatter(sp.st(idx1)./60,sp.spikeAmps(idx1)+addforamplitude,4,cols(uidx,:),'filled')

        xlims = get(gca,'xlim');
        % Other axis
        [h1,edges,binsz]=histcounts(sp.spikeAmps(idx1));
        %Normalize between 0 and 1
        h1 = ((h1-nanmin(h1))./(nanmax(h1)-nanmin(h1)))*10+max(sp.st./60);
        plot(h1,edges(1:end-1)+addforamplitude,'-','color',cols(uidx,:));
        addforamplitude = addforamplitude+edges(end-1);

        % compute ACG
        [ccg, t] = CCGBz([double(sp.st(idx1)); double(sp.st(idx1))], [ones(size(sp.st(idx1), 1), 1); ...
            ones(size(sp.st(idx1), 1), 1) * 2], 'binSize', param.ACGbinSize, 'duration', param.ACGduration, 'norm', 'rate'); %function
        ACG = ccg(:, 1, 1);

        subplot(3,3,7);
        hold on
        plot(t,ACG,'color',cols(uidx,:));
        title(['AutoCorrelogram'])
        makepretty
    end

    % make subplots pretty
    subplot(3,3,[1,4])
    subplot
    makepretty
    set(gca,'yticklabel',arrayfun(@(X) num2str(X./10),arrayfun(@(X) X,get(gca,'ytick')),'UniformOutput',0))
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    ylimcur = get(gca,'ylim');
    ylim([ylimcur(1) ylimcur(2)*1.005])
    xlimcur = get(gca,'xlim');
    xlim([min(Locs(:,1))*10-param.spikeWidth/2 max(Locs(:,1))*10+param.spikeWidth*1.5])
    set(gca,'xtick',get(gca,'xtick'),'xticklabel',arrayfun(@(X) num2str(X./10),cellfun(@(X) str2num(X),get(gca,'xticklabel')),'UniformOutput',0))

    legend(hleg,arrayfun(@(X) ['ID' num2str(OriClusID(X)) ', Rec' num2str(recsesGood(X))],Pairs{pairid},'Uni',0),'Location','best')
    Probs = cell2mat(arrayfun(@(X) [num2str(round(MatchProbability(Pairs{pairid}(X),Pairs{pairid}(X+1)).*100)) ','],1:length(Pairs{pairid})-1,'Uni',0));
    Probs(end)=[];
    title(['Probability=' Probs '%'])


    subplot(3,3,[2])
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    xlims = [WaveformInfo.ProjectedLocation(1,uid,cv)-25 WaveformInfo.ProjectedLocation(1,uid,cv)+25];
    ylims = [WaveformInfo.ProjectedLocation(2,uid,cv)-25 WaveformInfo.ProjectedLocation(2,uid,cv)+25];  
    set(gca,'xlim',xlims,'ylim',ylims)
    axis square
    %     legend([h(1),h(2)],{['Unit ' num2str(uid)],['Unit ' num2str(uid2)]})
    hc= colorbar;
    try
        hc.Label.String = 'timesample';
    catch ME
        disp(ME)
        keyboard
    end
    makepretty
    if exist('TrajDistSim')
        tmp = cell2mat(arrayfun(@(X) [num2str(round(TrajDistSim(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp(end)=[];
    else
        tmp = 'nan';
    end
    if exist('TrajAngleSim')
        tmp2 = cell2mat(arrayfun(@(X) [num2str(round(TrajAngleSim(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp2(end)=[];
    else
        tmp2 = 'nan';
    end
    title(['Trajectory length: ' tmp '%, angle: ' tmp2 '%'])

    subplot(3,3,5)
    xlabel('X position')
    ylabel('um from tip')
    makepretty

    if exist('CentroidDist')
        tmp = cell2mat(arrayfun(@(X) [num2str(round(CentroidDist(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp(end)=[];
    else
        tmp = 'nan';
    end
    if exist('CentroidVar')
        tmp2 = cell2mat(arrayfun(@(X) [num2str(round(CentroidVar(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp2(end)=[];
    else
        tmp2 = 'nan';
    end
    title(['Centroid Distance: ' tmp '%, Variance: ' tmp2 '%'])

    subplot(3,3,3)
    ylims = get(gca,'ylim');
    patch([param.waveidx(1) param.waveidx(end) param.waveidx(end) param.waveidx(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')

    makepretty
    if exist('WavformMSE')
        tmp = cell2mat(arrayfun(@(X) [num2str(round(WavformMSE(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp(end)=[];
    else
        tmp = 'nan';
    end
    if exist('WVCorr')
         tmp2 = cell2mat(arrayfun(@(X) [num2str(round(WVCorr(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp2(end)=[];
    else
        tmp2 = 'nan';
    end
     if exist('WavformSim')
         tmp5 = cell2mat(arrayfun(@(X) [num2str(round(WavformSim(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp5(end)=[];
    else
        tmp5 = 'nan';
    end
    if exist('AmplitudeSim')
        tmp3 = cell2mat(arrayfun(@(X) [num2str(round(AmplitudeSim(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp3(end)=[];
    else
        tmp3 = 'nan';
    end
    if exist('spatialdecaySim')
        tmp4 = cell2mat(arrayfun(@(X) [num2str(round(spatialdecaySim(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
        tmp4(end)=[];
    else
        tmp4 = 'nan';
    end
    axis square

    title(['MSE=' tmp '%, corr= ' tmp2 '%' 'SIM=, ' tmp5 '%'])

    subplot(3,3,6)
    xlabel('Time (min)')
    ylabel('Abs(Amplitude)')
    title(['Amplitude distribution'])
    ylabel('Amplitude')
    set(gca,'YTick',[])
    makepretty
    title(['Ampl=' tmp3 ', decay='  tmp4 '%'])

    subplot(3,3,8)
    idx1=find(ismember(sp.spikeTemplates, OriClusID((Pairs{pairid}))) & ismember(sp.RecSes, recsesGood(Pairs{pairid})));
    isitot = diff(sort([tmpst(idx1)]));
    histogram(isitot,'FaceColor',[0 0 0])
    hold on
    line([1.5/1000 1.5/1000],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--')
    title([num2str(round(sum(isitot*1000<1.5)./length(isitot)*1000)/10) '% ISI violations']); %The higher the worse (subtract this percentage from the Total score)
    xlabel('ISI (ms)')
    ylabel('Nr. Spikes')
    makepretty


    subplot(3,3,9)
    hold on
    for uidx = 1:length(Pairs{pairid})-1
        uid = Pairs{pairid}(uidx);
        uid2 = Pairs{pairid}(uidx+1);
        SessionCorrelations = AllSessionCorrelations{recsesGood(uid),recsesGood(uid2)};
        addthis3=-SessionSwitch(recsesGood(uid))+1;
        if recsesGood(uid2)>recsesGood(uid)
            addthis4=-SessionSwitch(recsesGood(uid2))+1+ncellsperrecording(recsesGood(uid));
        else
            addthis4=-SessionSwitch(recsesGood(uid2))+1;
        end
        plot(SessionCorrelations(uid+addthis3,:),'-','color',cols(uidx,:)); hold on; plot(SessionCorrelations(uid2+addthis4,:),'-','color',cols(uidx+1,:))

    end
    xlabel('Unit')
    ylabel('Cross-correlation')
    tmp = cell2mat(arrayfun(@(X) [num2str(round(FingerprintR(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
    tmp(end)=[];
    tmp2 = cell2mat(arrayfun(@(X) [num2str(round(RankScoreAll(Pairs{pairid}(X),Pairs{pairid}(X+1)).*10)./10) ','],1:length(Pairs{pairid})-1,'Uni',0));
    tmp2(end)=[];
    title(['Fingerprint r=' tmp ', rank=' tmp2])
    ylims = get(gca,'ylim');
    set(gca,'ylim',[ylims(1) ylims(2)*1.2])
    PosMain = get(gca,'Position');
    makepretty
    hold off

    axes('Position',[PosMain(1)+(PosMain(3)*0.8) PosMain(2)+(PosMain(4)*0.8) PosMain(3)*0.2 PosMain(4)*0.2])
    box on
    hold on
    for uidx = 1:length(Pairs{pairid})-1
        uid = Pairs{pairid}(uidx);
        uid2 = Pairs{pairid}(uidx+1);
        tmp1 = FingerprintR(uid,:);
        tmp1(uid2)=nan;
        tmp2 = FingerprintR(uid2,:);
        tmp2(uid)=nan;
        tmp = cat(2,tmp1,tmp2);
        histogram(tmp,'EdgeColor','none','FaceColor',[0.5 0.5 0.5])
        line([FingerprintR(uid,uid2) FingerprintR(uid,uid2)],get(gca,'ylim'),'color',[1 0 0])
    end
    xlabel('Finger print r')
    makepretty

    set(tmpfig,'units','normalized','outerposition',[0 0 1 1])

    fname = cell2mat(arrayfun(@(X) ['ID' num2str(OriClusID(X)) ', Rec' num2str(recsesGood(X))],Pairs{pairid},'Uni',0));
    saveas(tmpfig,fullfile(param.SaveDir,'MatchFigures',[fname '.fig']))
    saveas(tmpfig,fullfile(param.SaveDir,'MatchFigures',[fname '.bmp']))
    if strcmp(VisibleSetting,'off')
        delete(tmpfig)
    end

end
disp(['Plotting pairs took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])

return