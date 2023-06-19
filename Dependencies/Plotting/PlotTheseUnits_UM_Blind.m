function PlotTheseUnits_UM_Blind(Pairs,MatchTable,UniqueIDConversion,WaveformInfo,AllSessionCorrelations,param)
% Plot UM results for groups of units
% This time leave out functional scores and labels, save labels separately
% in table
if ~isfield(param,'TakeChannelRadius')
        param.TakeChannelRadius = 200;
        param.waveidx = 41-7:41+15;
end
%% Extracting all relevant data/parameters
if ~isfolder(fullfile(param.SaveDir,'BlindFigures'))
    mkdir(fullfile(param.SaveDir,'BlindFigures'))
end

% Load sp
disp('Loading spike information...')
nKSFiles = length(param.KSDir);
sp = cell(1,nKSFiles);
for did = 1:nKSFiles
    tmp = matfile(fullfile(param.KSDir{did},'PreparedData.mat'));
    sptmp = tmp.sp;
    clear tmp

    % Only keep parameters used
    sp{did}.st = sptmp.st;
    sp{did}.spikeTemplates = sptmp.spikeTemplates;
    sp{did}.spikeAmps = sptmp.spikeAmps;
    % Replace recsesid with subsesid
    if param.RunPyKSChronicStitched
        sp{did}.RecSes = sptmp.RecSes;
    else
        sp{did}.RecSes = repmat(did,size(sp{did}.st));
    end
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
ndays = length(param.AllRawPaths);
tmpst=sp.st;
maxtime = 0;
for did=1:ndays
    tmpst(sp.RecSes==did)= tmpst(sp.RecSes==did)+maxtime;
    maxtime = max(tmpst(sp.RecSes==did));
end


% Find session switch
if param.GoodUnitsOnly
    GoodId = logical(UniqueIDConversion.GoodID);
else
    GoodId = true(1,length(UniqueIDConversion.GoodID));
end
recsesGood = UniqueIDConversion.recsesAll(GoodId); %Rec session of these units
OriClusID = UniqueIDConversion.OriginalClusID(GoodId);
UniqueID = UniqueIDConversion.UniqueID(GoodId);
Path4Unit = UniqueIDConversion.Path4UnitNPY;
nclus = length(OriClusID);
DayOpt = unique(recsesGood);
ndays = length(DayOpt);

SessionSwitch = arrayfun(@(X) find(recsesGood==X,1,'first'),DayOpt,'Uni',0);
SessionSwitch(cellfun(@isempty,SessionSwitch))=[];
SessionSwitch = [cell2mat(SessionSwitch); nclus+1];

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

% Prepare table with information on plots
OriIDs = cell2mat(cellfun(@(X) OriClusID(X),Pairs,'Uni',0));
RecIDs = cell2mat(cellfun(@(X) recsesGood(X),Pairs,'Uni',0));
tbldat = cat(1,1:length(Pairs),OriIDs,RecIDs);
tbl = array2table(tbldat','VariableNames',{'BlindID','ClusID1','ClusID2','RecID1','RecID2'});
tbl.MatchProb = cell2mat(cellfun(@(X) MatchProbability(X(1),X(2)),Pairs,'Uni',0))';
save(fullfile(param.SaveDir,'BlindFigures','BlindTable.mat'),'tbl');

%% Flip trajectory if necessary?
% for which dimensions do we allow flipping?
channelpos_AllCat = cat(1,Allchannelpos{:});
AllowFlipping = false(size(channelpos_AllCat,2),nclus); % Dimension x channel
for uid = 1:nclus
    channelpos = Allchannelpos{recsesGood(uid)};
    %Load channels
    ChanIdx = find(cell2mat(arrayfun(@(Y) norm(channelpos(WaveformInfo.MaxChannel(uid,1),:)-channelpos(Y,:)),1:size(channelpos,1),'UniformOutput',0))<param.TakeChannelRadius); %Averaging over 10 channels helps with drift
    Locs = channelpos(ChanIdx,:);
    AllowFlipping(cell2mat(arrayfun(@(X) length(unique(Locs(:,X))),1:size(Locs,2),'Uni',0))<=2,:) = true;
end
FlipDim = find(any(AllowFlipping,2));

% Flip trajectory for flip dimensions
ProjectedLocationPerTPAllFlips = nan(size(WaveformInfo.ProjectedLocationPerTP,1),size(WaveformInfo.ProjectedLocationPerTP,2),size(WaveformInfo.ProjectedLocationPerTP,3),size(WaveformInfo.ProjectedLocationPerTP,4),length(FlipDim));
for flipid = 1:length(FlipDim)
    tmpdat = squeeze(WaveformInfo.ProjectedLocationPerTP(FlipDim(flipid),:,:,:));
    range = cat(2,nanmin(tmpdat,[],2), nanmax(tmpdat,[],2));

    % change values
    newvals = nanmin(tmpdat,[],2) + (nanmax(tmpdat,[],2) - tmpdat);
    ProjectedLocationPerTPAllFlips(:,:,:,:,flipid) = WaveformInfo.ProjectedLocationPerTP;
    ProjectedLocationPerTPAllFlips(FlipDim(flipid),:,:,:,flipid) = newvals;
end

ProjectedLocationPerTPAllFlips = cat(5,WaveformInfo.ProjectedLocationPerTP,ProjectedLocationPerTPAllFlips); % add them all together

%% Plot figures
timercounter = tic;
disp('Plotting pairs...')
cv=2;
for pairid=1:length(Pairs)
    tmpfig = figure('visible',param.VisibleSetting);
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
        try
            channelpos = Allchannelpos{recsesGood(uid)};
        catch
            channelpos = Allchannelpos{1}; % Assuming it's the same configuration
        end
        % Load raw data
        try
            spikeMap = readNPY(fullfile(Path4Unit{uid})); %0-indexed
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
         scatter(Locs(:,1)*10,Locs(:,2)*20,20,[0.5 0.5 0.5],'filled') % Indicate sites
        for id = 1:length(Locs)
            plot(Locs(id,1)*10+[1:size(spikeMap,1)],Locs(id,2)*20+spikeMap(:,ChanIdx(id),cv),'-','color',cols(uidx,:),'LineWidth',1)
        end
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv)*10,WaveformInfo.ProjectedLocation(2,uid,cv)*20,20,[0 0 0],'filled') %Indicate Centroid
        hleg(uidx) = plot(WaveformInfo.ProjectedLocation(1,uid,cv)*10+[1:size(spikeMap,1)],WaveformInfo.ProjectedLocation(2,uid,cv)*20+WaveformInfo.ProjectedWaveform(:,uid,cv),'-','color',cols(uidx,:),'LineWidth',2);


        subplot(3,3,[2])
        hold on
        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled')
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv),WaveformInfo.ProjectedLocation(2,uid,cv),20,[0 0 0],'filled')

        takesamples = param.waveidx;
        takesamples = unique(takesamples(~isnan(takesamples)));
        if uidx > 1 % To flip or not to flip?
            tmptr = squeeze(ProjectedLocationPerTPAllFlips(:,uid,takesamples,cv,:));
            [~,flipidx] = nanmin(nanmean(sqrt(nansum((tmptr-repmat(tmptmpl,1,1,size(tmptr,3))).^2,1)),2),[],3);

        else
            tmptmpl = squeeze(ProjectedLocationPerTPAllFlips(:,uid,takesamples,cv,1));
            flipidx = 1;
        end

        h(1) = plot(squeeze(ProjectedLocationPerTPAllFlips(1,uid,takesamples,cv,flipidx)),squeeze(ProjectedLocationPerTPAllFlips(2,uid,takesamples,cv,flipidx)),'-','color',cols(uidx,:));
        scatter(squeeze(ProjectedLocationPerTPAllFlips(1,uid,takesamples,cv,flipidx)),squeeze(ProjectedLocationPerTPAllFlips(2,uid,takesamples,cv,flipidx)),30,takesamples,'filled')
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
        idx1=find(sp.spikeTemplates == OriClusID(uid) & sp.RecSes == recsesGood(uid));
        if ~isempty(idx1)
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
            xlim([-0.1 0.1])
        end
        makepretty
    end
    set(tmpfig,'units','normalized','outerposition',[0 0 1 1])

    % make subplots pretty
    subplot(3,3,[1,4])
    makepretty
    set(gca,'yticklabel',arrayfun(@(X) num2str(X./20),get(gca,'ytick'),'UniformOutput',0))
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    ylimcur = get(gca,'ylim');
    ylim([ylimcur(1) ylimcur(2)*1.005])
    xlimcur = get(gca,'xlim');
    xlim([min(Locs(:,1))*10-param.spikeWidth/2 max(Locs(:,1))*10+param.spikeWidth*1.5])
    set(gca,'xtick',get(gca,'xtick'),'xticklabel',arrayfun(@(X) num2str(X./10),cellfun(@(X) str2num(X),get(gca,'xticklabel')),'UniformOutput',0))

    Probs = cell2mat(arrayfun(@(X) [num2str(round(MatchProbability(Pairs{pairid}(X),Pairs{pairid}(X+1)).*100)) ','],1:length(Pairs{pairid})-1,'Uni',0));
    Probs(end)=[];
    title(['BlindID=' num2str(pairid)])


    subplot(3,3,[2])
    xlabel('Xpos (um)')
    ylabel('Ypos (um)')
    xlims = [min(WaveformInfo.ProjectedLocation(1,Pairs{pairid},cv))-25 max(WaveformInfo.ProjectedLocation(1,Pairs{pairid},cv))+25];
    ylims = [min(WaveformInfo.ProjectedLocation(2,Pairs{pairid},cv))-25 max(WaveformInfo.ProjectedLocation(2,Pairs{pairid},cv))+25];  
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

    subplot(3,3,5)
    xlabel('X position')
    ylabel('um from tip')
    makepretty


    subplot(3,3,3)
    ylims = get(gca,'ylim');
    patch([param.waveidx(1) param.waveidx(end) param.waveidx(end) param.waveidx(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')

    makepretty

    axis square


    subplot(3,3,6)
    xlabel('Time (min)')
    ylabel('Abs(Amplitude)')
    title(['Amplitude distribution'])
    ylabel('Amplitude')
    set(gca,'YTick',[])
    makepretty

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

    tmpfp = FingerprintR(Pairs{pairid},Pairs{pairid}); % It's symmetric, take best cross-validation for illustration
    tmpfp(logical(eye(size(tmpfp)))) = nan;

    subplot(3,3,9)
    hold on
    tmp = [];
    tmp2 = [];
    for uidx = 1:length(Pairs{pairid})
        for uidx2 = 2:length(Pairs{pairid})
            if uidx2<=uidx
                continue
            end

            [r,c] = find(tmpfp == max([tmpfp(uidx,uidx2),tmpfp(uidx2,uidx)']),1,'first');

            uid = Pairs{pairid}(r);
            uid2 = Pairs{pairid}(c);
            tmp = [tmp round(FingerprintR(uid,uid2)*10)/10];
            tmp2 = [tmp2 round(RankScoreAll(uid,uid2)*10)/10];
            SessionCorrelations = AllSessionCorrelations{ismember(DayOpt,recsesGood(uid)),ismember(DayOpt,recsesGood(uid2))};
            addthis3=-SessionSwitch(ismember(DayOpt,recsesGood(uid)))+1;
            if recsesGood(uid2)>recsesGood(uid)
                addthis4=-SessionSwitch(find(ismember(DayOpt,recsesGood(uid2))))+1+ncellsperrecording(ismember(DayOpt,recsesGood(uid)));
            else
                addthis4=-SessionSwitch(find(ismember(DayOpt,recsesGood(uid2))))+1;
            end
            plot(SessionCorrelations(uid+addthis3,:),'-','color',cols(uidx,:)); hold on; plot(SessionCorrelations(uid2+addthis4,:),'-','color',cols(uidx2,:))
           
        end
    end
    xlabel('Unit')
    ylabel('Cross-correlation')
    tmp = cell2mat(arrayfun(@(X) [num2str(X) ','],tmp,'Uni',0));
    tmp(end) = [];
    tmp2 = cell2mat(arrayfun(@(X) [num2str(X) ','],tmp2,'Uni',0));
    tmp2(end) = [];

     title(['Fingerprint r=' tmp ', rank=' tmp2])
    ylims = get(gca,'ylim');
    set(gca,'ylim',[ylims(1) ylims(2)*1.2])
    PosMain = get(gca,'Position');
    makepretty
    hold off

    axes('Position',[PosMain(1)+(PosMain(3)*0.8) PosMain(2)+(PosMain(4)*0.8) PosMain(3)*0.2 PosMain(4)*0.2])
    box on
    hold on
    for uidx = 1:length(Pairs{pairid})
        for uidx2 = 2:length(Pairs{pairid})
            if uidx2<=uidx
                continue
            end
            [r,c] = find(tmpfp == max([tmpfp(uidx,uidx2),tmpfp(uidx2,uidx)']));

            uid = Pairs{pairid}(r);
            uid2 = Pairs{pairid}(c);

           
            tmp1 = FingerprintR(uid,:);
            tmp1(uid2)=nan;
            tmp2 = FingerprintR(uid2,:);
            tmp2(uid)=nan;
            tmp = cat(2,tmp1,tmp2);
            histogram(tmp,'EdgeColor','none','FaceColor',[0.5 0.5 0.5])
            line([FingerprintR(uid,uid2) FingerprintR(uid,uid2)],get(gca,'ylim'),'color',nanmean(cols([uidx,uidx2],:),1))
        end
    end
    xlabel('Finger print r')
    makepretty

    fname = ['BlindID_' num2str(pairid)];
    saveas(tmpfig,fullfile(param.SaveDir,'BlindFigures',[fname '.fig']))
    saveas(tmpfig,fullfile(param.SaveDir,'BlindFigures',[fname '.bmp']))
    if strcmp(param.VisibleSetting,'off')
        delete(tmpfig)
    end

end
disp(['Plotting pairs took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])

return