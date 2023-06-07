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

% Calculations for ISI!
ndays = length(param.AllRawPaths);

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
    
        subplot(2,3,[1,4])
        hold on
         scatter(Locs(:,1)*10,Locs(:,2)*20,20,[0.5 0.5 0.5],'filled') % Indicate sites
        for id = 1:length(Locs)
            plot(Locs(id,1)*10+[1:size(spikeMap,1)],Locs(id,2)*20+spikeMap(:,ChanIdx(id),cv),'-','color',cols(uidx,:),'LineWidth',1)
        end
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv)*10,WaveformInfo.ProjectedLocation(2,uid,cv)*20,20,[0 0 0],'filled') %Indicate Centroid
        hleg(uidx) = plot(WaveformInfo.ProjectedLocation(1,uid,cv)*10+[1:size(spikeMap,1)],WaveformInfo.ProjectedLocation(2,uid,cv)*20+WaveformInfo.ProjectedWaveform(:,uid,cv),'-','color',cols(uidx,:),'LineWidth',2);


        subplot(2,3,[2])
        hold on
        scatter(Locs(:,1),Locs(:,2),20,[0.5 0.5 0.5],'filled')
        scatter(WaveformInfo.ProjectedLocation(1,uid,cv),WaveformInfo.ProjectedLocation(2,uid,cv),20,[0 0 0],'filled')

        takesamples = param.waveidx;
        takesamples = unique(takesamples(~isnan(takesamples)));
        h(1) = plot(squeeze(WaveformInfo.ProjectedLocationPerTP(1,uid,takesamples,cv)),squeeze(WaveformInfo.ProjectedLocationPerTP(2,uid,takesamples,cv)),'-','color',cols(uidx,:));
        scatter(squeeze(WaveformInfo.ProjectedLocationPerTP(1,uid,takesamples,cv)),squeeze(WaveformInfo.ProjectedLocationPerTP(2,uid,takesamples,cv)),30,takesamples,'filled')
        colormap(hot)


        subplot(2,3,5)
        if uidx==1
            plot(channelpos(:,1),channelpos(:,2),'k.')
            hold on
        end
        h(1)=plot(channelpos(WaveformInfo.MaxChannel(uid,cv),1),channelpos(WaveformInfo.MaxChannel(uid,cv),2),'.','color',cols(uidx,:),'MarkerSize',15);


        subplot(2,3,3)
        hold on
        h(1)=plot(spikeMap(:,WaveformInfo.MaxChannel(uid,cv),cv),'-','color',cols(uidx,:));

        makepretty
    end
    set(tmpfig,'units','normalized','outerposition',[0 0 1 1])

    % make subplots pretty
    subplot(2,3,[1,4])
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


    subplot(2,3,[2])
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

    subplot(2,3,5)
    xlabel('X position')
    ylabel('um from tip')
    makepretty


    subplot(2,3,3)
    ylims = get(gca,'ylim');
    patch([param.waveidx(1) param.waveidx(end) param.waveidx(end) param.waveidx(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')

    makepretty

    axis square


    fname = ['BlindID_' num2str(pairid)];
    saveas(tmpfig,fullfile(param.SaveDir,'BlindFigures',[fname '.fig']))
    saveas(tmpfig,fullfile(param.SaveDir,'BlindFigures',[fname '.bmp']))
    if strcmp(param.VisibleSetting,'off')
        delete(tmpfig)
    end

end
disp(['Plotting pairs took ' num2str(round(toc(timercounter)./60)) ' minutes for ' num2str(nclus) ' units'])

return