function Depth2AreaPerUnit = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,removenoise,surfacefirst,LFPDir,treeversion,trackcoordinates)
% Enny van Beest, based on AP_histology from AJPeters & IBLAPP from Mayo
% Faulkner

%% Important
% This tool is purely to align ephys-data (functional) to histology data,
% as obtained by other methods. Preferably brainglobe (https://docs.brainglobe.info/) , or
% alternatively allenccf (https://github.com/cortex-lab/allenCCF)
% histinfo needs to come at a fine enough scale to include ALL areas, it
% will not interpollate in the Allen Brain and assume additional areas.

%% Inputs:
% histology info: probe track coordinates. If histinfo is a cell, every
% cell is expected to be from a shank of a multiple-shank probe
% AP_Histology pipeline (probe_ccf) or Brain Globe output (CSV file), 100a
% equally spaced points along probe track.
% AllenCCFPath: Path to AllenCCF (Github repository)
% Output from sp = loadKSdir(myKsDir); (Nick Steinmetz: https://github.com/cortex-lab/spikes)
% cluster information (KS/Phy output): channel (ID per cluster) and depth (also per Cluster)

%% Optional inputs:
% surfacefirst = 1: position with lowest index is the surface of the brain. Default zero: Position with heighest index deeper in the brain
% LFP folder, to also align using LFP. if empty don't use this
% treeversion: which Allen Brain Tree structure to use? (default 2, = 2017; 1 = older)
% trackcoordinates: for as many datapoints as in histinfo the X/Y/Z
% coordinates of the probe (e.g. the npy file from Brain Globe Output,
% using readNPY(fullfile(histofile(1).folder,strrep(histofile(1).name,'.csv','.npy')))

%% Outputs:
% Depth2AreaPerUnit: Table with Cluster_ID,Depth,Areaname (as defined by
% Allen Brain),Color (as defined by Allen Brain), Coordinates (as defined
% by Allen Brain) for every unit. Be aware that Cluster_ID change with
% merging/splitting, so this output is no longer valid after changes are
% made with e.g. phy --> re-run alignatlasdata

%% Check inputs
if ~iscell(histinfo)
    try
    histinfo{1} = histinfo;
    catch
        histinfo = {histinfo};
    end    
end
if nargin>8 && ~iscell(trackcoordinates)
    trackcoordinates{1} = trackcoordinates;
end

%% LFP?
LFP_On = 0;
if nargin>7 && exist(LFPDir) && ~isempty(LFPDir)
    freqBands = {[1.5 4], [4 12], [12 20], [20 30], [25 100],[150 200]};
    FreqNames = {'Delta','Theta','Alpha','Beta','Gamma','Ripples'};
    
    if ~isempty(strfind(LFPDir,'.lf')) % Saved out separately for NP1
        % Get information from meta file
        lfpD = dir(LFPDir);
        [Imecmeta] = ReadMeta2(lfpD.folder,'lf');
        lfpFs = str2num(Imecmeta.imSampRate);
        nChansInFile = strsplit(Imecmeta.acqApLfSy,',');  % neuropixels phase3a, from spikeGLX
        nChansInFile = str2num(nChansInFile{1})+1; %add one for sync
        
        try
            [lfpByChannel, allPowerEst, F, allPowerVar] = ...
                lfpBandPower(fullfile(lfpD.folder,lfpD.name), lfpFs, nChansInFile, freqBands);
            allPowerEst = allPowerEst(:,1:nChansInFile)'; % now nChans x nFreq
            
            LFP_On =1;
            %normalize LFP per frequency
            lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,1))./nanstd(lfpByChannel,[],1);
                %normalize LFP per channel
            lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,2))./nanstd(lfpByChannel,[],2);
        catch ME
            disp(ME)
            LFP_On =0;
            
        end
          
    elseif ~isempty(strfind(LFPDir,'.ap')) % Saved out in .ap for NP1
        [lfpByChannel, allPowerEst, F, allPowerVar] = lfpBandPowerNP2(LFPDir,freqBands);
        LFP_On=1;
        %normalize LFP per frequency
        lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,1))./nanstd(lfpByChannel,[],1);
        lfpByChannel = (lfpByChannel-nanmean(lfpByChannel,2))./nanstd(lfpByChannel,[],2);                %normalize LFP per channel

    else
        disp('No file found...')
        LFP_On=0;
    end
    clear allPowerEst allPowerVar
end

%% Extract all fields in sp
field=fieldnames(sp);
fnidx = find(ismember(field,{'clu','st','spikeDepths','RecSes'}));
for fn=fnidx'
    eval([field{fn} '= extractfield(sp,field{fn});'])
end
spikeCluster = clu;
spikeTimes = st;
spikeRecording = RecSes;


%% Extract cluster info
try
    cluster_id = clusinfo.id;
catch
    cluster_id = clusinfo.cluster_id; %No manual curation?
end
try
    Label = clusinfo.group;
catch
    Label = clusinfo.KSLabel; %No manual curation?
end
try
    ShankID = clusinfo.ShankID;
catch
    ShankID = ones(1,length(Label));
end
nshanks = length(histinfo);
spikeShank = nan(length(spikeCluster),1);
for shid = 1:nshanks
    spikeShank(ismember(spikeCluster,cluster_id(ShankID==shid)) & ismember(RecSes,clusinfo.RecSesID(ShankID==shid))) = shid;
end
depth = clusinfo.depth;
channel = clusinfo.ch;
RecSes_ID = clusinfo.RecSesID;
Good_ID = find(ismember(cellstr(Label),'good')); %Identify good clusters
if isempty(Good_ID)
    Good_ID = find(ismember(Label,'g')); %Identify good clusters
end
try

    Noise_ID = find(clusinfo.Noise_ID);
catch
    Noise_ID = find(~ismember(Label,'g'));
end
%% Select only good units to clean up MUA
if nargin>4 && removenoise
    spikeID = ~ismember(spikeCluster,Noise_ID-1); %0-indexed
else
    spikeID = true(length(spikeCluster),1);
    try
        % Some form of quality control
        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes<=quantile(clusinfo.n_spikes,0.01))))=0; %too little spikes
        spikeID(ismember(spikeCluster,cluster_id(clusinfo.n_spikes>=quantile(clusinfo.n_spikes,0.99))))=0; %too many spikes
        spikeID(ismember(spikeCluster,cluster_id(clusinfo.amp>=quantile(clusinfo.amp,0.99))))=0; %ridiculous amplitudes
        spikeID(ismember(spikeCluster,cluster_id(clusinfo.fr>=quantile(clusinfo.fr,0.99))))=0; %ridiculous firing rates
        spikeID(ismember(spikeCluster,cluster_id(clusinfo.ContamPct>=quantile(clusinfo.ContamPct,0.99))))=0; %ridiculous Contamination percentage
        spikeID(ismember(spikeCluster,cluster_id(ismember(cellstr(Label),'noise'))))=0; %noise should not count, only MUA and good unit
    catch
        disp('This is non curated data, using only good units from kilosort output')
%         spikeID = ismember(spikeCluster,Good_ID);
    end
end
if size(spikeID,2)==1
    spikeID=spikeID';
end
%% Surface first?
if nargin<6
    surfacefirst = 0;
end

%% Load structure tree allen brain
if nargin<8
    treeversion = 2;
end
if treeversion == 1
    tmp = readtable(fullfile(AllenCCFPath,'structure_tree_safe.csv'));
elseif treeversion == 2
    tmp = readtable(fullfile(AllenCCFPath,'structure_tree_safe_2017.csv'));
end
acronyms = lower(tmp.acronym);
color_hex = tmp.color_hex_triplet;

%% Actual coordinates known? - coordinates of track in Allen Brain space
coordinateflag =0;
if nargin>8
    coordinateflag = 1;
    figure
    % Sometimes trackcoordinates are saved in a weird order by

    
    cols = lines(length(trackcoordinates));
    DistProbe = nan(1,length(trackcoordinates));
    for trid = 1:length(trackcoordinates)
        % brainglobe. Try to correct that here:
        [~,maxid]=max(nanvar(trackcoordinates{trid},[],1)); %Sort in the dimension with largest variance
        [~,sortidx]=sort(trackcoordinates{trid}(:,maxid),'descend');
        
        trackcoordinates{trid} = trackcoordinates{trid}(sortidx,:);
        histinfo{trid} = histinfo{trid}(sortidx,:);
        
        X_ave=mean(trackcoordinates{trid},1);            % mean; line of best fit will pass through this point
        dX=bsxfun(@minus,trackcoordinates{trid},X_ave);  % residuals
        C=(dX'*dX)/(size(trackcoordinates{trid},1)-1);           % variance-covariance matrix of X
        [R,DProbe]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
        
        DProbe=diag(DProbe);
        R2=DProbe(1)/sum(DProbe);
        disp(['Linear fit R2 = ' num2str(round(R2*1000)/10) '%'])
        
        x=dX*R(:,1);    % project residuals on R(:,1)
        x_min=min(x);
        x_max=max(x);
        dx=x_max-x_min;
        Xa=(x_min-0.05*dx)*R(:,1)' + X_ave;
        Xb=(x_max+0.05*dx)*R(:,1)' + X_ave;
        X_end=[Xa;Xb];
        Fits{trid} = X_end;
        h(trid)=plot3(X_end(:,3),X_end(:,1),X_end(:,2),'-','color',cols(trid,:),'LineWidth',3); % best fit line;
        hold on
        plot3(trackcoordinates{trid}(:,3),trackcoordinates{trid}(:,1),trackcoordinates{trid}(:,2),'.k','MarkerSize',13)           % simulated noisy data
        
        % Distance line:
        DistProbe(trid) = norm(trackcoordinates{trid}(1,:)-trackcoordinates{trid}(end,:));
        hold on
    end
    legend([h(:)],arrayfun(@(X) ['Shank ' num2str(X)],1:length(trackcoordinates),'UniformOutput',0))
    title('These Shanks should be in the right order, if not order histinfo accordingly')
    DProbe = nanmax(DistProbe);
end
%% Open gui figure
gui_fig = figure('color','w');
flag = 0; %to keep track of finishing this loop
depthunique = unique(spikeDepths(spikeID));

%
% %Find 'gaps' of low activity:
% gaps = depthunique(find(nrspikesperdepth<thresh));
endpoint = max(depthunique);
startpoint = min(depthunique);

%% Make LFP plot
if LFP_On
    LFP_axis = subplot(3,9,[1:2,10:11,19:20]);
    
    %Infer depth per channel
    [sortedchannels,sortid] = unique(channel);
    sorteddepth = depth(sortid);
    
   
    imagesc(1:length(FreqNames),sorteddepth, lfpByChannel(sortedchannels+1,:),[-2 2])

%     imagesc(1:length(FreqNames),[0:(nChansInFile-1)]*10,lfpByChannel,[-2 2])
    xlim([1,length(FreqNames)]);
    set(gca, 'YDir', 'normal','XTick',1:length(FreqNames),'XTicklabel',FreqNames,'XTickLabelRotation',25);
    ylabel('depth on probe (µm)');
    % h = colorbar;
    % h.Label.String = 'power (dB)';
    colormap hot
    title('LFP')
    makepretty
    
end

%% Get multiunit correlation - Copied from Petersen github
n_corr_groups = 80;
depth_group_edges = linspace(startpoint,endpoint,n_corr_groups+1);
depth_group = discretize(spikeDepths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges)-1;

spike_binning = 0.5; % seconds
corr_edges = nanmin(spikeTimes(spikeID==1)):spike_binning:nanmax(spikeTimes(spikeID==1));
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

mua_corr = cell(1,nshanks);
for shid = 1:nshanks
    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
    parfor curr_depth = 1:length(unique_depths)
        binned_spikes_depth(curr_depth,:) = histcounts(spikeTimes(depth_group == unique_depths(curr_depth) & spikeID==1 & spikeShank'==shid), corr_edges);
    end
    %     % Z-score
    %     binned_spikes_depth = (binned_spikes_depth - nanmean(binned_spikes_depth(:)))./nanstd(binned_spikes_depth(:));
    
    binned_spikes_depth(:,nansum(binned_spikes_depth,1)>quantile(nansum(binned_spikes_depth,1),0.95))=0;
    mua_corr{shid} = smooth2a(corrcoef(binned_spikes_depth'),3);
end
mua_corr = cat(3,mua_corr{:});
mua_corr = reshape(mua_corr,size(mua_corr,1),[]);
limup = [quantile(mua_corr(:),0.1) quantile(mua_corr(:),0.95)];
% Plot multiunit correlation
if LFP_On
    multiunit_ax = subplot(3,9,[3:5,12:14,21:23]);
else
    multiunit_ax = subplot(3,9,[1:5,10:14,19:23]);
end
h=imagesc(1:length(depth_group_centers)*nshanks,depth_group_centers,mua_corr,limup);
caxis([0,max(mua_corr(mua_corr~=1))]); colormap(hot);
set(h,'Alphadata',~isnan(mua_corr))
set(gca,'Color',[0.5 0.5 0.5])
hold on
for shid = 1:nshanks
    line([size(mua_corr,1)*shid size(mua_corr,1)*shid],[startpoint endpoint],'color',[0.4 0.6 0],'LineStyle','--','LineWidth',3)
end
ylim([startpoint,endpoint]);
xlim([1,length(depth_group_centers)*nshanks]);
set(gca,'XTickLabel','')
DChannels = endpoint-startpoint;

if coordinateflag
    ScaleChannelsToProbe = endpoint./DProbe;
    disp(['Distance between 0 and top channel: ' num2str(round(endpoint)) 'micron'])
    disp(['Penetration into brain: ' num2str(round(DProbe)) 'micron'])
    disp(['Take ' num2str(round((ScaleChannelsToProbe)*100)) '% deepest areas along track'])
end
set(multiunit_ax,'YDir','normal');
title('MUA correlation');
xlabel(multiunit_ax,'Multiunit depth X Probe');
ylabel(multiunit_ax,'Multiunit depth');

if LFP_On
    set(LFP_axis,'ylim',[startpoint,endpoint])
end
%% Put histinfo in new shap with all shanks below each other
for shid = 1:length(histinfo)
    npoints = size(histinfo{shid},1);
    histinfo{shid}.shank = repmat(shid,npoints,1);
    if shid==1
        histinfonew = histinfo{shid};
    else
        histinfonew = cat(1,histinfonew,histinfo{shid});
    end
end
histinfo = histinfonew;
clear histinfonew
if coordinateflag
    DifPr = (DProbe-endpoint);
else
    DifPr = (startpoint-endpoint);
    
end
if DifPr<0
    DifPr=0;
end

while ~flag
    %Now divide position of probe along this track
    if istable(histinfo)&& any(histinfo.Position)
        histinfo.RegionAcronym(ismember(histinfo.RegionAcronym,'Not found in brain')| ismember(histinfo.RegionAcronym,'void')) = {'root'};
           for shid = 1:nshanks
            if ~surfacefirst
                if coordinateflag
                    areapoints{shid} = linspace(0-DifPr,endpoint,sum(histinfo.shank==shid));
                    trackcoordinates{shid} = [linspace(Fits{shid}(1,1),Fits{shid}(2,1),length( areapoints{shid}));linspace(Fits{shid}(1,2),Fits{shid}(2,2),length( areapoints{shid}));linspace(Fits{shid}(1,3),Fits{shid}(2,3),length( areapoints{shid}))]';
                else
                    areapoints{shid} = linspace(0-DifPr,endpoint,sum(histinfo.shank==shid));
                end
            else
                if coordinateflag
                    areapoints{shid} = linspace(endpoint,0-DifPr,sum(histinfo.shank==shid));
                    trackcoordinates{shid} = [linspace(Fits{shid}(2,1),Fits{shid}(1,1),length( areapoints{shid}));linspace(Fits{shid}(2,2),Fits{shid}(1,2),length( areapoints{shid}));linspace(Fits{shid}(2,3),Fits{shid}(1,3),length( areapoints{shid}))]';
                else
                    areapoints{shid} = linspace(endpoint,0-DifPr,sum(histinfo.shank==shid));
                end
            end
            [UniqueAreas{shid},IA{shid},IC{shid}] = unique((histinfo.RegionAcronym(histinfo.shank==shid)),'stable');
        end
    elseif isstruct(histinfo)&& isfield(histinfo,'probe_ccf')
        Error('Not yet adapted for multiple shanks')
        if ~surfacefirst
            areapoints = (linspace(startpoint,endpoint,length(histinfo.probe_ccf.trajectory_coords)));
        else
            areapoints = (linspace(endpoint,startpoint,length(histinfo.probe_ccf.trajectory_coords)));
        end
        histinfo.probe_ccf.trajectory_acronyms = acronyms(histinfo.probe_ccf.trajectory_areas);
        [UniqueAreas,IA,IC] = unique(histinfo.probe_ccf.trajectory_acronyms,'stable');
        histinfo.RegionAcronym = histinfo.probe_ccf.trajectory_acronyms;
    else
        Error('Not yet adapted for multiple shanks')
        areapoints = nan(1,max(channel)+1);
        histinfo.RegionAcronym  = cell(1,max(channel)+1);
        for chid = 1:max(channel)+1
            eval(['areapoints(1,' num2str(chid) ')= abs(histinfo.channel_' num2str(chid-1) '.z);'])
            eval([' histinfo.RegionAcronym {1,' num2str(chid) '}= histinfo.channel_' num2str(chid-1) '.brain_region;'])
        end
        % not always in the correct order, align
        if sum(ismember([-1,1],unique(sign(diff(areapoints)))))==2
            [areapoints, sortid] = sort(areapoints,'descend');
            histinfo.RegionAcronym = histinfo.RegionAcronym(sortid);
        end
        histinfo.RegionAcronym(ismember(histinfo.RegionAcronym,'void')) = {'root'};
        
        [UniqueAreas,IA,IC] = unique(fliplr(histinfo.RegionAcronym),'stable');
    end
    UniqueAreas = cellfun(@(X) lower(X),UniqueAreas,'UniformOutput',0); % case insensitive
    switchpoints = cellfun(@(X) [1; find(diff(X)~=0)+1],IC,'UniformOutput',0); %Find points where region changes
    %     switchpoints = cellfun(@(X) [1; find(diff(X)~=0)+1; length(X)],IC,'UniformOutput',0); %Find points where region changes
    AllAreas = cell(1,nshanks);
    for shid = 1:nshanks
        AllAreas{shid} = histinfo.RegionAcronym(switchpoints{shid}+(npoints*(shid-1)));
    end
    
    if exist('probe_areas_ax')
        delete(probe_areas_ax)
    end
    % To see entire probe as reference
    Entireprobe_areas_ax  = subplot(3,9,[8,9,17,18,26,27]);
    yyaxis right
    
    oripatchobj = gobjects;
    oritextobj = gobjects;
    for shid=1:nshanks
        for i=2:length(switchpoints{shid})
            oripatchobj(shid,i-1) = patch([shid-1 shid shid shid-1],[areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(i))],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i-1))}))));
            oritextobj(shid,i-1) = text(double(shid-0.5),double(nanmean([areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i))])),UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i-1))},'HorizontalAlignment','center');
        end
        oripatchobj(shid,i) = patch([shid-1 shid shid shid-1],[areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(end)) areapoints{shid}(switchpoints{shid}(end))],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i))}))));
        oritextobj(shid,i) = text(double(shid-0.5),double(nanmean([areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(end))])),UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i))},'HorizontalAlignment','center');
    end
    title('Reference')
    ylim([-inf inf])
    
    % THe one to manipulate
    probe_areas_ax  =subplot(3,9,[6,7,15,16,24,25]);
    patchobj = gobjects;
    textobj = gobjects;
    for shid=1:nshanks
        for i=2:length(switchpoints{shid})
            patchobj(shid,i-1) = patch([shid-1 shid shid shid-1],[areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(i))],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i-1))}))));
            textobj(shid,i-1) = text(double(shid-0.5),double(nanmean([areapoints{shid}(switchpoints{shid}(i-1)) areapoints{shid}(switchpoints{shid}(i))])),UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i-1))},'HorizontalAlignment','center');
        end
        patchobj(shid,i) = patch([shid-1 shid shid shid-1],[areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(end)) areapoints{shid}(switchpoints{shid}(end))],hex2rgb(color_hex(ismember(acronyms,UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i))}))));
        textobj(shid,i) = text(double(shid-0.5),double(nanmean([areapoints{shid}(switchpoints{shid}(i)) areapoints{shid}(switchpoints{shid}(end))])),UniqueAreas{shid}{IC{shid}(switchpoints{shid}(i))},'HorizontalAlignment','center');
    end
    ylim([startpoint,endpoint]);
    
    % Instructions
    disp(['Use arrow keys to move up/down, left arrow to shrink, right arrow to enlarge'])
    disp(['Use ''a/d'' to add/delete a reference line, (first) click on the probe, (then on the MUA)'])
    disp(['Use "f" to flip probe orientation'])
    disp(['To select specific shank input number'])
    disp(['To apply input to all shanks, input 0'])
    disp(['Press "i/k" to increase/decrease stepsize'])
    disp(['Press "q" to save and quite'])
    disp(['Press "r" to reset'])
    
    f = msgbox({'arrow keys: move up/down, left arrow to shrink, right arrow to enlarge';...
        '''a/d'': add/delete a reference line, (first) click on the probe, (then on the MUA)';...
        '"f": flip probe orientation';'"[Number]": To select specific shank input number';'"0": apply input to all shanks';...
        '"i/k": increase/decrease stepsize';'"q": save and quite';'"r": reset'});
    
    title('See instructions in command window');
    
    %     title({'Probe areas','(ws keys to move, a to add ref line, d to delete ref line, f for flip probe ori, 123 for factor)','(q: save & quit, r: reset)'});
    if coordinateflag
        yyaxis right
        ylim([startpoint,endpoint]);
        %Find corresponding trackcoordinates
        tmplabel=trackcoordinates{end}(cell2mat(arrayfun(@(X) find(abs(areapoints{end}-X)==min(abs(areapoints{end}-X)),1,'first'),get(gca,'YTick'),'UniformOutput',0)),:);
        tmplabel = num2cell(tmplabel,2);
        tmplabel = cellfun(@(X) ['[' num2str(round(X(1))) ';', num2str(round(X(2))), ';', num2str(round(X(3))),']'],tmplabel,'UniformOutput',0);
        set(gca,'YTickLabel',tmplabel)
        yyaxis left
    end
    
    % Draw corresponding area lines on multi unit
    subplot(multiunit_ax)
    % Draw boundary lines at borders (and undo clipping to extend across all)
    if exist('boundary_lines')
        delete(boundary_lines)
    end
    boundary_lines = gobjects;
    for shid=1:nshanks
        for curr_boundary = 1:length(switchpoints{shid})
            boundary_lines(curr_boundary,1,shid) = line(probe_areas_ax,[shid-1 shid], ...
                repmat(areapoints{shid}(switchpoints{shid}(curr_boundary)),1,2),'color','b','linewidth',1);
            boundary_lines(curr_boundary,2,shid) = line(multiunit_ax,[size(mua_corr,1)*(shid-1) size(mua_corr,1)*shid], ...
                repmat(areapoints{shid}(switchpoints{shid}(curr_boundary)),1,2),'color','y','linewidth',1,'LineStyle','--');
        end
    end
    %% Interface
    matchedswitchpoints = cell(nshanks,1);%nan(2,length(switchpoints),nshanks);
    for shid  =1:nshanks
        matchedswitchpoints{shid} = nan(2,length(switchpoints{shid}));
        matchedswitchpoints{shid}(1,:)=areapoints{shid}(switchpoints{shid});
    end
    newswitchpoints = switchpoints;
    if coordinateflag
        newtrackcoordinates = cellfun(@(X) nan(size(X)),trackcoordinates,'UniformOutput',0);
    end
    newareapoints = areapoints; %new s
    oristartpoint = startpoint;
    oriendpoint = endpoint;
    
    okay = 0;
    key = '';
    y_change = 10;
    selectedshank = [1:nshanks];
    while ~okay
        switch key
            % Set amounts to move by with/without shift
            case 'i'
                y_change = y_change*10
            case 'k'
                y_change = y_change/10
                % up/down: move probe areas
            case 'uparrow'
                disp('Moving UP')
                for shid=selectedshank
                    if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                        matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)]) + y_change;
                    else
                        matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) + y_change;
                    end
                end
            case 'downarrow'
                disp('Moving Down')
                for shid=selectedshank
                    if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                        matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)]) - y_change;
                    else
                        matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) - y_change;
                    end
                end
            case 'rightarrow' %Enlarge
                disp('Stretching')
                if surfacefirst
                    for shid=selectedshank
                        if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])+[adddif -adddif];
                        else
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])+[adddif -adddif];
                        end
                    end
                else
                    for shid=selectedshank
                        if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])+[-adddif adddif];
                        else
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])+[-adddif adddif];
                        end
                    end
                    
                    
                end
            case 'leftarrow' %Shrink
                disp('Shrinking')
                if surfacefirst
                    for shid=selectedshank
                        if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])+[-adddif adddif];
                        else
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])+[-adddif adddif];
                        end
                    end
                else
                    for shid=selectedshank
                        if isnan(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]))
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(1,[1 size(matchedswitchpoints{shid},2)])+[adddif -adddif];
                        else
                            if shid==selectedshank(1)
                                adddif = abs(diff(matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])))/2*y_change/1000;
                            end
                            matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)]) = matchedswitchpoints{shid}(2,[1 size(matchedswitchpoints{shid},2)])+[adddif -adddif];
                        end
                    end
                end
            case 'r'
                disp('Reset')
                newswitchpoints = switchpoints;
                for shid  =1:nshanks
                    matchedswitchpoints{shid}(2,:) = nan;
                    for curr_boundary = 1:length(newswitchpoints{shid})
                        set(boundary_lines(curr_boundary,1,shid),'color','b')
                        set(boundary_lines(curr_boundary,2,shid),'color','y')
                    end
                end
            case 'a' %Add reference line
                disp('Click to add reference line on probe')
                roi1 = drawpoint;
                for shid=selectedshank
                    [~,minidx] = nanmin(abs(areapoints{shid}(newswitchpoints{shid})-roi1.Position(2)));
                    set(boundary_lines(minidx,1,shid),'color','r')
                end
                delete(roi1)
                %find closest line;
                disp('Click to add reference line on MUA')
                roi2 = drawpoint;
                for shid=selectedshank
                    set(boundary_lines(minidx,2,shid),'color','r','YData',[roi2.Position(2) roi2.Position(2)])
                    matchedswitchpoints{shid}(2,minidx) = roi2.Position(2);
                end
                delete(roi2)
            case 'd' %Remove reference line
                disp('Click to remove reference line on probe')
                roi1 = drawpoint;
                for shid=selectedshank
                    [~,minidx] = nanmin(abs(areapoints{shid}(newswitchpoints{shid})-roi1.Position(2)));
                    set(boundary_lines(minidx,1,shid),'color','b')
                end
                delete(roi1)
                parfor shid=selectedshank
                    set(boundary_lines(minidx,2,shid),'color','y','YData',[areapoints{shid}(newswitchpoints{shid}(minidx)) areapoints{shid}(newswitchpoints{shid}(minidx))])
                    matchedswitchpoints{shid}(2,minidx) = nan;
                end
                delete(roi2)
                % q: save and quit
            case 'f' %Flip orientation of probe
                disp('Flipping orientation')
                surfacefirst = abs(surfacefirst-1);
                break
            case 'q'
                disp('Happy :)')
                okay = 1;
                flag = 1;
                break
            case '1'
                selectedshank = 1
            case '2'
                selectedshank =2
            case '3'
                selectedshank = 3
            case '4'
                selectedshank = 4
            case '5'
                selectedshank =5 %etc. add more if needed
            case '0'
                selectedshank = [1:nshanks]
        end
        key = '';
        if any(selectedshank>nshanks)
            warning('Selected non-existing shank. Set to select all shanks')
            selectedshank = [1:nshanks];
        end
        for shid=1:nshanks
            newswitchpoints{shid}(newswitchpoints{shid}<1) = nan;
            newswitchpoints{shid}(newswitchpoints{shid}>length(newareapoints{shid}))=nan;
            % Update figure
            if sum(~isnan(matchedswitchpoints{shid}(2,:))) > 1
                %         match the two
                nonnanidx = find(~isnan(matchedswitchpoints{shid}(2,:)));
                if nonnanidx(1)~=1
                    %         1 to first matchedswitchpoint
                    newvals = matchedswitchpoints{shid}(2,1:nonnanidx(1));
                    oldvals = matchedswitchpoints{shid}(1,1:nonnanidx(1));
                    proportion = (oldvals-oldvals(1))./(oldvals(end)-oldvals(1)); %proportion of areas in between
                    %New switchpoints, keep proportion of in between areas the same
                    newswitchpoints{shid}(1:nonnanidx(1)) = cell2mat(arrayfun(@(X) find(abs(newareapoints{shid}-X)==nanmin(abs(newareapoints{shid}-X)),1,'first'),(newvals(end)-oldvals(1))*proportion+oldvals(1),'UniformOutput',0));
                    if coordinateflag
                        newtrackcoordinates{shid}(newswitchpoints{shid}(1):newswitchpoints{shid}(nonnanidx(1)),:)= [linspace(trackcoordinates{shid}(switchpoints{shid}(1),1),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(1)),1),length(newswitchpoints{shid}(1):newswitchpoints{shid}(nonnanidx(1)))); ...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(1),2),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(1)),2),length(newswitchpoints{shid}(1):newswitchpoints{shid}(nonnanidx(1)))); ...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(1),3),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(1)),3),length(newswitchpoints{shid}(1):newswitchpoints{shid}(nonnanidx(1))))]';
                    end
                end
                for i = 1:length(nonnanidx)-1
                    newvals = matchedswitchpoints{shid}(2,nonnanidx(i):nonnanidx(i+1));
                    oldvals = matchedswitchpoints{shid}(1,nonnanidx(i):nonnanidx(i+1));
                    proportion = (oldvals-oldvals(1))./(oldvals(end)-oldvals(1)); %proportion of areas in between
                    %New switchpoints, keep proportion of in between areas the same
                    newswitchpoints{shid}(nonnanidx(i):nonnanidx(i+1)) = cell2mat(arrayfun(@(X) find(abs(newareapoints{shid}-X)==nanmin(abs(newareapoints{shid}-X)),1,'first'),(newvals(end)-newvals(1))*proportion+newvals(1),'UniformOutput',0));
                    if coordinateflag
                        newtrackcoordinates{shid}(newswitchpoints{shid}(nonnanidx(i)):newswitchpoints{shid}(nonnanidx(i+1)),:)= [linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i)),1),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),1),length(newswitchpoints{shid}(nonnanidx(i)):newswitchpoints{shid}(nonnanidx(i+1)))); ...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i)),2),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),2),length(newswitchpoints{shid}(nonnanidx(i)):newswitchpoints{shid}(nonnanidx(i+1))));...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i)),3),trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),3),length(newswitchpoints{shid}(nonnanidx(i)):newswitchpoints{shid}(nonnanidx(i+1))))]';
                    end
                end
                %Now the bit after
                if nonnanidx(i+1)<size(matchedswitchpoints{shid},2)
                    newvals = matchedswitchpoints{shid}(2,nonnanidx(i+1):end);
                    oldvals = matchedswitchpoints{shid}(1,nonnanidx(i+1):end);
                    proportion = ((oldvals-oldvals(1))./(oldvals(end)-oldvals(1))); %proportion of areas in between
                    %New switchpoints, keep proportion of in between areas the same
                    newswitchpoints{shid}(nonnanidx(i+1):end) = cell2mat(arrayfun(@(X) find(abs(newareapoints{shid}-X)==nanmin(abs(newareapoints{shid}-X)),1,'first'),(oldvals(end)-newvals(1))*proportion+newvals(1),'UniformOutput',0));
                    if coordinateflag
                        newtrackcoordinates{shid}(newswitchpoints{shid}(nonnanidx(i+1)):length(newtrackcoordinates{shid}),:)= [linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),1),trackcoordinates{shid}(end,1),length(newswitchpoints{shid}(nonnanidx(i+1)):length(newtrackcoordinates{shid}))); ...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),2),trackcoordinates{shid}(end,2),length(newswitchpoints{shid}(nonnanidx(i+1)):length(newtrackcoordinates{shid}))); ...
                            linspace(trackcoordinates{shid}(switchpoints{shid}(nonnanidx(i+1)),3),trackcoordinates{shid}(end,3),length(newswitchpoints{shid}(nonnanidx(i+1)):length(newtrackcoordinates{shid})))]';
                    end
                end
                
                if any(diff(newswitchpoints{shid})<0)
                    disp('Something went wrong, Reset')
                    newswitchpoints{shid} = switchpoints{shid};
                    matchedswitchpoints{shid}(2,:) = nan;
                    for curr_boundary = 1:length(newswitchpoints{shid})
                        set(boundary_lines(curr_boundary,1,shid),'color','b')
                        set(boundary_lines(curr_boundary,2,shid),'color','y')
                    end
                end
                
                for i=2:length(newswitchpoints{shid})
                    patchobj(shid,i-1).YData=[newareapoints{shid}(newswitchpoints{shid}(i-1)) newareapoints{shid}(newswitchpoints{shid}(i-1)) newareapoints{shid}(newswitchpoints{shid}(i)) newareapoints{shid}(newswitchpoints{shid}(i))];
                    textobj(shid,i-1).Position(2) = nanmean([newareapoints{shid}(newswitchpoints{shid}(i-1)) newareapoints{shid}(newswitchpoints{shid}(i))]);
                end
                patchobj(shid,i).YData=[newareapoints{shid}(newswitchpoints{shid}(i)) newareapoints{shid}(newswitchpoints{shid}(i)) newareapoints{shid}(newswitchpoints{shid}(end)) newareapoints{shid}(newswitchpoints{shid}(end))];
                textobj(shid,i).Position(2) = nanmean([newareapoints{shid}(newswitchpoints{shid}(i)) newareapoints{shid}(newswitchpoints{shid}(end))]);
                
                ylim([oristartpoint oriendpoint])
                
                % update boundary lines at borders (and undo clipping to extend across all)
                for curr_boundary = 1:length(newswitchpoints{shid})
                    set(boundary_lines(curr_boundary,1,shid),'YData',repmat(newareapoints{shid}(newswitchpoints{shid}(curr_boundary)),1,2))
                    set(boundary_lines(curr_boundary,2,shid),'YData',repmat(newareapoints{shid}(newswitchpoints{shid}(curr_boundary)),1,2))
                end
                
                if coordinateflag
                    subplot(probe_areas_ax)
                    yyaxis right
                    hold on
                    
                    ylim([oristartpoint,oriendpoint]);
                    %Find corresponding trackcoordinates
                    tmplabel=newtrackcoordinates{shid}(cell2mat(arrayfun(@(X) find(abs(newareapoints{shid}-X)==min(abs(newareapoints{shid}-X)),1,'first'),get(gca,'YTick'),'UniformOutput',0)),:);
                    tmplabel = num2cell(tmplabel,2);
                    tmplabel = cellfun(@(X) ['[' num2str(round(X(1))) ';', num2str(round(X(2))), ';', num2str(round(X(3))),']'],tmplabel,'UniformOutput',0);
                    set(gca,'YTickLabel',tmplabel)
                    yyaxis left
                    
                end
            end
        end
        
        %Input?
        waitforbuttonpress
        key = get(gcf,'CurrentKey');
    end
end
try
close(f)
catch
end
%% Shift area
areasaligned = cell(nshanks,length(histinfo.RegionAcronym)/nshanks);
for shid=1:nshanks
    for i = 1:length(newswitchpoints{shid})-1
        areasaligned(shid,newswitchpoints{shid}(i):newswitchpoints{shid}(i+1)) = lower(AllAreas{shid}(i));
        %          areasaligned(shid,newswitchpoints{shid}(i):newswitchpoints{shid}(i+1)) = lower(UniqueAreas{shid}{IC{shid}(i)})
    end
end
areasaligned(cell2mat(cellfun(@isempty,areasaligned,'UniformOutput',0)))={'root'};

%% Make a depth to area conversion table - per unit
for shid=1:nshanks
    tmp = (arrayfun(@(X) find(abs(areapoints{shid}-X)==min(abs(areapoints{shid}-X)),1,'first'),depth(ShankID==shid),'UniformOutput',0));
    clusterarea = repmat({'unknown'},1,length(cluster_id(ShankID==shid)));
    idx = 1:length(clusterarea);
    idx(cell2mat(cellfun(@isempty,tmp,'UniformOutput',0))) = [];
    clusterarea(idx) = areasaligned(shid,cell2mat(tmp(idx)));
    clustercolor = repmat({'#808080'},1,length(cluster_id(ShankID==shid)));
    clustercolor(idx) = cellfun(@(X) color_hex(ismember(acronyms,X)),clusterarea(idx),'UniformOutput',0);
    
    if coordinateflag
        clustercoord = repmat({[nan,nan,nan]},1,length(cluster_id(ShankID==shid)));
        clustercoord(idx) = num2cell(newtrackcoordinates{shid}(cell2mat(tmp(idx)),:),2);
        Depth2AreaPerUnit{shid} = table(cluster_id(ShankID==shid)',depth(ShankID==shid),ShankID(ShankID==shid),clusterarea',clustercolor',clustercoord','VariableNames',{'Cluster_ID','Depth','Shank','Area','Color','Coordinates'});
    else
        Depth2AreaPerUnit{shid} = table(cluster_id(ShankID==shid)',depth(ShankID==shid),ShankID(ShankID==shid),clusterarea',clustercolor','VariableNames',{'Cluster_ID','Depth','Shank','Area','Color'});
    end
end
Depth2AreaPerUnit=cat(1,Depth2AreaPerUnit{:});

end







