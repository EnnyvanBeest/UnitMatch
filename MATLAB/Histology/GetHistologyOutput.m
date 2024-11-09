histodone=0;
histoflag = 0;
Depth2Area = [];
AllenCCFPath = fullfile(GithubDir,'allenCCF');
if ~exist(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
    mkdir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe))
end
tmphistfile = dir(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,'*HistoEphysAlignment.mat'));
if ~isempty(tmphistfile) && ~NewHistologyNeeded
    tmpfile = load(fullfile(tmphistfile.folder,tmphistfile.name));
    try
        Depth2Area = tmpfile.Depth2Area;
        histodone=1;
        disp('Data already aligned to histology')
        histoflag=1;
    catch ME
        histodone=0;
    end
end
if ~histodone %Just in case it's not done yet
    if iscell(myKsDir)
        myKsDir = myKsDir{1};
    end
    histofile = dir(fullfile(myKsDir,'channel_locations.json'));
    if isempty(histofile)
        histofile = dir(fullfile(myKsDir,'*.csv'));
        if isempty(histofile)
            % Try the histology folder
            histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',thisdate,thisprobe,'*.csv'));
        end
        if isempty(histofile)
            % Maybe chronic?
            MetaFileDir = dir(fullfile(myLFDir,'**',['*' strrep(thisprobe,'Probe','imec')]));
            ImecMeta =ReadMeta2(fullfile(MetaFileDir(1).folder,MetaFileDir(1).name));
            ProbeSN = ImecMeta.imDatPrb_sn;

            histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',['*' ProbeSN],'*.csv')); % Use Probe Serial Number is safest!
            if isempty(histofile)
                histofile = dir(fullfile(HistoFolder,MiceOpt{midx},'ProbeTracks',thisprobe,'*.csv'));
            end
        end
        if length(histofile)>1
            disp('Detecting multiple files... npix2 probe?')
        end
        if isempty(histofile)
            histofile = dir(fullfile(myKsDir,'probe_ccf.mat'));
            if isempty(histofile)
                disp('No channel locations known yet, do histology first')
                histoflag=0;
            else
                disp('This doesn''t work yet')
                keyboard
                % alignment with petersen probe output
                histoflag = 1;
                histinfo =load(fullfile(histofile(1).folder,histofile(1).name));
                fullfile(fullfile(histofile(1).folder,histofile(1).name))

                % Align ephys data with probe
                Depth2Area  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,1,1,2);

            end
        else % Automatic alignment with brain globe output
            histoflag = 1;
            histinfo = arrayfun(@(X) readtable(fullfile(histofile(X).folder,histofile(X).name),'ReadVariableNames',1,'Delimiter',','),1:length(histofile),'UniformOutput',0);

            % If available; find actual depth in brain
            % corresponding to this:
            trackcoordinates = arrayfun(@(X) readNPY(fullfile(histofile(X).folder,strrep(histofile(X).name,'.csv','.npy'))),1:length(histofile),'UniformOutput',0);
            % Align ephys data with probe
            if exist('lfpD') && ~isempty(lfpD)
                Depth2Area  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,0,0,fullfile(lfpD.folder,lfpD.name),2,trackcoordinates);
            else
                Depth2Area  = alignatlasdata(histinfo,AllenCCFPath,sp,clusinfo,0,0,[],2,trackcoordinates);
            end
        end
    else
        histoflag=1;
        histinfo = fileread(fullfile(histofile(1).folder,histofile(1).name)); %Read json file
        histinfo = jsondecode(histinfo);% Decode json text

        % Make new type table of this
        AllCh = fieldnames(histinfo);
        histinfonew = cell(length(AllCh)-1,4);
        for id = 1:length(AllCh)
            if strcmp(AllCh{id},'origin')
                continue
            end
            eval(['histinfonew(id,1) = {' num2str(id) '};'])
            eval(['histinfonew(id,2) = {histinfo.' AllCh{id} '.brain_region_id};'])
            eval(['histinfonew(id,3) = {histinfo.' AllCh{id} '.brain_region};'])
            eval(['histinfonew(id,4) = {histinfo.' AllCh{id} '.brain_region};'])
        end
        histinfonew = cell2table(histinfonew);
        histinfonew.Properties.VariableNames = {'Position','RegionID','RegionAcronym','RegionName'};

        % Align ephys data with probe
        Depth2Area  = alignatlasdata(histinfonew,AllenCCFPath,sp,clusinfo,0,0,fullfile(lfpD.folder,lfpD.name),2);
    end
end
if histoflag && ~histodone
    saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.fig']))
    save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisprobe,[num2str(SN) '_HistoEphysAlignment.mat']),'Depth2Area')
end

% Add information to clusinfo
if exist('clusinfo','var') & ~isempty(Depth2Area)
    [mindist,DistIdx] = nanmin(abs(clusinfo.depth - Depth2Area.Depth')+abs(clusinfo.Shank-Depth2Area.Shank')*1000,[],2);
    clusinfo.Area = Depth2Area.Area(DistIdx);
    clusinfo.Color = cell2mat(arrayfun(@(X) hex2rgb(Depth2Area.Color{X}),DistIdx,'Uni',0));
    clusinfo.Coordinates = cell2mat(Depth2Area.Coordinates(DistIdx));


    %% I only care about certain levels of histology - not layer specificity at this point
    % Group some areas
    Area= strrep(clusinfo.Area,'/','');

    areaopt = unique(Area,'stable');
    uniquearean = length(areaopt);
    % Use Allen Brain Atlas to find all areas in the recordings
    % Add units to area specific cell
    % Area = Area(~cellfun(@isempty,Area));
    atlastable = readtable('structure_tree_safe_2017.csv');
    AutomaticAREASOfInterest = {};
    AutomaticAREASOfInterestFullName = {};
    ccfIdx = {};
    ColPerFullArea = nan(0,3);
    automaticareasofinterestid = nan(1,uniquearean); %Index for which area
    newareaname = cell(1,length(Area));
    newareaabrev = cell(1,length(Area));
    % Make names in table similar to those used here
    atlastable.acronym = lower(atlastable.acronym);
    atlastable.acronym= strrep(atlastable.acronym,'/','');
    for areaid=1:length(areaopt)
        %Find structure_id_path of area
        structure_id_path = atlastable.structure_id_path{find(ismember(lower(atlastable.acronym),areaopt{areaid}))};
        % Cut off last part
        parts = strsplit(structure_id_path,'/');
        %     parts(cellfun(@isempty,parts))=[];
        if length(parts)<=4 || any(strfind(areaopt{areaid},'ca'))%if we cannot go further up or hippocampal area
            tmpareaname = areaopt{areaid}; %if we cannot go further up
            newstructure_id_path = structure_id_path;
        else
            newstructure_id_path = fullfile(parts{2:end-2});
            newstructure_id_path = ['/' strrep(newstructure_id_path,'\','/') '/'];
            % Find area with this new structure id
            tmpareaname = atlastable.acronym{find(ismember(atlastable.structure_id_path,newstructure_id_path))};
        end
        %Does this one already exist in Automatic AREAS Of Interest?
        if ~any(ismember(AutomaticAREASOfInterest,tmpareaname))
            %no, then create
            AutomaticAREASOfInterest = {AutomaticAREASOfInterest{:} tmpareaname};
            AutomaticAREASOfInterestFullName = {AutomaticAREASOfInterestFullName{:} atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}};
            % Use average color scheme
            ColPerFullArea = cat(1,ColPerFullArea,nanmean(clusinfo.Color(ismember(Area,areaopt{areaid}),:),1));
            ccfIdx = {ccfIdx{:} find(ismember(lower(atlastable.acronym),areaopt{areaid}))};

        else
            ccfIdx{find(ismember(AutomaticAREASOfInterest,tmpareaname))} = [ccfIdx{find(ismember(AutomaticAREASOfInterest,tmpareaname))} find(ismember(lower(atlastable.acronym),areaopt{areaid}))];
        end
        % Index correctly
        automaticareasofinterestid(areaid) = find(ismember(AutomaticAREASOfInterest,tmpareaname));
        newareaabrev(ismember(Area,areaopt{areaid})) = {tmpareaname};
        newareaname(ismember(Area,areaopt{areaid}))={atlastable.name{find(ismember(atlastable.structure_id_path,newstructure_id_path))}}; %Save out per unit


    end
    if any(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'})) % make these void
        AutomaticAREASOfInterestFullName(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'}))={'root'};
        AutomaticAREASOfInterest(ismember(AutomaticAREASOfInterest,{'root','vs','fiber tracts','cc','fxs'}))={'root'};
    end
    [~,colidx] = ismember(newareaname,AutomaticAREASOfInterestFullName);
    [~,sortareaidx] = sortrows(ColPerFullArea);

    clusinfo.Area = newareaabrev';
    clusinfo.AreaFN = newareaname';

end