function [channelPos, probeSN, recordingduration] = ChannelIMROConversion(datapath,saveout,drawthis)
%Extract actual channelPositions from metafile


% Read IMRO from meta files
if nargin==0
    [datafile datapath] = uigetfile();
end
if nargin<2
    saveout = 0;
end
if nargin<3
    drawthis = 0;
end

%Extract meta file
meta = ReadMeta2(datapath);

%Extract probe serial number
probeSN = str2num(meta.imDatPrb_sn);

% Recording duration in minutes
if isfield(meta,'fileTimeSecs')
    recordingduration = str2num(meta.fileTimeSecs)./60;
else
    recordingduration = nan;
end

% To make sure the order of recording is correct
APRecordingOrder = meta.snsChanMap;
APRecordingOrder = strsplit(APRecordingOrder,'(');
APRecordingOrder(cellfun(@isempty,APRecordingOrder)) = [];
NrChanRecorded = APRecordingOrder{1};
APRecordingOrder = APRecordingOrder(2:end);
APRecordingOrder = cellfun(@(X) strsplit(X,';'),APRecordingOrder,'UniformOutput',0);
APRecordingOrder = cellfun(@(X) X{1},APRecordingOrder,'UniformOutput',0);
APRecordingOrder = cellfun(@(X) strsplit(X,'AP'),APRecordingOrder,'UniformOutput',0);
APRecordingOrder(cell2mat(cellfun(@(X) length(X)<2,APRecordingOrder,'UniformOutput',0)))=[];
APRecordingOrder = cell2mat(cellfun(@(X) str2num(X{2}),APRecordingOrder,'UniformOutput',0));

%Extract Shank map
if isfield(meta,'snsShankMap')
    Shankmap = meta.snsShankMap;
    % Split in channels
    Shankmap = strsplit(Shankmap,'(');
    Shankmap(cellfun(@isempty,Shankmap)) = [];
    Shankmap = cellfun(@(X) strsplit(X,')'),Shankmap,'UniformOutput',0);
    Shankmap = cellfun(@(X) X{1},Shankmap,'UniformOutput',0);
    % LayOut, read number shanks, columns per shank and channels per shank
    LayOut = Shankmap{1};
    LayOut = strsplit(LayOut,',');
    LayOut = cellfun(@(X) str2num(X),LayOut,'UniformOutput',0);
    LayOut(cellfun(@isempty,LayOut))=[];
    nShanks = LayOut{1};
    nCols = LayOut{2};
    nChan = LayOut{3};
    if nChan>length(Shankmap)
        nChan = length(Shankmap)-1;
    end


    % Read channel positions - stand for shank number, column, and row
    % (channel)
    Shankmap = Shankmap(2:end);
    Shankmap = cellfun(@(X) strsplit(X,':'),Shankmap,'UniformOutput',0);
    Shank = cell2mat(cellfun(@(X) str2num(X{1}),Shankmap,'UniformOutput',0));
    Col = cell2mat(cellfun(@(X) str2num(X{2}),Shankmap,'UniformOutput',0));
    Row = cell2mat(cellfun(@(X) str2num(X{3}),Shankmap,'UniformOutput',0));
    connected = cell2mat(cellfun(@(X) str2num(X{4}),Shankmap,'UniformOutput',0));

    if nShanks==1 
        % NP 1_phase 3B
        hSep = 32;
        shankSep = 0;
        basey = 20;
        if ~any(strfind(meta.imDatPrb_pn,'4s')) % Phase 3B
            vSep = 20;      % in um
            basex = zeros(1,numel(Row));
            basex(mod(Row,2)==0) = 27;% Even staggered
            basex(mod(Row,2)==1) = 11; % odd staggered
        else
            % Npix 2 is not staggered
            basex = 27;
            vSep = 15;
        end

    elseif nShanks == 4
        % NP 2.0 MS (4 shank), probe type 24 electrode positions
        vSep = 15;      % in um
        hSep = 32;
        shankSep = 250;
        basex = 27; % shift by this much in x
        basey = 0;
    else
        error('No layout known')
    end
   
  
    % Make channelMapToPos for conversion
    channelPos = nan(length(Shank),2);
    channelPos(APRecordingOrder+1,1) = Shank*shankSep+Col*hSep+basex; %x-position
    channelPos(APRecordingOrder+1,2) = Row*vSep+basey; %y-position

elseif isfield(meta,'snsGeomMap')
    Shankmap = meta.snsGeomMap;
    % Split in channels
    Shankmap = strsplit(Shankmap,'(');
    Shankmap(cellfun(@isempty,Shankmap)) = [];
    Shankmap = cellfun(@(X) strsplit(X,')'),Shankmap,'UniformOutput',0);
    Shankmap = cellfun(@(X) X{1},Shankmap,'UniformOutput',0);
    % LayOut, read number shanks, columns per shank and channels per shank
    LayOut = Shankmap{1};
    LayOut = strsplit(LayOut,',');
    LayOut = cellfun(@(X) str2num(X),LayOut,'UniformOutput',0);
    LayOut(cellfun(@isempty,LayOut))=[];
    nShanks = LayOut{1};
    nChan = length(Shankmap)-1;
    nCols = 2;


    if nShanks==1

        % NP 1_phase 3B
        vSep = 20;      % in um
        hSep = 15;
        shankSep = 15;
        basex = 0; % shift by this much in x

    elseif nShanks == 4
        % NP 2.0 MS (4 shank), probe type 24 electrode positions
        vSep = 15;      % in um
        hSep = 32;
        shankSep = 250;
        basex = 27; % shift by this much in x

    else
        disp('No layout known')
    end

    % Read channel positions - stand for shank number, column, and row
    % (channel)
    Shankmap = Shankmap(2:end);
    Shankmap = cellfun(@(X) strsplit(X,':'),Shankmap,'UniformOutput',0);
    Shank = cell2mat(cellfun(@(X) str2num(X{1}),Shankmap,'UniformOutput',0));
    Col = cell2mat(cellfun(@(X) str2num(X{2}),Shankmap,'UniformOutput',0));
    Row = cell2mat(cellfun(@(X) str2num(X{3}),Shankmap,'UniformOutput',0));
    connected = cell2mat(cellfun(@(X) str2num(X{4}),Shankmap,'UniformOutput',0));

    % Make channelMapToPos for conversion
    channelPos = nan(length(Shank),2);
    channelPos(APRecordingOrder+1,1) = Col+(Shank*shankSep); %x-position
    channelPos(APRecordingOrder+1,2) = Row; %y-position

end


if drawthis

   
    % Draw probe layout
    xpos = [];
    ypos = [];
    for shid = 1:nShanks
        for colid = 1:nCols
            for chanid = 1:nChan
                xpos = [xpos (shid-1)*shankSep+(colid-1)*hSep+basex];
                ypos = [ypos (chanid-1)*vSep];
            end
        end
    end
    figure('name',['Probe Layout ' datapath])
    scatter(xpos,ypos,4,[0 0 0],'filled')
    hold on
    scatter(channelPos(:,1),channelPos(:,2),4,connected,'filled')

    xlim([-shankSep-0.5 nShanks*(shankSep+0.5)])
end

%% 
if saveout
    fname = (strsplit(datapath,'.ap.meta'));
    fname = fname{1};
    newName = [fname,'_kilosortChanMap.mat'];
    chanMap = (1:nChan)';
    chanMap0ind = chanMap-1;
    connected = logical(connected);
    xcoords = channelPos(:,1);   %KS2 not yet using kcoords, so x coord includes shank sep
    ycoords = channelPos(:,2); % variable names need to  match KS standard
    kcoords = Shank + 1;     %KS1 uses kcoords to force templates to be on one shank
    name = fname;
    save(newName, 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );
end
