function channelPos = ChannelIMROConversion(datapath,drawthis)
%Extract actual channelPositions from metafile


% Read IMRO from meta files
if nargin==0
    [datafile datapath] = uigetfile();
end
%Extract meta file
meta = ReadMeta2(datapath);

%Extract Shank map
Shankmap = meta.snsShankMap;
% Split in channels
Shankmap = strsplit(Shankmap,'(');
Shankmap(cellfun(@isempty,Shankmap)) = [];
Shankmap = cellfun(@(X) strsplit(X,')'),Shankmap,'UniformOutput',0);
Shankmap = cellfun(@(X) X{1},Shankmap,'UniformOutput',0);

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


% LayOut, read number shanks, columns per shank and channels per shank
LayOut = Shankmap{1};
LayOut = strsplit(LayOut,',');
LayOut = cellfun(@(X) str2num(X),LayOut,'UniformOutput',0);
nShanks = LayOut{1};
nCols = LayOut{2};
nChan = LayOut{3};
if nChan>length(Shankmap)
    nChan = length(Shankmap)
end

if nShanks==1
    
    % NP 1_phase 3B
    vSep = 20;      % in um
    hSep = 15;
    shankSep = 0;
elseif nShanks == 4
    % NP 2.0 MS (4 shank), probe type 24 electrode positions
    vSep = 15;      % in um
    hSep = 32;
    shankSep = 250;
    
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
draw = cell2mat(cellfun(@(X) str2num(X{4}),Shankmap,'UniformOutput',0));

% Draw probe layout
xpos = [];
ypos = [];
for shid = 1:nShanks
    for colid = 1:nCols
        for chanid = 1:nChan
            xpos = [xpos (shid-1)*shankSep+(colid-1)*hSep];
            ypos = [ypos (chanid-1)*vSep];
        end
    end
end
if drawthis
figure('name',['Probe Layout ' datapath])
scatter(xpos,ypos,4,[0 0 0],'filled')
hold on
scatter(Shank*shankSep+Col*hSep,Row*vSep,4,draw,'filled')

xlim([-shankSep-0.5 nShanks*shankSep+0.5])
end
% Make channelMapToPos for conversion
channelPos = nan(length(Shank),2);
channelPos(APRecordingOrder+1,1) = Shank*shankSep+Col*hSep; %x-position
channelPos(APRecordingOrder+1,2) = Row*vSep; %y-position

