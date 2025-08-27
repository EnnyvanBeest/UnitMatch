function FigureFlick(UMDir,user,recompute, loadMATsToSave)
%%% Go through pairs to manually curate matches. Will save a .mat file
%%% with the label of the pairs (0: uncurated, 1: match, -1: non-match)
%
% Inputs:
%   pathFiles: path to where the .fig files are (mandatory)
%   user: user's name (for saving results) (default = 'default')
%   recompute: recompute the labels from scratch (default = 0)
%
% Key presses:
%   Right arrow: next pair
%   Left arrow: previous pair
%   Up arrow: label as match
%   Down arrow: label as non-match
%   o: label as I don't know (=uncurated)
%   m: go to first uncurated pair
%   p: select pair number
%   s: SAVE

mf = msgbox(["Key presses:";"Right arrow: next pair";"Left arrow: previous pair";...
    "Up arrow: label as match"; "Down arrow: label as non-match"; "o: label as I don't know (=uncurated)";...
    "m: go to first uncurated pair";"p: select pair number";"s: SAVE";"A warning will appear if you've seen all pairs"]);

if ~exist('user','var')
    warning('No user specified. Will use default computername.')
    [ret,user] = system('hostname');%'default';
    user = (user(1:end-1));
end

if ~exist('recompute','var')
    recompute = 0;
end

% Load matchtable
UMFile = dir(fullfile(UMDir,'UnitMatch.mat'));
if nargin > 3 && ~isempty(loadMATsToSave)
    TmpFile = load(fullfile(UMFile.folder,UMFile.name));
else
    TmpFile = matfile(fullfile(UMFile.folder,UMFile.name));
    loadMATsToSave = '';
end
MatchTable = TmpFile.MatchTable;

if ~any(ismember(MatchTable.Properties.VariableNames,user))
    eval(['MatchTable.' user ' = zeros(height(MatchTable),1);']) % Default is what UM output was
end


% Get list of all figures
d = dir(fullfile(UMDir,'MatchFigures','*.fig'));
% fname = cell2mat(arrayfunh(@(X) ['UID' num2str(X)],unique(UniqueID(Pairs{pairid})),'Uni',0));

%     d(cellfun(@(X) datetime(X)<datetime(UMFile.date),{d.date})) = []; % Remove figures that were made prior to this matchtable
if isempty(d)
    disp(['No figures were created after date of matchtable output: ' UMFile.date newline ...
        'Run the function DrawPairsUnitMatch to create them'])
    return
end
% filenames
fnames = {d.name};

% extract OriIDs and RecIDs for each file
file_info = cellfun(@(s) regexp(s,'OriID(\d+)Rec(\d+)','tokens'), fnames, 'uni', false);

% convert from string tokens to numeric arrays
OriIDs = cellfun(@(c) cellfun(@(x) str2double(x{1}), c), file_info, 'uni', false);
RecIDs = cellfun(@(c) cellfun(@(x) str2double(x{2}), c), file_info, 'uni', false);

% Example:
% OriIDs{i} and RecIDs{i} are numeric vectors for the i-th file

% Fill in GUI data
guiData = struct;
guiData.d = d;
guiData.name = user;
eval(['tmpMatch = MatchTable.' user ';']);
% Assumes you already have:
%  - fnames = {d.name};
%  - OriIDs, RecIDs as cell arrays per file (from earlier regex extraction)
%  - MatchTable with numeric columns: ID1, ID2, RecSes1, RecSes2

% Precompute row matrices (both orientations) for fast row-wise matching
rowsMat  = [MatchTable.ID1, MatchTable.ID2, MatchTable.RecSes1, MatchTable.RecSes2];

% Per file: find all pairwise combos and match to MatchTable rows
fileRowIdx = cell(numel(OriIDs),1);   % each cell: row indices in MatchTable (0 = not found)
filePairs  = cell(numel(OriIDs),1);   % each cell: [ID1 Rec1 ID2 Rec2] for each query row

for i = 1:numel(OriIDs)
    ids = OriIDs{i}(:);
    rs  = RecIDs{i}(:);

    if isempty(ids) || isempty(rs) || numel(ids) ~= numel(rs)
        warning('File %d: empty or mismatched IDs/Recs. Skipping.', i);
        fileRowIdx{i} = [];
        filePairs{i}  = [];
        continue
    end

    % all pairwise combinations among entries in the filename
    if numel(ids) == 1
        comb = [1 1]; % degenerate; will not match anything
    else
        comb = nchoosek(1:numel(ids), 2);
    end

    % build query as [ID1 ID2 Rec1 Rec2] in the filename’s orientation
    query = [ ids(comb(:,1)) , ids(comb(:,2)) , rs(comb(:,1)) , rs(comb(:,2)) ];
    query2 = [ ids(comb(:,2)) , ids(comb(:,1)) , rs(comb(:,2)) , rs(comb(:,1)) ];
    query = cat(1,query,query2);

    % order-invariant matching: try direct, then swapped
    [tf1, rowIdx] = ismember(query, rowsMat,  'rows');

    % store results
    fileRowIdx{i} = rowIdx;     % 0 means "no matching row found"
    filePairs{i}  = query;      % for auditing which pairs were looked up
end

guiData.match = cellfun(@(X) round(nanmean(tmpMatch(X))), fileRowIdx);
guiData.OriIDs = OriIDs;
guiData.RecIDs = RecIDs;
guiData.fileRowIdx = fileRowIdx;
guiData.filePairs = filePairs;
guiData.curr.pair = 1; % Start at first

guiData.currPairPosition = 1;
guiData.curr.updateFig = 1;
guiData.curr.match = guiData.match(guiData.curr.pair);
guiData.showFinishBox = 1;
guiData.loadMATsToSave = loadMATsToSave;

% Create figure
blindFlickGUI = figure('color','w','name',user);

guidata(blindFlickGUI, guiData);
updatePlot(blindFlickGUI);
end

function updatePlot(blindFlickGUI)
% Get guidata
guiData = guidata(blindFlickGUI);

% Slow -- to change
if guiData.curr.updateFig == 1
    % Clear previous fig
    clf(blindFlickGUI)

    % Load and copy new one
    currPairPosition = guiData.currPairPosition;
    fig = openfig(fullfile(guiData.d(currPairPosition).folder,guiData.d(currPairPosition).name),'invisible');
    axNew = findobj(fig,'type','axes');
    copyobj(axNew,blindFlickGUI)
    close(fig)
end

% Update title with label -- hacky, should be predefined (but clear
% fig prevents from doing that...)
if ~isfield(guiData,'titleMain')
    guiData.titleMain = annotation(blindFlickGUI,'textbox', [0, 1, 1, 0], 'string', 'My Text', 'EdgeColor', 'none',...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
tmpname = strrep(guiData.d(currPairPosition).name,'_',' ');
tmpname = strsplit(tmpname,'.fig');
set(guiData.titleMain, 'String', sprintf('No. %d: %s: match? %d (%d/%d left)', ...
    currPairPosition,tmpname{1},guiData.curr.match,sum(guiData.match == 0),numel(guiData.OriIDs)));

% Set functions for key presses
set(blindFlickGUI,'WindowKeyPressFcn',@keyPress,'DeleteFcn',@guiClosure);

end

function keyPress(blindFlickGUI,eventdata)
% Get guidata
guiData = guidata(blindFlickGUI);

guiData.curr.updateFig = 0;
switch eventdata.Key
    case 'rightarrow' % Next pair
        guiData.currPairPosition = guiData.currPairPosition+1;
        if guiData.currPairPosition>numel(guiData.OriIDs)
            guiData.currPairPosition = 1;
        end
        guiData.curr.pair = guiData.currPairPosition;
        guiData.curr.updateFig = 1;
    case 'leftarrow' % Previous pair
        guiData.currPairPosition = guiData.currPairPosition-1;
        if guiData.currPairPosition<1
            guiData.currPairPosition = numel(guiData.OriIDs);
        end
        guiData.curr.pair = guiData.currPairPosition;
        guiData.curr.updateFig = 1;
    case 'uparrow' % It's a match
        guiData.match(guiData.currPairPosition) = 1;
        guiData.currPairPosition = guiData.currPairPosition+1;
        guiData.curr.pair = guiData.currPairPosition;
        guiData.curr.updateFig = 1;
    case 'downarrow' % It's not a match
        guiData.match(guiData.currPairPosition) = -1;
        guiData.currPairPosition = guiData.currPairPosition+1;
        guiData.curr.pair = guiData.currPairPosition;
        guiData.curr.updateFig = 1;
    case 'o' % Back to uncurated
        guiData.match(guiData.currPairPosition) = 0;
    case 'm' % Go to first uncurated
        firstUncuratedPair = find(guiData.match == 0,1);
        if ~isempty(firstUncuratedPair)
            guiData.currPairPosition = firstUncuratedPair;
            guiData.curr.pair = guiData.currPairPosition;
            guiData.curr.updateFig = 1;
        else
            guiData.showFinishBox = 1;
        end
    case 'p' % Go to specific ID
        newPair = str2double(cell2mat(inputdlg('Go to pair:')));
        if newPair > numel(guiData.OriIDs)
            error(['Pair ' num2str(newPair) ' not present'])
        end
        guiData.curr.pair = newPair;
        guiData.currPairPosition  = newPair;
        guiData.curr.updateFig = 1;
    case 's' % save
        savedata(guiData)

end
guiData.curr.match = guiData.match(guiData.currPairPosition );


if ~any(guiData.match == 0) && guiData.showFinishBox
    msgbox('All pairs have been curated. Don''t forget to save with ''s''')
    guiData.showFinishBox = 0;
end

guidata(blindFlickGUI,guiData);
updatePlot(blindFlickGUI);
end

function savedata(guiData)

SaveDir = strsplit(guiData.d(1).folder,'MatchFigures');
SaveDir = SaveDir{1};
% Load MatchTable
load(fullfile(SaveDir,'UnitMatch.mat'))
if ~isempty(guiData.loadMATsToSave)
    MatchTable = TmpFile.MatchTable;
end

if ~any(ismember(MatchTable.Properties.VariableNames,guiData.name))
    eval(['MatchTable.' guiData.name  ' = zeros(height(MatchTable),1);'])
end

match = guiData.match;
for pairid = 1:numel(guiData.match)
    eval(['MatchTable.' guiData.name '(guiData.fileRowIdx{pairid}) = match(pairid);']);
end

disp(['Saving curation to ' fullfile(SaveDir,'UnitMatch.mat')])
save(fullfile(SaveDir,'UnitMatch.mat'),'MatchTable','-append')


end

function guiClosure(blindFlickGUI,eventdata)
answer = questdlg('Saving curation in matchtable?','Save?','Yes','No','Yes');
if strcmpi(answer,'yes')
    % Get guidata
    guiData = guidata(blindFlickGUI);

    savedata(guiData)
end


end