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
        eval(['MatchTable.' user ' = zeros(height(MatchTable),1);'])
    end



    % Get list of all figures
    d = dir(fullfile(UMDir,'MatchFigures','*.fig'));
    %     d(cellfun(@(X) datetime(X)<datetime(UMFile.date),{d.date})) = []; % Remove figures that were made prior to this matchtable
    if isempty(d)
        disp(['No figures were created after date of matchtable output: ' UMFile.date newline ...
            'Run the function DrawPairsUnitMatch to create them'])
        return
    end
    UIDs = cell2mat(cellfun(@(y) str2num(y{end}(1:strfind(y{end},'_')-1)), cellfun(@(x) strsplit(x,'UID'), {d.name}, 'uni',0),'uni',0));
    [~,sortIdx] = sort(UIDs,'ascend');
    d = d(sortIdx);
    UIDs = UIDs(sortIdx);
   
    % Fill in GUI data
    guiData = struct;
    guiData.d = d;
    guiData.name = user;
    eval(['tmpMatch = MatchTable.' user ';']);
    try
    tmpMatchIdx = arrayfun(@(X)find(MatchTable.UID1 == X & MatchTable.UID2 == X,1,'first'),UIDs);
    catch ME
        disp(ME)
        disp('Maybe you didn''t yet make the images, or older images are present in this folder, remove these')
    end
    guiData.match = tmpMatch(tmpMatchIdx)';
    guiData.pairIDs = UIDs;
    guiData.curr.pair = guiData.pairIDs(1); % Start at first
    guiData.currPairPosition = 1;
    guiData.curr.updateFig = 1;
    guiData.curr.match = guiData.match(guiData.pairIDs == guiData.curr.pair);
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
    tmpname = strsplit(guiData.d(currPairPosition).name,'_');
    tmpname = tmpname{1};
    set(guiData.titleMain, 'String', sprintf('No. %d: %s: match? %d (%d/%d left)', ...
        currPairPosition,tmpname,guiData.curr.match,sum(guiData.match == 0),numel(guiData.pairIDs)));

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
            if guiData.currPairPosition>length(guiData.pairIDs)
                guiData.currPairPosition = 1;
            end
            guiData.curr.pair = guiData.pairIDs(guiData.currPairPosition);
            guiData.curr.updateFig = 1;
        case 'leftarrow' % Previous pair
            guiData.currPairPosition = guiData.currPairPosition-1;
            if guiData.currPairPosition<1
                guiData.currPairPosition = length(guiData.pairIDs);
            end
            guiData.curr.pair = guiData.pairIDs(guiData.currPairPosition);
            guiData.curr.updateFig = 1;
        case 'uparrow' % It's a match
            guiData.match(guiData.currPairPosition) = 1;
            guiData.currPairPosition = guiData.currPairPosition+1;
            guiData.curr.pair = guiData.pairIDs(guiData.currPairPosition);
            guiData.curr.updateFig = 1;
        case 'downarrow' % It's not a match
            guiData.match(guiData.currPairPosition) = -1;
            guiData.currPairPosition = guiData.currPairPosition+1;
            guiData.curr.pair = guiData.pairIDs(guiData.currPairPosition);
            guiData.curr.updateFig = 1;
        case 'o' % Back to uncurated
            guiData.match(guiData.currPairPosition) = 0;
        case 'm' % Go to first uncurated
            firstUncuratedPair = find(guiData.match == 0,1);
            if ~isempty(firstUncuratedPair)
                guiData.currPairPosition = firstUncuratedPair;
                guiData.curr.pair = guiData.pairIDs(guiData.currPairPosition);
                guiData.curr.updateFig = 1;
            else
                guiData.showFinishBox = 1;
            end
        case 'p' % Go to specific UID
            newPair = str2double(cell2mat(inputdlg('Go to pair:')));
            if ~ismember(newPair,unique(guiData.pairIDs))
                error(['Pair ' num2str(newPair) ' not present'])
            end
            guiData.curr.pair = newPair;
            guiData.currPairPosition  = guiData.pairIDs == guiData.curr.pair;          
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
tmpMatchIdx = arrayfun(@(X)find(MatchTable.UID1 == X & MatchTable.UID2 == X,1,'first'),guiData.pairIDs);
eval(['MatchTable.' guiData.name '(tmpMatchIdx) = match;']);
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