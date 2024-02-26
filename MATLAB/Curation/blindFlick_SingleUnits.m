function blindFlick_SingleUnits(pathFiles,user,recompute)
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
    %   f: toggle functional scores

    if ~exist('user','var')
        warning('No user specified. Will use default.')
        user = 'default';
    end

    if ~exist('recompute','var')
        recompute = 0;
    end
    
    % Get list of pair IDs
    d = dir(fullfile(pathFiles,'*.fig'));
    IDs = cell2mat(cellfun(@(y) str2num(y{1}(2:end-1)), cellfun(@(x) regexp(x,'_\d+.','match'), {d.name}, 'uni',0),'uni',0));
    [~,sortIdx] = sort(IDs,'ascend');
    d = d(sortIdx);
    IDs = IDs(sortIdx);
    
    % Fill in GUI data
    guiData = struct;
    guiData.d = d;
    guiData.matchPath = fullfile(pathFiles,sprintf('manualCuration_%s.mat',user));
    if exist(guiData.matchPath,'file') && ~recompute
        tmp = load(guiData.matchPath);
        guiData.match = tmp.match;
        if isfield(tmp,'IDs') && ((numel(tmp.IDs) ~= numel(IDs)) || ~all(tmp.IDs == IDs))
            error('Saved IDs don''t match current IDs in the folder. Recheck?')
        else
            guiData.pairIDs = tmp.IDs;
        end
    else
        guiData.match = zeros(1,numel(d));
        guiData.pairIDs = IDs;
        save(guiData.matchPath,'IDs')
    end
    guiData.curr.pair = guiData.pairIDs(1); % Start at first
    guiData.curr.updateFig = 1;
    guiData.curr.match = guiData.match(guiData.pairIDs == guiData.curr.pair);
    guiData.showFinishBox = 1;

    % Create figure
    blindFlickGUI = figure('color','w');
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
        currPairPosition = guiData.pairIDs == guiData.curr.pair;
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
    set(guiData.titleMain, 'String', sprintf('pair ID %d: match? %d (%d/%d left)', ...
        guiData.curr.pair,guiData.curr.match,sum(guiData.match == 0),numel(guiData.pairIDs)));

    % Set functions for key presses
    set(blindFlickGUI,'WindowKeyPressFcn',@keyPress);
end

function keyPress(blindFlickGUI,eventdata)
    % Get guidata
    guiData = guidata(blindFlickGUI);
    
    guiData.curr.updateFig = 0;
    switch eventdata.Key
        case 'rightarrow' % Next pair
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'leftarrow' % Previous pair
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,-1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'uparrow' % It's a match
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = 1;
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'downarrow' % It's not a match
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = -1;
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            newPair = guiData.pairIDs(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'o' % Back to uncurated
            currPairPosition = guiData.pairIDs == guiData.curr.pair;
            guiData.match(currPairPosition) = 0;
        case 'm' % Go to first uncurated
            firstUncuratedPair = find(guiData.match == 0,1);
            if ~isempty(firstUncuratedPair)
                guiData.curr.pair = firstUncuratedPair;
                guiData.curr.updateFig = 1;
            else
                guiData.showFinishBox = 1;
            end
        case 'p' % Go to specific pair
            newPair = str2double(cell2mat(inputdlg('Go to pair:')));
            if ~ismember(newPair,unique(guiData.pairIDs))
                error(['Pair ' num2str(newPair) ' not present'])
            end
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'f' % Toggle functional scores
           keyboard %To be made
    end
    currPairPosition = guiData.pairIDs == guiData.curr.pair;
    guiData.curr.match = guiData.match(currPairPosition);
    match = guiData.match;
    save(guiData.matchPath,'match','-append')

    if ~any(guiData.match == 0) && guiData.showFinishBox
        msgbox('All pairs have been curated.')
        guiData.showFinishBox = 0;
    end

    guidata(blindFlickGUI,guiData);
    updatePlot(blindFlickGUI);
end
