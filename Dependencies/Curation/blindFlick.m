function blindFlick(pathFiles,user,recompute)
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

    if ~exist('user','var')
        warning('No user specified. Will use default.')
        user = 'default';
    end

    if ~exist('recompute','var')
        recompute = 0;
    end
    
    d = dir(fullfile(pathFiles,'*.fig'));
    
    blindFlickGUI = figure('color','w');
    guiData = struct;

    guiData.d = d; 
    guiData.curr.pairList = 1:numel(d);
    guiData.curr.pair = 1;
    guiData.curr.updateFig = 1;
    guiData.matchPath = fullfile(pathFiles,sprintf('manualCuration_%s.mat',user));
    if exist(guiData.matchPath,'file') && ~recompute
        tmp = load(guiData.matchPath);
        guiData.match = tmp.match;
    else
        guiData.match = zeros(1,numel(d));
    end
    guiData.curr.match = guiData.match(1);

    guidata(blindFlickGUI, guiData);
    updatePlot(blindFlickGUI);
end

function updatePlot(blindFlickGUI)
    % Get guidata
    guiData = guidata(blindFlickGUI);

    if guiData.curr.updateFig == 1
        % Clear previous fig
        clf(blindFlickGUI)

        % Load and copy new one (slow -- to change asap)
        fig = openfig(fullfile(guiData.d(guiData.curr.pair).folder,guiData.d(guiData.curr.pair).name),'invisible');
        axNew = findobj(fig,'type','axes');
        copyobj(axNew,blindFlickGUI)
        close(fig)
    end

    % Update title with label -- hacky
    if ~isfield(guiData,'titleMain')
        guiData.titleMain = annotation('textbox', [0, 1, 1, 0], 'string', 'My Text', 'EdgeColor', 'none',...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    set(guiData.titleMain, 'String', sprintf('pair %d: match? %d',guiData.curr.pair,guiData.curr.match))

    % Set functions for key presses
    set(blindFlickGUI,'WindowKeyPressFcn',@keyPress);
end

function keyPress(blindFlickGUI,eventdata)
    % Get guidata
    guiData = guidata(blindFlickGUI);
    
    guiData.curr.updateFig = 0;
    switch eventdata.Key
        case 'rightarrow' % Next pair
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            newPair = guiData.curr.pairList(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'leftarrow' % Previous pair
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            newPair = guiData.curr.pairList(circshift(currPairPosition,-1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'uparrow' % It's a match
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            guiData.match(currPairPosition) = 1;
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            newPair = guiData.curr.pairList(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'downarrow' % It's not a match
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            guiData.match(currPairPosition) = -1;
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            newPair = guiData.curr.pairList(circshift(currPairPosition,1));
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
        case 'o'
            currPairPosition = guiData.curr.pairList == guiData.curr.pair;
            guiData.match(currPairPosition) = 0;
        case 'm'
            guiData.curr.pair = find(guiData.match == 0,1);
            guiData.curr.updateFig = 1;
        case 'p'
            newPair = str2double(cell2mat(inputdlg('Go to pair:')));
            if ~ismember(newPair,unique(guiData.curr.pairList))
                error(['Pair ' num2str(newPair) ' not present'])
            end
            guiData.curr.pair = newPair;
            guiData.curr.updateFig = 1;
    end
    currPairPosition = guiData.curr.pairList == guiData.curr.pair;
    guiData.curr.match = guiData.match(currPairPosition);
    match = guiData.match;
    save(guiData.matchPath,'match')

    guidata(blindFlickGUI,guiData);
    updatePlot(blindFlickGUI);
end
