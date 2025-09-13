function PlotUIDAcrossDaysGUI(SaveDir)
% PlotUIDAcrossDaysGUI Launch a GUI to visualize waveforms and ISI histograms
% across days for a given unique ID (UID).
%
%   PlotUIDAcrossDaysGUI(SaveDir) loads UnitMatch output from SaveDir and
%   opens a simple graphical interface where a UID can be entered. The
%   waveforms and inter-spike interval (ISI) histograms for that UID are
%   plotted for all recording sessions in which it was detected.
%
% This utility relies on the UnitMatch output structure and expects the
% following variables to be present in SaveDir/UnitMatch.mat:
%   - UniqueIDConversion: struct containing fields UniqueID, OriginalClusID,
%     recsesAll and Path4UnitNPY.
%   - WaveformInfo: struct containing field ProjectedWaveform (time x unit x cv).
%   - UMparam: struct containing field KSDir with paths to recordings.
%
% Example:
%   PlotUIDAcrossDaysGUI('C:\Data\UnitMatchOutput');
%
% coded with Codex

%% Get save directory
if nargin < 1 || ~isfolder(SaveDir)
    [file, path] = uigetfile('UnitMatch.mat');
    if isequal(file,0)
        error('No UnitMatch.mat selected');
    end
    SaveDir = path;
end

UIDList = [];
currentIdx = 1;

%% Load data from UnitMatch output
loadData();

%% Create GUI components
f = figure('Name','UID across days', 'NumberTitle','off', 'Color','w', ...
    'KeyPressFcn', @keyPress);

uicontrol('Style','text', 'Position',[10 400 60 20], 'String','SaveDir');
saveDirEdit = uicontrol('Style','edit', 'Position',[70 400 320 25], ...
    'String',SaveDir);
uicontrol('Style','pushbutton', 'String','Browse', 'Position',[395 400 60 25], ...
    'Callback',@browseSaveDir);

uicontrol('Style','text', 'Position',[10 370 40 20], 'String','UID');
uidEdit = uicontrol('Style','edit', 'Position',[50 370 80 25]);
plotBtn = uicontrol('Style','pushbutton','String','Plot','Position',[140 370 60 25], ...
    'Callback',@plotCallback);
uicontrol('Style','pushbutton','String','<','Position',[205 370 40 25], ...
    'Callback',@prevCallback);
uicontrol('Style','pushbutton','String','>','Position',[250 370 40 25], ...
    'Callback',@nextCallback);

axWave = axes('Parent', f, 'Position', [0.1 0.55 0.85 0.35]);
axISI  = axes('Parent', f, 'Position', [0.1 0.1 0.85 0.35]);
xlabel(axWave, 'Time (samples)'); ylabel(axWave, 'Amplitude');
xlabel(axISI, 'ISI (ms)');      ylabel(axISI, 'Probability');
set(axISI,'XScale','log');

%% Nested functions
    function loadData()
        umFile = fullfile(SaveDir, 'UnitMatch.mat');
        if ~isfile(umFile)
            error('UnitMatch.mat not found in %s', SaveDir);
        end
        load(umFile, 'UniqueIDConversion', 'WaveformInfo', 'UMparam');
        UIDList = unique(UniqueIDConversion.UniqueID(logical(UniqueIDConversion.GoodID)));
        currentIdx = 1;
    end

    function browseSaveDir(~,~)
        [file, path] = uigetfile('UnitMatch.mat');
        if isequal(file,0); return; end
        SaveDir = path;
        saveDirEdit.String = SaveDir;
        loadData();
    end

    function plotCallback(~,~)
        uid = str2double(uidEdit.String);
        idx = find(UIDList == uid,1);
        if ~isempty(idx)
            currentIdx = idx;
        else
            title(axWave, sprintf('UID %d not found', uid));
            return
        end
        plotCurrent();
    end

    function prevCallback(~,~)
        currentIdx = mod(currentIdx-2, numel(UIDList)) + 1;
        uidEdit.String = num2str(UIDList(currentIdx));
        plotCurrent();
    end

    function nextCallback(~,~)
        currentIdx = mod(currentIdx, numel(UIDList)) + 1;
        uidEdit.String = num2str(UIDList(currentIdx));
        plotCurrent();
    end

    function keyPress(~,event)
        switch event.Key
            case 'rightarrow'
                nextCallback([],[]);
            case 'leftarrow'
                prevCallback([],[]);
        end
    end

    function plotCurrent()
        uid = UIDList(currentIdx);
        cla(axWave); cla(axISI);

        idx = find(UniqueIDConversion.UniqueID(logical(UniqueIDConversion.GoodID)) == uid);
        RecAllGood = UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID));
        OriClusGood = UniqueIDConversion.OriginalClusID(logical(UniqueIDConversion.GoodID));

        colors = flipud(copper(numel(idx)));
        hold(axWave, 'on');
        hold(axISI, 'on');
        for k = 1:numel(idx)
            rec = RecAllGood(idx(k));
            dateStr = getRecDate(rec);
            % Waveform: average across cross-validations
            wf = squeeze(nanmean(WaveformInfo.ProjectedWaveform(:, idx(k), :), 3));
            plot(axWave, wf, 'Color', colors(k,:), 'DisplayName', dateStr);

            % Load spikes for ISI calculation
            try
                tmp = load(fullfile(UMparam.KSDir{rec}, 'PreparedData.mat'), 'sp');
            catch
                warning('PreparedData.mat missing for session %d', rec);
                continue
            end
            sp = tmp.sp;
            if isfield(sp,'clu')
                clu = sp.clu;
            else
                clu = sp.spikeTemplates; % fallback
            end
            st = double(sp.st(clu == OriClusGood(idx(k))));
            if numel(st) < 2
                continue
            end
            isi = diff(sort(st));
            ISIbins = 5*10.^(-4:0.1:0);
            histISI = histcounts(isi, [ISIbins Inf], 'Normalization','probability');
            stairs(axISI, ISIbins*1000, histISI, 'Color', colors(k,:), ...
                'DisplayName', dateStr);
        end
        hold(axWave, 'off');
        hold(axISI, 'off');
        legend(axWave, 'show');
        legend(axISI, 'show');

        axes(axWave); makepretty; offsetAxes(axWave);
        axes(axISI); makepretty; offsetAxes(axISI);
    end

    function dateStr = getRecDate(rec)
        dateStr = sprintf('Rec %d', rec);
        if isfield(UMparam,'KSDir') && numel(UMparam.KSDir) >= rec
            parts = strsplit(UMparam.KSDir{rec}, filesep);
            mask = cellfun(@(p) ~isempty(regexp(p,'\d{4}-\d{2}-\d{2}','once')), parts);
            if any(mask)
                dateStr = parts{find(mask,1,'last')};
            end
        end
    end
end
