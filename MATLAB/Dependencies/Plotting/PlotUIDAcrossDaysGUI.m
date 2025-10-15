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
UniqueIDConversion = struct();
WaveformInfo = struct();
UMparam = struct();

%% Load data from UnitMatch output
loadData();

%% Create GUI components
f = figure('Name','UID across days', 'NumberTitle','off', 'Color','w', ...
    'KeyPressFcn', @keyPress);

uicontrol('Style','text', 'Units','normalized', 'Position',[0.01 0.94 0.06 0.04], 'String','SaveDir');
saveDirEdit = uicontrol('Style','edit', 'Units','normalized', 'Position',[0.08 0.94 0.58 0.05], ...
    'String',SaveDir);
uicontrol('Style','pushbutton', 'String','Browse', 'Units','normalized', 'Position',[0.67 0.94 0.1 0.05], ...
    'Callback',@browseSaveDir);

uicontrol('Style','text', 'Units','normalized', 'Position',[0.01 0.88 0.05 0.04], 'String','UID');
uidEdit = uicontrol('Style','edit', 'Units','normalized', 'Position',[0.07 0.88 0.1 0.05]);
plotBtn = uicontrol('Style','pushbutton','String','Plot','Units','normalized','Position',[0.18 0.88 0.1 0.05], ...
    'Callback',@plotCallback);
uicontrol('Style','pushbutton','String','<','Units','normalized','Position',[0.29 0.88 0.05 0.05], ...
    'Callback',@prevCallback);
uicontrol('Style','pushbutton','String','>','Units','normalized','Position',[0.35 0.88 0.05 0.05], ...
    'Callback',@nextCallback);

axWave = axes('Parent', f, 'Position', [0.1 0.55 0.85 0.35]);
axISI  = axes('Parent', f, 'Position', [0.1 0.1 0.85 0.35]);
xlabel(axWave, 'Time (ms)'); ylabel(axWave, 'Amplitude');
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
        % Time vector in milliseconds for waveforms
        try
            nTime = size(WaveformInfo.ProjectedWaveform, 1);
        catch
            nTime = numel(WaveformInfo.ProjectedWaveform(:,1,1));
        end
        fs = getWaveformFs(); % samples per second
        tms = (0:nTime-1) / fs * 1000; % milliseconds
        % Track data ranges to autoscale axes
        wfYMin = inf; wfYMax = -inf; hasWave = false;
        ISIbins = 5*10.^(-4:0.1:0); % seconds
        isiXms = ISIbins * 1000; hasISI = false; isiYMax = 0;
        for k = 1:numel(idx)
            rec = RecAllGood(idx(k));
            dateStr = getRecDate(rec);
            % Waveform: average across cross-validations
            wf = squeeze(nanmean(WaveformInfo.ProjectedWaveform(:, idx(k), :), 3));
            plot(axWave, tms, wf, 'Color', colors(k,:), 'DisplayName', dateStr);
            wfValid = wf(~isnan(wf));
            if ~isempty(wfValid)
                wfYMin = min(wfYMin, min(wfValid));
                wfYMax = max(wfYMax, max(wfValid));
                hasWave = true;
            end

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
            histISI = histcounts(isi, [ISIbins Inf], 'Normalization','probability');
            stairs(axISI, ISIbins*1000, histISI, 'Color', colors(k,:), ...
                'DisplayName', dateStr);
            if ~isempty(histISI)
                isiYMax = max(isiYMax, max(histISI));
                hasISI = true;
            end
        end
        hold(axWave, 'off');
        hold(axISI, 'off');
        % Autoscale axes based on plotted data
        if hasWave
            xlim(axWave, [tms(1) tms(end)]);
            yrng = wfYMax - wfYMin;
            if ~isfinite(yrng) || yrng == 0
                pad = max(1e-12, 0.05*max(1, abs(wfYMax)));
            else
                pad = 0.05 * yrng;
            end
            ylim(axWave, [wfYMin - pad, wfYMax + pad]);
        end
        if hasISI
            xlim(axISI, [isiXms(1) isiXms(end)]);
            ymax = max(1e-6, isiYMax*1.1);
            ylim(axISI, [0 ymax]);
        end
        legend(axWave, 'show');
        legend(axISI, 'show');

        axes(axWave); makepretty; offsetAxes(axWave);
        axes(axISI); makepretty; offsetAxes(axISI);
    end

    function fs = getWaveformFs()
        % Determine sampling rate (Hz) for waveforms; default to 30 kHz
        fs = 30000; % default
        if exist('UMparam','var') && isstruct(UMparam)
            candidates = {'SampleRate','SamplingRate','APSampleRate','Fs','fs','spikeSampleRate','WaveformFs'};
            for c = 1:numel(candidates)
                f = candidates{c};
                if isfield(UMparam, f)
                    val = UMparam.(f);
                    if isnumeric(val) && isscalar(val) && isfinite(val) && val > 0
                        fs = val; return
                    end
                end
            end
        end
    end

    function dateStr = getRecDate(rec)
        dateStr = sprintf('Rec %d', rec);
        if isfield(UMparam,'KSDir') && numel(UMparam.KSDir) >= rec
            parts = strsplit(UMparam.KSDir{rec}, filesep);
            mask = cellfun(@(p) ~isempty(regexp(p,'\d{4}-\d{2}-\d{2}','once')), parts);
            if any(mask)
                lastIdx = find(mask,1,'last');
                % Extract only the yyyy-mm-dd substring from the matching part
                m = regexp(parts{lastIdx}, '(\d{4}-\d{2}-\d{2})', 'match', 'once');
                if ~isempty(m)
                    dateStr = m;
                end
            end
        end
    end
end
