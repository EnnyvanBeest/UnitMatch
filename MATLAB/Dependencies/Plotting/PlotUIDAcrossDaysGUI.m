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

if nargin < 1 || ~isfolder(SaveDir)
    error('SaveDir must be a valid directory');
end

%% Load data from UnitMatch output
umFile = fullfile(SaveDir, 'UnitMatch.mat');
if ~isfile(umFile)
    error('UnitMatch.mat not found in %s', SaveDir);
end
load(umFile, 'UniqueIDConversion', 'WaveformInfo', 'UMparam');

%% Create GUI components
f = figure('Name','UID across days', 'NumberTitle','off', 'Color','w');

uicontrol('Style','text', 'Position',[10 370 40 20], 'String','UID');
uidEdit = uicontrol('Style','edit', 'Position',[50 370 80 25]);
plotBtn = uicontrol('Style','pushbutton','String','Plot','Position',[140 370 60 25], ...
    'Callback',@plotCallback);

axWave = axes('Parent', f, 'Position', [0.1 0.55 0.85 0.35]);
axISI  = axes('Parent', f, 'Position', [0.1 0.1 0.85 0.35]);
xlabel(axWave, 'Time (samples)'); ylabel(axWave, 'Amplitude');
xlabel(axISI, 'ISI (ms)');      ylabel(axISI, 'Probability');

%% Callback function
    function plotCallback(~,~)
        uid = str2double(uidEdit.String);
        cla(axWave); cla(axISI);
        if isnan(uid)
            title(axWave, 'Enter a numeric UID');
            return
        end

        idx = find(UniqueIDConversion.UniqueID(logical(UniqueIDConversion.GoodID)) == uid);
        if isempty(idx)
            title(axWave, sprintf('UID %d not found', uid));
            return
        end
        RecAllGood = UniqueIDConversion.recsesAll(logical(UniqueIDConversion.GoodID));
        OriClusGood = UniqueIDConversion.OriginalClusID(logical(UniqueIDConversion.GoodID));

        colors = flipud(copper(numel(idx)));
        hold(axWave, 'on');
        hold(axISI, 'on');
        for k = 1:numel(idx)
            rec = RecAllGood(idx(k));
            % Waveform: average across cross-validations
            wf = squeeze(nanmean(WaveformInfo.ProjectedWaveform(:, idx(k), :), 3));
            plot(axWave, wf, 'Color', colors(k,:), 'DisplayName', sprintf('Rec %d', rec));

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
            ISIbins = [0 5*10.^(-4:0.1:0)];
            histISI = histcounts(isi, ISIbins, 'Normalization','probability');
            stairs(axISI, ISIbins(1:end-1)*1000, histISI, 'Color', colors(k,:), ...
                'DisplayName', sprintf('Rec %d', rec));
        end
        hold(axWave, 'off');
        hold(axISI, 'off');
        legend(axWave, 'show');
        legend(axISI, 'show');
    end

end
