function fig = AlignmentCurationGUI(matFile, allenCCFPath, histoRoot)
% ALIGNMENTCURATIONGUI  Standalone interactive GUI to review and correct
% automated histology-to-ephys alignment stored in *_HistoEphysAlignment_Auto.mat.
%
% Usage:
%   AlignmentCurationGUI(matFile)
%   AlignmentCurationGUI(matFile, allenCCFPath)
%   AlignmentCurationGUI(matFile, allenCCFPath, histoRoot)
%
% Inputs:
%   matFile      - full path to *_HistoEphysAlignment_Auto.mat
%   allenCCFPath - (optional) folder containing:
%                    annotation_volume_10um_by_index.npy
%                    structure_tree_safe_2017.csv
%                  If omitted, searches the MATLAB path automatically.
%   histoRoot    - (optional) root folder for per-mouse brainreg output;
%                  default 'H:\Histology'. Expected structure:
%                  <histoRoot>/<Mouse>/brainreg_output/downsampled_standard.tiff
%
% What the GUI shows:
%   - Probe depth bar (left): area-coloured patches, draggable white boundary lines
%   - Coronal + Sagittal brain slices: probe track overlaid, coloured by area
%     Uses per-mouse brainreg histology when available, falls back to Allen CCF atlas
%   - Functional features panel (right, only when AlignmentFeatures was saved)
%
% Interactions:
%   Drag boundary lines  - reassign areas at probe depth transitions
%   Shift slider         - rigid up/down offset to all depths (µm)
%   Scale slider         - stretch/compress probe depth range around midpoint
%   Reset button         - restore automated alignment
%   Accept & Save        - write updated Depth2Area + ManuallyCurated=true
%   Flag for Review      - write ManuallyCurated=false + FlaggedForReview=true
%
% Saved variables (appended to matFile on Accept or Flag):
%   Depth2Area       - updated table
%   ManuallyCurated  - true (Accept) or false (Flag)
%   FlaggedForReview - false (Accept) or true (Flag)
%   CurationNotes    - string (from Flag dialog), else ''
%   CurationDate     - timestamp string
%
% ManuallyCurated missing or false in an existing file means the alignment
% has not yet been manually reviewed.

if nargin < 2, allenCCFPath = []; end
if nargin < 3 || isempty(histoRoot), histoRoot = 'H:\Histology'; end

%% ── Load alignment mat file ──────────────────────────────────────────────────
raw = load(matFile);
if ~isfield(raw,'Depth2Area') || isempty(raw.Depth2Area)
    error('AlignmentCurationGUI: %s does not contain Depth2Area.', matFile);
end
D2A_orig  = sortrows(raw.Depth2Area, {'Shank','Depth'});
D2A       = D2A_orig;

wasCurated       = isfield(raw,'ManuallyCurated') && raw.ManuallyCurated;
wasFlagged       = isfield(raw,'FlaggedForReview') && raw.FlaggedForReview;
AlignmentFeatures = {};
if isfield(raw,'AlignmentFeatures') && ~isempty(raw.AlignmentFeatures)
    AF = raw.AlignmentFeatures;
    % Normalise storage format to a cell vector {nShanks × 1} of numeric matrices.
    % The field may be saved as: a plain numeric matrix (single shank),
    % or a cell array {nShanks}.
    if isnumeric(AF)
        AlignmentFeatures = {AF};
    elseif iscell(AF)
        AlignmentFeatures = AF(:);   % ensure column cell
    end
end
clear raw

% Infer mouse name from path structure:  .../SaveDir/<Mouse>/<Date>/<Probe>/
[matDir, ~, ~] = fileparts(matFile);
parts = strsplit(matDir, filesep);
mouseName = '';
if numel(parts) >= 3, mouseName = parts{end-2}; end

%% ── Load Allen CCF annotation volume ─────────────────────────────────────────
av = []; avSz = [0 0 0]; atlasCmap = []; st = [];
if isempty(allenCCFPath)
    tmp = which('annotation_volume_10um_by_index.npy');
    if ~isempty(tmp), allenCCFPath = fileparts(tmp); end
end
if ~isempty(allenCCFPath)
    avF = fullfile(allenCCFPath, 'annotation_volume_10um_by_index.npy');
    stF = fullfile(allenCCFPath, 'structure_tree_safe_2017.csv');
    if exist(avF,'file')
        disp('AlignmentCurationGUI: loading Allen CCF...')
        av   = readNPY(avF);
        avSz = size(av);
        disp('Done.')
    end
    if exist(stF,'file')
        st = readtable(stF);
        nSt = height(st);
        atlasCmap = 0.12 * ones(nSt+1, 3);
        atlasCmap(1,:) = [0 0 0];
        for ri = 1:nSt
            h6 = local_getHex(st, ri);
            if numel(h6) >= 6
                atlasCmap(ri+1,:) = [hex2dec(h6(1:2)) hex2dec(h6(3:4)) hex2dec(h6(5:6))]/255;
            end
        end
    end
end
if isempty(atlasCmap) && ~isempty(av)
    atlasCmap = gray(double(max(av(:)))+1);
end

%% ── Load per-mouse histology volume ──────────────────────────────────────────
hVol = local_loadHistoVol(mouseName, histoRoot);
useHistology = ~isempty(hVol) && ~isempty(hVol.vol);
useAtlas     = ~isempty(av);

%% ── Mutable state ────────────────────────────────────────────────────────────
boundaries    = computeBoundaries(D2A);  % per-shank cell array of boundary structs
dragState     = struct('active',false,'shankIdx',0,'boundaryIdx',0);
lastShift     = 0;    % tracks cumulative slider shift so changes are incremental
lastScale     = 1.0;  % tracks cumulative slider scale
allShanks     = unique(D2A.Shank);
selectedShank = allShanks(1);   % shank whose track is shown in the brain slices
atlasAlpha    = 0.4;            % transparency of atlas colour overlay (0 = off)
brainZoom     = 0.25;            % 1 = full slice view; <1 = zoomed in around probe

%% ── Build figure ─────────────────────────────────────────────────────────────
% hasFeatures is true only if at least one shank has a non-empty numeric matrix
hasFeatures = ~isempty(AlignmentFeatures) && ...
    any(cellfun(@(c) isnumeric(c) && ~isempty(c), AlignmentFeatures));
nRightCols  = 2 + hasFeatures;  % coronal + sagittal [+ features]
nCols       = 1 + nRightCols;

fig = uifigure('Name', ['Alignment Curation  —  ' matFile], ...
    'Position', [30 50 min(1800, 380*nCols) 860]);
fig.WindowButtonDownFcn   = @onMouseDown;
fig.WindowButtonMotionFcn = @onMouseMove;
fig.WindowButtonUpFcn     = @onMouseUp;
fig.Pointer               = 'arrow';

colWidths = ['0.65x', repmat({'1x'}, 1, nRightCols)];
outerG = uigridlayout(fig, [2, nCols], ...
    'RowHeight',   {'7x','1x'}, ...
    'ColumnWidth', colWidths, ...
    'Padding', [4 4 4 4], 'ColumnSpacing', 5, 'RowSpacing', 4);

% ── Probe depth bar ──────────────────────────────────────────────────────────
axProbe = uiaxes(outerG);
axProbe.Layout.Row = 1; axProbe.Layout.Column = 1;
hold(axProbe,'on');
set(axProbe, 'YDir','normal', 'XTick',[], 'Box','off', 'FontSize',8);
title(axProbe, 'Probe  —  drag boundary lines here');
ylabel(axProbe, 'Depth on probe (µm)');

% ── Brain slice axes ─────────────────────────────────────────────────────────
axCoronal = uiaxes(outerG);
axCoronal.Layout.Row = 1; axCoronal.Layout.Column = 2;
hold(axCoronal,'on'); axis(axCoronal,'image');
set(axCoronal, 'YDir','reverse', 'XTickLabel',{}, 'YTickLabel',{});
title(axCoronal,'Coronal'); xlabel(axCoronal,'ML'); ylabel(axCoronal,'DV');

axSagittal = uiaxes(outerG);
axSagittal.Layout.Row = 1; axSagittal.Layout.Column = 3;
hold(axSagittal,'on'); axis(axSagittal,'image');
set(axSagittal, 'YDir','reverse', 'XTickLabel',{}, 'YTickLabel',{});
title(axSagittal,'Sagittal'); xlabel(axSagittal,'AP'); ylabel(axSagittal,'DV');

% ── Features axis (optional) ─────────────────────────────────────────────────
axFeat = [];
if hasFeatures
    axFeat = uiaxes(outerG);
    axFeat.Layout.Row = 1; axFeat.Layout.Column = 4;
    hold(axFeat,'on');
    set(axFeat, 'YDir','normal', 'FontSize',8);
    title(axFeat,'Functional features');
    xlabel(axFeat,'Feature'); ylabel(axFeat,'Depth (µm)');
end

%% ── Control panel ────────────────────────────────────────────────────────────
% Layout: 2 rows × 11 cols
%   Cols 1-8  (both rows): Shift | Scale | Shank | Atlas-overlay sliders+labels
%   Cols 9-10 (both rows): View zoom slider
%   Cols 11-13 (row 1):    Reset | Flag | Accept buttons
%   Cols 11-13 (row 2):    status label
ctrlG = uigridlayout(outerG, [2, 13], ...
    'ColumnWidth', {'0.9x','3x','0.9x','3x','0.7x','1.2x','0.7x','2.5x','0.7x','2.5x','1.4x','1.6x','1.6x'}, ...
    'RowHeight',   {'1x','1x'}, ...
    'Padding', [2 2 2 2], 'ColumnSpacing', 5, 'RowSpacing', 2);
ctrlG.Layout.Row = 2; ctrlG.Layout.Column = [1 nCols];

% Shift
lblShift = uilabel(ctrlG, 'Text','Shift (µm):', 'HorizontalAlignment','right', 'FontWeight','bold');
lblShift.Layout.Row = [1 2]; lblShift.Layout.Column = 1;
slShift = uislider(ctrlG, 'Limits',[-500 500], 'Value',0, ...
    'MajorTicks',[-500 -250 0 250 500], 'ValueChangedFcn',@onShiftChanged);
slShift.Layout.Row = [1 2]; slShift.Layout.Column = 2;

% Scale
lblScale = uilabel(ctrlG, 'Text','Scale:', 'HorizontalAlignment','right', 'FontWeight','bold');
lblScale.Layout.Row = [1 2]; lblScale.Layout.Column = 3;
slScale = uislider(ctrlG, 'Limits',[0.5 2.0], 'Value',1.0, ...
    'MajorTicks',[0.5 1.0 1.5 2.0], 'ValueChangedFcn',@onScaleChanged);
slScale.Layout.Row = [1 2]; slScale.Layout.Column = 4;

% Shank selector (only interactive when >1 shank)
lblShank = uilabel(ctrlG, 'Text','Shank:', 'HorizontalAlignment','right', 'FontWeight','bold');
lblShank.Layout.Row = [1 2]; lblShank.Layout.Column = 5;
shankItems = arrayfun(@(s) sprintf('S%d',s), allShanks, 'UniformOutput',false);
ddShank = uidropdown(ctrlG, 'Items',shankItems, 'Value',shankItems{1}, ...
    'ValueChangedFcn',@onShankChanged, 'Enable', matlab.lang.OnOffSwitchState(numel(allShanks)>1));
ddShank.Layout.Row = [1 2]; ddShank.Layout.Column = 6;

% Atlas overlay alpha slider (only shown when both histology and atlas are available)
lblAtlas = uilabel(ctrlG, 'Text','Atlas:', 'HorizontalAlignment','right', 'FontWeight','bold');
lblAtlas.Layout.Row = [1 2]; lblAtlas.Layout.Column = 7;
canOverlay = useHistology && useAtlas;
slAtlas = uislider(ctrlG, 'Limits',[0 0.8], 'Value', 0.4 * canOverlay, ...
    'MajorTicks',[0 0.4 0.8], 'ValueChangedFcn',@onAtlasAlphaChanged, ...
    'Enable', matlab.lang.OnOffSwitchState(canOverlay));
slAtlas.Layout.Row = [1 2]; slAtlas.Layout.Column = 8;

% Brain view zoom
lblZoom = uilabel(ctrlG, 'Text','View:', 'HorizontalAlignment','right', 'FontWeight','bold');
lblZoom.Layout.Row = [1 2]; lblZoom.Layout.Column = 9;
slZoom = uislider(ctrlG, 'Limits',[0.05 1.0], 'Value',0.25, ...
    'MajorTicks',[0.05 0.25 0.5 0.75 1.0], 'ValueChangedFcn',@onZoomChanged);
slZoom.Layout.Row = [1 2]; slZoom.Layout.Column = 10;

% Buttons (row 1)
btnReset = uibutton(ctrlG, 'Text','Reset', ...
    'BackgroundColor',[0.45 0.45 0.45], 'FontColor','w', 'ButtonPushedFcn',@onReset);
btnReset.Layout.Row = 1; btnReset.Layout.Column = 11;
btnFlag = uibutton(ctrlG, 'Text','Flag for Review', ...
    'BackgroundColor',[0.80 0.55 0.0], 'FontColor','w', 'ButtonPushedFcn',@onFlag);
btnFlag.Layout.Row = 1; btnFlag.Layout.Column = 12;
btnAccept = uibutton(ctrlG, 'Text','Accept & Save', ...
    'BackgroundColor',[0.10 0.60 0.10], 'FontColor','w', 'ButtonPushedFcn',@onAccept);
btnAccept.Layout.Row = 1; btnAccept.Layout.Column = 13;

% Status label (row 2, under buttons)
if wasCurated
    initStatus = sprintf('Previously curated  (%s)', ...
        getfield_safe(load(matFile),'CurationDate','unknown date'));
elseif wasFlagged
    initStatus = 'Previously flagged for review.';
else
    initStatus = 'Automated alignment — not yet reviewed.';
end
lblStatus = uilabel(ctrlG, 'Text',initStatus, 'WordWrap','on', ...
    'HorizontalAlignment','left', 'FontSize',8, 'FontColor',[0.35 0.35 0.35]);
lblStatus.Layout.Row = 2; lblStatus.Layout.Column = [11 13];

%% ── Initial draw ─────────────────────────────────────────────────────────────
drawAll();

%% ════════════════════════════════════════════════════════════════════════════
%  NESTED FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════

    % ── Draw all panels ───────────────────────────────────────────────────────
    function drawAll()
        drawProbeBar();
        drawSlices();
        if hasFeatures, drawFeatures(); end
    end

    % ── Probe depth bar ───────────────────────────────────────────────────────
    function drawProbeBar()
        cla(axProbe);
        hold(axProbe,'on');

        shanks = unique(D2A.Shank);
        nS  = numel(shanks);
        gap = 0.05;
        sw  = (1 - gap*(nS+1)) / nS;

        for si = 1:nS
            shk  = shanks(si);
            mask = D2A.Shank == shk;
            idx  = find(mask);
            [~, ord] = sort(D2A.Depth(idx));
            idx = idx(ord);
            depths_s = D2A.Depth(idx);
            hex_s    = D2A.Color(idx);

            x0   = gap + (si-1)*(sw+gap);
            x1   = x0 + sw;
            xMid = (x0+x1)/2;

            % Area-coloured filled regions between boundaries.
            % Build segment list: each segment spans from one depth boundary
            % to the next, filled with the colour of the points in that range.
            bnd_si = boundaries{si};
            segEdges = [-Inf; bnd_si.depths(:); Inf];   % N+1 edges for N segments
            for seg = 1:numel(segEdges)-1
                lo = segEdges(seg);
                hi = segEdges(seg+1);
                % Points belonging to this segment
                inSeg = depths_s >= lo & depths_s < hi;
                if ~any(inSeg), continue; end
                % Use the colour of the first (shallowest) point in the segment
                c = hex2rgbSafe(hex_s{find(inSeg,1)});
                % Actual depth extent: clamp to the real probe range
                yLo = max(lo, depths_s(1));
                yHi = min(hi, depths_s(end));
                if isinf(yLo), yLo = depths_s(1); end
                if isinf(yHi), yHi = depths_s(end); end
                patch(axProbe, [x0 x1 x1 x0], [yLo yLo yHi yHi], c, ...
                    'EdgeColor','none', 'HitTest','off', 'PickableParts','none');
                % Area acronym label centred within this segment, to the right of the bar.
                % Pick a text colour that contrasts with the background.
                if mean(c) < 0.5
                    txtCol = [1 1 1];   % light text on dark background
                else
                    txtCol = [0 0 0];   % dark text on light background
                end
                segLabel = getParentAreaLocal(D2A.Area{idx(find(inSeg,1))});
                text(axProbe, xMid, (yLo + yHi)/2, segLabel, ...
                    'FontSize', 7, 'FontWeight', 'bold', 'Color', txtCol, ...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
                    'Clipping', 'on', 'HitTest', 'off');
            end

            % Shank label above column
            text(axProbe, xMid, max(depths_s)+30, sprintf('S%d',shk), ...
                'HorizontalAlignment','center', 'FontSize',7, 'FontWeight','bold', ...
                'HitTest','off');

            % Boundary lines with labels
            bnd = boundaries{si};
            for bi = 1:numel(bnd.depths)
                % Dark outline (drawn first, slightly thicker) for contrast
                plot(axProbe, [x0 x1], [bnd.depths(bi) bnd.depths(bi)], ...
                    '-', 'Color',[0 0 0], 'LineWidth', 5, 'HitTest','off', 'PickableParts','none');
                % White draggable line on top
                plot(axProbe, [x0 x1], [bnd.depths(bi) bnd.depths(bi)], ...
                    'w-', 'LineWidth', 2.5, 'HitTest','off', 'PickableParts','none');
                % Handle circle: dark ring + white fill
                plot(axProbe, xMid, bnd.depths(bi), 'ko', ...
                    'MarkerSize',11, 'MarkerFaceColor','k', 'HitTest','off', 'PickableParts','none');
                plot(axProbe, xMid, bnd.depths(bi), 'wo', ...
                    'MarkerSize',7, 'MarkerFaceColor','w', 'HitTest','off', 'PickableParts','none');
            end
        end

        set(axProbe, 'XLim',[0 1], 'YLim', ...
            [min(D2A.Depth)-30, max(D2A.Depth)+60]);

        % Highlight the shank whose track is shown in the brain slices
        % (drawn last so the yellow box is always on top)
        selIdx = find(allShanks == selectedShank, 1);
        if ~isempty(selIdx)
            x0sel = gap + (selIdx-1)*(sw+gap);
            x1sel = x0sel + sw;
            yl    = axProbe.YLim;
            plot(axProbe, [x0sel x1sel x1sel x0sel x0sel], ...
                [yl(1) yl(1) yl(2) yl(2) yl(1)], ...
                'Color',[1 1 0], 'LineWidth',2.5, 'HitTest','off', 'PickableParts','none');
        end

        title(axProbe, sprintf('Probe  shift=%.0fµm  scale=%.2fx  —  drag boundaries here', ...
            slShift.Value, slScale.Value));
    end

    % ── Brain slices ──────────────────────────────────────────────────────────
    function drawSlices()
        % Use the selected shank's coordinates to centre the slices.
        shankMask = D2A.Shank == selectedShank;
        coords    = D2A.Coordinates(shankMask, :);
        valid     = ~any(isnan(coords), 2);
        if ~any(valid)   % fallback: all shanks
            coords = D2A.Coordinates;
            valid  = ~any(isnan(coords), 2);
        end
        if ~any(valid), return; end
        medCoord = median(coords(valid,:), 1);   % [AP DV ML] in µm

        % Track overlay shows only the selected shank
        selSi = find(allShanks == selectedShank, 1);
        tc = buildTrackVox(D2A(shankMask,:), useHistology, hVol, avSz);

        % Boundary brain-space positions (interpolated along probe track)
        bnd_s = boundaries{selSi};
        bm    = buildBoundaryVox(D2A(shankMask,:), bnd_s.depths, useHistology, hVol, avSz);
        % Probe entry (top) and tip (bottom) as extent markers on the brain plots
        shankDepths_vm = D2A.Depth(shankMask);
        bm_extent = buildBoundaryVox(D2A(shankMask,:), ...
            [min(shankDepths_vm); max(shankDepths_vm)], useHistology, hVol, avSz);

        % ── Coronal slice (y=DV, x=ML) ───────────────────────────────────────
        cla(axCoronal);
        if useHistology
            vs   = hVol.voxSize_um;
            ap_h = max(1, min(size(hVol.vol,3), round(medCoord(1)/vs)));
            sl   = local_makeSlice(hVol, 'coronal', ap_h);
            imagesc(axCoronal, sl);
            if ~hVol.hasRGB, colormap(axCoronal, gray(256)); clim(axCoronal,[0 1]); end
            title(axCoronal, sprintf('Coronal  AP=%.0fµm  (S%d)', ap_h*vs, selectedShank));
            % Atlas colour overlay
            if atlasAlpha > 0 && useAtlas
                drawAtlasOverlay(axCoronal, 'coronal', ap_h, vs);
            end
        elseif useAtlas
            ap_a = max(1, min(avSz(1), round(medCoord(1)/10)));
            sl   = fliplr(squeeze(double(av(ap_a,:,:))));
            imagesc(axCoronal, sl);
            colormap(axCoronal, atlasCmap); clim(axCoronal,[0 size(atlasCmap,1)-1]);
            title(axCoronal, sprintf('Coronal  AP=%.0fµm  (S%d)', medCoord(1), selectedShank));
        end
        hold(axCoronal,'on');
        if ~isempty(tc)
            scatter(axCoronal, tc.ml, tc.dv, 20, tc.colors, 'filled', ...
                'MarkerFaceAlpha',0.85, 'HitTest','off');
        end
        if ~isempty(bm)
            xl_c = axCoronal.XLim;
            for bmi_ = 1:numel(bm.dv)
                plot(axCoronal, xl_c, [bm.dv(bmi_) bm.dv(bmi_)], 'w--', ...
                    'LineWidth',1.5, 'HitTest','off', 'PickableParts','none');
            end
            scatter(axCoronal, bm.ml, bm.dv, 200, 'k', 'd', 'filled', ...
                'HitTest','off', 'PickableParts','none');
            scatter(axCoronal, bm.ml, bm.dv, 100, [1 1 0], 'd', 'filled', ...
                'HitTest','off', 'PickableParts','none');
        end
        % Probe entry and tip: dotted lines spanning the full slice width
        if ~isempty(bm_extent)
            xl_ext = axCoronal.XLim;
            for ei = 1:numel(bm_extent.dv)
                plot(axCoronal, xl_ext, [bm_extent.dv(ei) bm_extent.dv(ei)], 'w:', ...
                    'LineWidth', 1.5, 'HitTest','off', 'PickableParts','none');
            end
        end
        axis(axCoronal,'image');
        set(axCoronal,'YDir','reverse','XTickLabel',{},'YTickLabel',{});
        xlabel(axCoronal,'ML (L→R)'); ylabel(axCoronal,'DV');
        % Apply zoom around probe centre (ML × DV)
        if brainZoom < 1 && ~isempty(tc)
            applySliceZoom(axCoronal, median(double(tc.ml)), median(double(tc.dv)));
        end

        % ── Sagittal slice (y=DV, x=AP) ──────────────────────────────────────
        cla(axSagittal);
        if useHistology
            ml_r = max(1, min(size(hVol.vol,2), round(medCoord(3)/vs)));
            sl   = local_makeSlice(hVol, 'sagittal', ml_r);
            imagesc(axSagittal, sl);
            if ~hVol.hasRGB, colormap(axSagittal, gray(256)); clim(axSagittal,[0 1]); end
            title(axSagittal, sprintf('Sagittal  ML=%.0fµm  (S%d)', ml_r*vs, selectedShank));
            % Atlas colour overlay
            if atlasAlpha > 0 && useAtlas
                drawAtlasOverlay(axSagittal, 'sagittal', ml_r, vs);
            end
        elseif useAtlas
            ml_a = max(1, min(avSz(3), round(medCoord(3)/10)));
            sl   = squeeze(double(av(:,:,ml_a)))';
            imagesc(axSagittal, sl);
            colormap(axSagittal, atlasCmap); clim(axSagittal,[0 size(atlasCmap,1)-1]);
            title(axSagittal, sprintf('Sagittal  ML=%.0fµm  (S%d)', medCoord(3), selectedShank));
        end
        hold(axSagittal,'on');
        if ~isempty(tc)
            scatter(axSagittal, tc.ap, tc.dv, 20, tc.colors, 'filled', ...
                'MarkerFaceAlpha',0.85, 'HitTest','off');
        end
        if ~isempty(bm)
            xl_s = axSagittal.XLim;
            for bmi_ = 1:numel(bm.dv)
                plot(axSagittal, xl_s, [bm.dv(bmi_) bm.dv(bmi_)], 'w--', ...
                    'LineWidth',1.5, 'HitTest','off', 'PickableParts','none');
            end
            scatter(axSagittal, bm.ap, bm.dv, 200, 'k', 'd', 'filled', ...
                'HitTest','off', 'PickableParts','none');
            scatter(axSagittal, bm.ap, bm.dv, 100, [1 1 0], 'd', 'filled', ...
                'HitTest','off', 'PickableParts','none');
        end
        % Probe entry and tip: dotted lines spanning the full slice width
        if ~isempty(bm_extent)
            xl_ext = axSagittal.XLim;
            for ei = 1:numel(bm_extent.dv)
                plot(axSagittal, xl_ext, [bm_extent.dv(ei) bm_extent.dv(ei)], 'w:', ...
                    'LineWidth', 1.5, 'HitTest','off', 'PickableParts','none');
            end
        end
        axis(axSagittal,'image');
        set(axSagittal,'YDir','reverse','XTickLabel',{},'YTickLabel',{});
        xlabel(axSagittal,'AP'); ylabel(axSagittal,'DV');
        % Apply zoom around probe centre (AP × DV)
        if brainZoom < 1 && ~isempty(tc)
            applySliceZoom(axSagittal, median(double(tc.ap)), median(double(tc.dv)));
        end
    end

    % ── Brain slice zoom helper ───────────────────────────────────────────────
    % Narrows the X-axis view to brainZoom fraction of its full width,
    % centred on cx (probe position in image x-coordinates).
    % Y-axis (DV) is always left at the full slice extent.
    function applySliceZoom(ax, cx, ~)
        xl = ax.XLim;
        hw = brainZoom * diff(xl) / 2;
        newXL = [cx - hw, cx + hw];
        % Clamp so we never pan outside the image
        if newXL(1) < xl(1), newXL = newXL - (newXL(1) - xl(1)); end
        if newXL(2) > xl(2), newXL = newXL - (newXL(2) - xl(2)); end
        xlim(ax, newXL);
    end

    % ── Atlas colour overlay on a histology background ────────────────────────
    % hVox is the histology voxel index (integer) along the slice-selecting axis.
    % vs   is hVol.voxSize_um (25 µm for brainreg).
    % The histology and atlas are both in CCF space, so physical extent matches;
    % we just need to map atlas pixel coordinates into histology pixel coordinates.
    function drawAtlasOverlay(ax, plane, hVox, vs)
        scale = 10 / vs;   % atlas pixels per histology pixel (e.g. 10/25 = 0.4)
        switch plane
            case 'coronal'
                ap_a    = max(1, min(avSz(1), round(hVox / scale)));
                rawSl   = fliplr(squeeze(double(av(ap_a,:,:))));  % [DV × ML]
                % map atlas [1..avSz(3)] pixels to histology x-range
                xd = [0.5,  avSz(3)*scale + 0.5];
                yd = [0.5,  avSz(2)*scale + 0.5];
            case 'sagittal'
                ml_a    = max(1, min(avSz(3), round(hVox / scale)));
                rawSl   = squeeze(double(av(:,:,ml_a)))';          % [DV × AP]
                xd = [0.5,  avSz(1)*scale + 0.5];
                yd = [0.5,  avSz(2)*scale + 0.5];
        end
        % Convert annotation indices to RGB (index 0 = outside brain → row 1 = black)
        atlasRGB = ind2rgb(rawSl + 1, atlasCmap);
        % Per-pixel alpha: outside brain fully transparent, inside at atlasAlpha
        alphaMap = atlasAlpha * double(rawSl > 0);
        hi = image(ax, 'XData',xd, 'YData',yd, 'CData',atlasRGB);
        hi.AlphaData    = alphaMap;
        hi.HitTest      = 'off';
        hi.PickableParts = 'none';
    end

    % ── Functional features ───────────────────────────────────────────────────
    function drawFeatures()
        cla(axFeat);

        % AlignmentFeatures is a cell {nShanks} of [nChan × nFeatures] matrices.
        % Pick the cell that corresponds to the currently selected shank.
        selSi = find(allShanks == selectedShank, 1);
        if isempty(selSi) || selSi > numel(AlignmentFeatures)
            selSi = 1;
        end
        F = AlignmentFeatures{selSi};

        % Validate: must be a non-empty numeric matrix
        if ~isnumeric(F) || isempty(F)
            text(axFeat, 0.5, 0.5, 'No features for this shank', ...
                'Units','normalized', 'HorizontalAlignment','center');
            title(axFeat, sprintf('Features  (S%d)', selectedShank));
            return;
        end
        F = double(F);

        % Depths for this shank, sorted ascending
        shankMask        = D2A.Shank == selectedShank;
        depths_s         = D2A.Depth(shankMask);
        [sortedD, sOrd]  = sort(depths_s, 'ascend');

        % Align feature rows to depth bins: trim or pad to match
        nDepths  = numel(sortedD);
        nFeatRows = size(F, 1);
        if nFeatRows >= nDepths
            % More (or equal) feature rows than depth bins → take first nDepths rows
            F = F(1:nDepths, :);
        else
            % Fewer feature rows than depth bins → pad with NaN at the bottom
            F = [F; NaN(nDepths - nFeatRows, size(F,2))];
        end
        % Reorder rows to match ascending depth order
        F = F(sOrd, :);

        % Per-column min-max normalisation (ignore NaNs)
        colMin = min(F, [], 1, 'omitnan');
        colMax = max(F, [], 1, 'omitnan');
        colRng = colMax - colMin;
        colRng(colRng < eps) = 1;   % avoid divide-by-zero for flat features
        F = (F - colMin) ./ colRng;

        imagesc(axFeat, 1:size(F,2), sortedD, F);
        set(axFeat, 'YDir','normal');
        colormap(axFeat, hot);

        % Match Y-axis exactly to the probe bar and fix X-axis
        xl_f = [0.5, size(F,2) + 0.5];
        set(axFeat, 'YLim', axProbe.YLim, 'XLim', xl_f, 'XTick', 1:size(F,2));

        % Overlay parent-area boundary lines with area-pair labels
        hold(axFeat, 'on');
        bnd_f = boundaries{selSi};
        for bi = 1:numel(bnd_f.depths)
            bd = bnd_f.depths(bi);
            plot(axFeat, xl_f, [bd bd], '-', 'Color','k', 'LineWidth',5, ...
                'HitTest','off', 'PickableParts','none');
            plot(axFeat, xl_f, [bd bd], 'w-', 'LineWidth',2, ...
                'HitTest','off', 'PickableParts','none');
            text(axFeat, xl_f(2), bd, ...
                sprintf('  %s | %s', bnd_f.parentBelow{bi}, bnd_f.parentAbove{bi}), ...
                'FontSize',6, 'Color','w', 'FontWeight','bold', ...
                'VerticalAlignment','middle', 'HitTest','off', 'Clipping','off');
        end

        xlabel(axFeat,'Feature index');
        ylabel(axFeat,'Depth (µm)');
        title(axFeat, sprintf('Features  (S%d)', selectedShank));
    end

    %% ── Mouse drag callbacks ─────────────────────────────────────────────────

    % Convert figure pixel coordinates to axProbe data coordinates.
    % In uifigures, axProbe.CurrentPoint is not reliably updated inside
    % WindowButton* callbacks, so we transform fig.CurrentPoint manually.
    function [cx, cy, inside] = figPtToAxProbe()
        axPos   = getpixelposition(axProbe, true); % [left bottom w h] in fig pixels
        figPx = fig.CurrentPoint;               % [x y] in figure pixels (from bottom-left)
        ti    = axProbe.TightInset;              % [left bottom right top] in pixels
        % Map into the inner data area so the transform is accurate at all depths
        innerLeft   = axPos(1) + ti(1);
        innerBottom = axPos(2) + ti(2);
        innerW      = axPos(3) - ti(1) - ti(3);
        innerH      = axPos(4) - ti(2) - ti(4);
        px = figPx(1) - innerLeft;
        py = figPx(2) - innerBottom;
        xl = axProbe.XLim;
        yl = axProbe.YLim;
        cx = xl(1) + px / innerW * diff(xl);
        cy = yl(1) + py / innerH * diff(yl);
        inside = (px >= 0) && (px <= innerW) && (py >= 0) && (py <= innerH);
    end

    function onMouseDown(~,~)
        [cx, cy, inside] = figPtToAxProbe();
        if ~inside, return; end

        yl = axProbe.YLim;
        nS  = numel(allShanks);
        gap = 0.05;
        sw  = (1 - gap*(nS+1)) / nS;
        snapTol = diff(yl) * 0.03;   % 3% of depth range

        % Only allow dragging boundaries of the currently selected shank
        si = find(allShanks == selectedShank, 1);
        if isempty(si), return; end
        x0 = gap + (si-1)*(sw+gap);
        x1 = x0 + sw;
        if cx < x0 || cx > x1, return; end
        bnd = boundaries{si};
        if isempty(bnd.depths), return; end
        [minDist, bi] = min(abs(bnd.depths - cy));
        if minDist < snapTol
            dragState.active      = true;
            dragState.shankIdx    = si;
            dragState.boundaryIdx = bi;
            fig.Pointer           = 'crosshair';
        end
    end

    function onMouseMove(~,~)
        [cx, cy, inside] = figPtToAxProbe();

        % Update cursor: show hand when hovering near a boundary of the selected shank
        if ~dragState.active
            nS      = numel(allShanks);
            gap     = 0.05;
            sw      = (1 - gap*(nS+1)) / nS;
            yl      = axProbe.YLim;
            snapTol = diff(yl) * 0.03;
            nearBnd = false;
            if inside
                si = find(allShanks == selectedShank, 1);
                if ~isempty(si)
                    x0 = gap + (si-1)*(sw+gap);
                    x1 = x0 + sw;
                    if cx >= x0 && cx <= x1
                        bnd = boundaries{si};
                        if ~isempty(bnd.depths) && min(abs(bnd.depths - cy)) < snapTol
                            nearBnd = true;
                        end
                    end
                end
            end
            fig.Pointer = conditional_ptr(nearBnd);
            return;
        end

        newY = cy;
        si   = dragState.shankIdx;
        bi   = dragState.boundaryIdx;
        bnd  = boundaries{si};

        % Clamp: keep ordering of adjacent boundaries (min 10µm gap)
        yLow  = axProbe.YLim(1);
        yHigh = axProbe.YLim(2);
        if bi > 1,              yLow  = bnd.depths(bi-1) + 10; end
        if bi < numel(bnd.depths), yHigh = bnd.depths(bi+1) - 10; end
        newY = max(yLow, min(yHigh, newY));

        % Preview: update boundary position without committing area labels
        boundaries{si}.depths(bi) = newY;
        drawProbeBar();
    end

    function onMouseUp(~,~)
        if ~dragState.active, return; end
        si       = dragState.shankIdx;
        bi       = dragState.boundaryIdx;
        newDepth = boundaries{si}.depths(bi);
        pa_below = boundaries{si}.parentBelow{bi};
        pa_above = boundaries{si}.parentAbove{bi};
        dragState.active = false;
        fig.Pointer = 'arrow';

        % Commit: reassign area labels for the selected shank
        D2A        = applyBoundaryMove(D2A, si, bi, newDepth);
        boundaries = computeBoundaries(D2A);

        % Sync the same parent-area boundary to all other shanks
        for osi = 1:numel(allShanks)
            if osi == si, continue; end
            bndO = boundaries{osi};
            for obi = 1:numel(bndO.depths)
                if strcmp(bndO.parentBelow{obi}, pa_below) && ...
                        strcmp(bndO.parentAbove{obi}, pa_above)
                    D2A        = applyBoundaryMove(D2A, osi, obi, newDepth);
                    boundaries = computeBoundaries(D2A);
                    break;
                end
            end
        end

        drawAll();
        lblStatus.Text = sprintf('Boundary moved.  %s', ...
            char(datetime('now','Format','HH:mm:ss')));
    end

    %% ── Slider / dropdown callbacks ──────────────────────────────────────────

    function onShankChanged(~,~)
        % Parse selected shank number from dropdown label (e.g. 'S2' → 2)
        selectedShank = sscanf(ddShank.Value, 'S%d');
        drawAll();   % redraw everything: probe bar (yellow box) + slices
    end

    function onAtlasAlphaChanged(~,~)
        atlasAlpha = slAtlas.Value;
        drawSlices();
    end

    function onZoomChanged(~,~)
        brainZoom = slZoom.Value;
        drawSlices();
    end

    function onShiftChanged(~,~)
        % Incremental shift preserves any boundary edits already applied
        newShift = slShift.Value;
        delta    = newShift - lastShift;
        lastShift = newShift;
        D2A.Depth  = D2A.Depth + delta;
        boundaries = computeBoundaries(D2A);
        drawAll();
    end

    function onScaleChanged(~,~)
        % Incremental scale around current midpoint, preserving boundary edits
        newScale = slScale.Value;
        if abs(lastScale) < 1e-9, return; end
        ratio    = newScale / lastScale;
        lastScale = newScale;
        mid       = mean([min(D2A.Depth), max(D2A.Depth)]);
        D2A.Depth  = mid + (D2A.Depth - mid) * ratio;
        boundaries = computeBoundaries(D2A);
        drawAll();
    end

    %% ── Button callbacks ─────────────────────────────────────────────────────

    function onReset(~,~)
        D2A        = D2A_orig;
        lastShift  = 0;
        lastScale  = 1.0;
        slShift.Value = 0;
        slScale.Value = 1.0;
        boundaries = computeBoundaries(D2A);
        drawAll();
        lblStatus.Text = 'Reset to automated alignment.';
    end

    function onAccept(~,~)
        doSave(true, false, '');
        lblStatus.Text = sprintf('Saved. ManuallyCurated=true  (%s)', ...
            char(datetime('now','Format','yyyy-MM-dd HH:mm')));
    end

    function onFlag(~,~)
        answer = inputdlg('Reason for flagging (optional):', 'Flag for Review', 1, {''});
        if isempty(answer), return; end
        doSave(false, true, answer{1});
        lblStatus.Text = sprintf('Flagged for review.  (%s)', ...
            char(datetime('now','Format','yyyy-MM-dd HH:mm')));
    end

    function doSave(curated, flagged, notes)
        Depth2Area       = D2A;
        ManuallyCurated  = curated;
        FlaggedForReview = flagged;
        CurationNotes    = notes;
        CurationDate     = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        save(matFile, 'Depth2Area','ManuallyCurated','FlaggedForReview', ...
            'CurationNotes','CurationDate', '-append');
        fprintf('AlignmentCurationGUI: saved to %s\n', matFile);
    end

    %% ── Boundary helpers ─────────────────────────────────────────────────────

    function bnd = computeBoundaries(D2A_in)
        shanks = unique(D2A_in.Shank);
        bnd    = cell(numel(shanks), 1);
        for si_ = 1:numel(shanks)
            mask = D2A_in.Shank == shanks(si_);
            idx  = find(mask);
            [~, ord] = sort(D2A_in.Depth(idx));
            idx  = idx(ord);
            areas       = D2A_in.Area(idx);
            parentAreas = cellfun(@getParentAreaLocal, areas, 'UniformOutput', false);
            % Only parent-area transitions are editable boundaries
            trans = find(~strcmp(parentAreas(1:end-1), parentAreas(2:end)));
            b.depths      = (D2A_in.Depth(idx(trans)) + D2A_in.Depth(idx(trans+1))) / 2;
            b.areaBelow   = areas(trans);
            b.areaAbove   = areas(trans+1);
            b.parentBelow = parentAreas(trans);
            b.parentAbove = parentAreas(trans+1);
            b.rowBelow    = idx(trans);
            b.rowAbove    = idx(trans+1);
            bnd{si_} = b;
        end
    end

    function D2A_out = applyBoundaryMove(D2A_in, si, bi, newBoundaryDepth)
        % Reassign rows that cross the parent-area boundary, and proportionally
        % remap child-area labels within each parent's new depth range so that
        % the sub-area sequence is preserved rather than collapsed to one label.
        D2A_out = D2A_in;
        shanks  = unique(D2A_in.Shank);
        shk     = shanks(si);
        bnd_    = boundaries{si};

        paBelow = bnd_.parentBelow{bi};
        paAbove = bnd_.parentAbove{bi};

        % Zone: depth range between the two neighbouring parent-area boundaries
        depthLow  = -Inf;
        depthHigh =  Inf;
        if bi > 1,                  depthLow  = bnd_.depths(bi-1); end
        if bi < numel(bnd_.depths), depthHigh = bnd_.depths(bi+1); end

        inZone  = D2A_in.Shank == shk & ...
                  D2A_in.Depth > depthLow & D2A_in.Depth < depthHigh;
        zoneIdx = find(inZone);
        if isempty(zoneIdx), return; end

        zonePAs    = cellfun(@getParentAreaLocal, D2A_in.Area(zoneIdx), 'UniformOutput', false);
        zoneDepths = D2A_in.Depth(zoneIdx);

        % ── Helper: remap child areas of one parent proportionally into newRange ──
        % templateIdx : row indices (into D2A_in) of all current members of this
        %               parent in the zone, sorted by depth — used as the child-area
        %               sequence template.
        % targetIdx   : row indices whose sub-area labels need to be assigned.
        % newRange    : [lo hi] depth interval the target rows now occupy.
        function remapChildren(templateIdx, targetIdx, newRange)
            if isempty(targetIdx), return; end
            % Sort both sets by depth
            [tplDepths, tOrd] = sort(D2A_in.Depth(templateIdx));
            tplAreas  = D2A_in.Area(templateIdx(tOrd));
            tplColors = D2A_in.Color(templateIdx(tOrd));
            tgtDepths = D2A_out.Depth(targetIdx);  % already committed depths

            % Map each target depth to a proportional position in the template range
            tplLo = tplDepths(1);   tplHi = tplDepths(end);
            if tplHi <= tplLo
                % Degenerate template: use only its first child label
                D2A_out.Area(targetIdx)  = tplAreas(1);
                D2A_out.Color(targetIdx) = tplColors(1);
                return;
            end
            for ti = 1:numel(targetIdx)
                % Normalised position within the target's new range → template depth
                frac     = (tgtDepths(ti) - newRange(1)) / max(diff(newRange), 1);
                frac     = max(0, min(1, frac));
                tplQuery = tplLo + frac * (tplHi - tplLo);
                % Find nearest template point
                [~, nearest] = min(abs(tplDepths - tplQuery));
                D2A_out.Area{targetIdx(ti)}  = tplAreas{nearest};
                D2A_out.Color{targetIdx(ti)} = tplColors{nearest};
            end
        end

        % ── paAbove rows that moved DOWN into paBelow territory ──────────────
        crossedDown = find(strcmp(zonePAs, paAbove) & (zoneDepths <= newBoundaryDepth));
        if ~isempty(crossedDown)
            % Template: all current paBelow rows in the zone
            tplMask = strcmp(zonePAs, paBelow);
            tplIdx  = zoneIdx(tplMask);
            tgtIdx  = zoneIdx(crossedDown);
            remapChildren(tplIdx, tgtIdx, [depthLow, newBoundaryDepth]);
        end

        % ── paBelow rows that moved UP into paAbove territory ────────────────
        crossedUp = find(strcmp(zonePAs, paBelow) & (zoneDepths > newBoundaryDepth));
        if ~isempty(crossedUp)
            % Template: all current paAbove rows in the zone
            tplMask = strcmp(zonePAs, paAbove);
            tplIdx  = zoneIdx(tplMask);
            tgtIdx  = zoneIdx(crossedUp);
            remapChildren(tplIdx, tgtIdx, [newBoundaryDepth, depthHigh]);
        end
    end

    function col = getAreaColor(D2A_in, area)
        idx = find(strcmp(D2A_in.Area, area), 1);
        if ~isempty(idx)
            col = D2A_in.Color{idx};
        else
            col = 'aaaaaa';
        end
    end

    %% ── Boundary brain-space coordinate helper ───────────────────────────────
    % Returns voxel coords (ap, dv, ml) for each boundary depth by linearly
    % interpolating along the probe track coordinates.
    function bm = buildBoundaryVox(D2A_shank, bndDepths, useHist, hV, avSz_)
        bm = [];
        if isempty(bndDepths), return; end
        coords = D2A_shank.Coordinates;
        depths = D2A_shank.Depth;
        valid  = ~any(isnan(coords), 2);
        coords = coords(valid, :);
        depths = depths(valid);
        if numel(depths) < 2, return; end
        [depths, sOrd] = sort(depths);
        coords = coords(sOrd, :);
        % Average duplicate-depth entries so interp1 gets unique knots
        [depths, ~, uIdx] = unique(depths, 'stable');
        nU = numel(depths);
        coordsU = zeros(nU, 3);
        for ui = 1:nU
            coordsU(ui,:) = mean(coords(uIdx == ui,:), 1);
        end
        coords = coordsU;
        if numel(depths) < 2, return; end
        % Clamp to track range then interpolate
        bndDepths = max(depths(1), min(depths(end), bndDepths(:)));
        ap_um = interp1(depths, coords(:,1), bndDepths, 'linear');
        dv_um = interp1(depths, coords(:,2), bndDepths, 'linear');
        ml_um = interp1(depths, coords(:,3), bndDepths, 'linear');
        if useHist && ~isempty(hV)
            vs   = hV.voxSize_um;
            sz   = size(hV.vol);      % [DV ML AP]
            ap   = max(1, min(sz(3), round(ap_um / vs)));
            dv   = max(1, min(sz(1), round(dv_um / vs)));
            ml_r = max(1, min(sz(2), round(ml_um / vs)));
            ml   = sz(2) - ml_r + 1;
        else
            ap   = max(1, min(avSz_(1), round(ap_um / 10)));
            dv   = max(1, min(avSz_(2), round(dv_um / 10)));
            ml_r = max(1, min(avSz_(3), round(ml_um / 10)));
            ml   = avSz_(3) - ml_r + 1;
        end
        bm.ap = ap;  bm.dv = dv;  bm.ml = ml;
    end

    %% ── Track overlay voxel helper ───────────────────────────────────────────
    function tc = buildTrackVox(D2A_in, useHist, hV, avSz_)
        tc     = [];
        coords = D2A_in.Coordinates;
        valid  = ~any(isnan(coords), 2);
        coords = coords(valid, :);
        if isempty(coords), return; end

        % Colour each scatter point with the exact same representative colour that
        % drawProbeBar paints onto the corresponding filled block.
        % Method: for each boundary segment, find the first (shallowest) data point
        % in the FULL shank depth sequence (identical to what drawProbeBar uses) and
        % assign its colour to every valid track point whose depth falls in that segment.
        allDepths_t  = D2A_in.Depth(valid);
        shk_         = D2A_in.Shank(1);          % all rows same shank
        si_          = find(allShanks == shk_, 1);
        fullIdx_     = find(D2A.Shank == shk_);  % all rows for this shank in D2A
        [~, fOrd_]   = sort(D2A.Depth(fullIdx_));
        fullDepths_  = D2A.Depth(fullIdx_(fOrd_));
        fullColors_  = D2A.Color(fullIdx_(fOrd_));
        segEdges_    = [-Inf; boundaries{si_}.depths(:); Inf];
        hexcols      = D2A_in.Color(valid);      % fallback: use stored colour
        for seg_ = 1:numel(segEdges_)-1
            lo_ = segEdges_(seg_);
            hi_ = segEdges_(seg_+1);
            inFull = fullDepths_ >= lo_ & fullDepths_ < hi_;
            if ~any(inFull), continue; end
            segCol = fullColors_{find(inFull, 1)};   % same colour drawProbeBar uses
            hexcols(allDepths_t >= lo_ & allDepths_t < hi_) = {segCol};
        end

        if useHist && ~isempty(hV)
            vs   = hV.voxSize_um;
            sz   = size(hV.vol);      % [DV ML AP]
            ap   = max(1, min(sz(3), round(coords(:,1)/vs)));
            dv   = max(1, min(sz(1), round(coords(:,2)/vs)));
            ml_r = max(1, min(sz(2), round(coords(:,3)/vs)));
            ml   = sz(2) - ml_r + 1; % flip to L→R
        else
            ap   = max(1, min(avSz_(1), round(coords(:,1)/10)));
            dv   = max(1, min(avSz_(2), round(coords(:,2)/10)));
            ml_r = max(1, min(avSz_(3), round(coords(:,3)/10)));
            ml   = avSz_(3) - ml_r + 1;
        end

        nPts   = size(coords,1);
        colRGB = zeros(nPts,3);
        for pi = 1:nPts
            colRGB(pi,:) = hex2rgbSafe(hexcols{pi});
        end
        tc.ap = ap; tc.dv = dv; tc.ml = ml; tc.colors = colRGB;
    end

    %% ── Allen CCF parent-area lookup ─────────────────────────────────────────
    % Uses `st` from the enclosing workspace (empty → returns area unchanged).
    % Mirrors getParentArea() in alignatlasdata_automated.m: one level up the
    % structure_id_path hierarchy, matching what the automated alignment uses
    % to define editable boundaries.
    function pa = getParentAreaLocal(area)
        if isempty(st)
            pa = area; return;
        end
        idx = find(strcmpi(st.acronym, area), 1);
        if isempty(idx)
            pa = area; return;
        end
        pathParts = strsplit(st.structure_id_path{idx}, '/');
        if numel(pathParts) > 3
            parentID  = str2double(pathParts{end-2});
            parentIdx = find(st.id == parentID, 1);
            if ~isempty(parentIdx)
                pa = st.acronym{parentIdx};
            else
                pa = area;
            end
        else
            pa = area;
        end
    end

end  % end of main function

%% ════════════════════════════════════════════════════════════════════════════
%  LOCAL (non-nested) helpers  — no access to GUI state
%% ════════════════════════════════════════════════════════════════════════════

function hVol = local_loadHistoVol(mouse, histoRoot)
% Load brainreg downsampled_standard.tiff for this mouse.
hVol = [];
if isempty(mouse), return; end
candidates = {'downsampled_standard.tiff','downsampled_standard.tif', ...
              'downsampled_standard_free.tiff','downsampled_standard_free.tif'};
tiffPath = '';
for ci = 1:numel(candidates)
    p = fullfile(histoRoot, mouse, 'brainreg_output', candidates{ci});
    if exist(p,'file'), tiffPath = p; break; end
end
if isempty(tiffPath)
    hits = dir(fullfile(histoRoot, mouse, 'brainreg_output', 'downsampled_standard*.tif*'));
    if ~isempty(hits)
        ok   = cellfun(@(n) ~contains(n,'additional_channel'), {hits.name});
        hits = hits(ok);
    end
    if ~isempty(hits), tiffPath = fullfile(hits(1).folder, hits(1).name); end
end
if isempty(tiffPath)
    fprintf('AlignmentCurationGUI: no brainreg TIFF for mouse "%s" under %s\n', ...
        mouse, histoRoot);
    return
end
try
    info   = imfinfo(tiffPath);
    nPg    = numel(info);
    h_     = info(1).Height;
    w_     = info(1).Width;
    vol    = local_loadTiffStack(tiffPath, h_, w_, nPg);
    hVol.vol        = vol;
    hVol.volR       = vol;
    hVol.volG       = [];
    hVol.volB       = [];
    hVol.voxSize_um = 25;
    hVol.source     = 'histology';
    hVol.hasRGB     = false;
    fprintf('AlignmentCurationGUI: loaded histology [DV%d ML%d AP%d] at 25µm\n', ...
        size(vol,1), size(vol,2), size(vol,3));
catch ME
    fprintf('AlignmentCurationGUI: could not load TIFF for %s: %s\n', mouse, ME.message);
end
end

function img = local_makeSlice(hVol, plane, idx)
switch plane
    case 'coronal'
        R = local_norm(fliplr(squeeze(hVol.volR(:,:,idx))));
    case 'sagittal'
        R = local_norm(squeeze(hVol.volR(:,idx,:)));
    case 'transverse'
        R = local_norm(fliplr(squeeze(hVol.volR(idx,:,:))'));
end
if hVol.hasRGB
    switch plane
        case 'coronal'
            G = local_norm(fliplr(squeeze(hVol.volG(:,:,idx))));
            B = local_norm(fliplr(squeeze(hVol.volB(:,:,idx))));
        case 'sagittal'
            G = local_norm(squeeze(hVol.volG(:,idx,:)));
            B = local_norm(squeeze(hVol.volB(:,idx,:)));
        case 'transverse'
            G = local_norm(fliplr(squeeze(hVol.volG(idx,:,:))'));
            B = local_norm(fliplr(squeeze(hVol.volB(idx,:,:))'));
    end
    img = cat(3, R, G, B);
else
    img = R;
end
end

function img = local_norm(raw)
raw = double(raw);
lo  = prctile(raw(:), 1);
hi  = prctile(raw(:), 99);
if hi > lo
    img = max(0, min(1, (raw-lo)/(hi-lo)));
else
    img = zeros(size(raw));
end
end

function vol = local_loadTiffStack(path, h, w, nPg)
if exist('tiffreadVolume','builtin') || exist('tiffreadVolume','file')
    vol = tiffreadVolume(path);
else
    vol = zeros(h, w, nPg, 'uint16');
    for zi = 1:nPg
        vol(:,:,zi) = imread(path, zi);
    end
end
end

function hex = local_getHex(st, ri)
hex = '';
try
    if iscell(st.color_hex_triplet)
        raw = st.color_hex_triplet{ri};
    else
        raw = char(st.color_hex_triplet(ri));
    end
    hex = strrep(strtrim(raw), '#', '');
catch
end
end

function rgb = hex2rgbSafe(hex)
hex = strrep(char(hex), '#', '');
if numel(hex) >= 6
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] / 255;
else
    rgb = [0.5 0.5 0.5];
end
end

function val = getfield_safe(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end

function ptr = conditional_ptr(condition)
% Return a valid uifigure Pointer string based on a boolean condition.
if condition
    ptr = 'hand';
else
    ptr = 'arrow';
end
end
