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
av = []; avSz = [0 0 0]; atlasCmap = [];
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
%   Cols 9-11 (row 1):     Reset | Flag | Accept buttons
%   Cols 9-11 (row 2):     status label
ctrlG = uigridlayout(outerG, [2, 11], ...
    'ColumnWidth', {'0.9x','3x','0.9x','3x','0.7x','1.2x','0.7x','2.5x','1.4x','1.6x','1.6x'}, ...
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

% Buttons (row 1)
btnReset = uibutton(ctrlG, 'Text','Reset', ...
    'BackgroundColor',[0.45 0.45 0.45], 'FontColor','w', 'ButtonPushedFcn',@onReset);
btnReset.Layout.Row = 1; btnReset.Layout.Column = 9;
btnFlag = uibutton(ctrlG, 'Text','Flag for Review', ...
    'BackgroundColor',[0.80 0.55 0.0], 'FontColor','w', 'ButtonPushedFcn',@onFlag);
btnFlag.Layout.Row = 1; btnFlag.Layout.Column = 10;
btnAccept = uibutton(ctrlG, 'Text','Accept & Save', ...
    'BackgroundColor',[0.10 0.60 0.10], 'FontColor','w', 'ButtonPushedFcn',@onAccept);
btnAccept.Layout.Row = 1; btnAccept.Layout.Column = 11;

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
lblStatus.Layout.Row = 2; lblStatus.Layout.Column = [9 11];

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
            hw   = 5;
            if numel(depths_s) > 1
                hw = median(diff(depths_s))/2 + 0.5;
            end

            % Area-coloured patches
            for di = 1:numel(depths_s)
                c = hex2rgbSafe(hex_s{di});
                patch(axProbe, [x0 x1 x1 x0], depths_s(di)+[-hw -hw hw hw], c, ...
                    'EdgeColor','none', 'HitTest','off', 'PickableParts','none');
            end

            % Shank label above column
            text(axProbe, xMid, max(depths_s)+hw*2.5, sprintf('S%d',shk), ...
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
                % Area label on right edge
                text(axProbe, x1+0.01, bnd.depths(bi), ...
                    sprintf('%s  |  %s', bnd.areaBelow{bi}, bnd.areaAbove{bi}), ...
                    'FontSize', 6, 'VerticalAlignment','middle', ...
                    'Clipping','on', 'HitTest','off');
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
        tc = buildTrackVox(D2A(shankMask,:), useHistology, hVol, avSz);

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
        axis(axCoronal,'image');
        set(axCoronal,'YDir','reverse','XTickLabel',{},'YTickLabel',{});
        xlabel(axCoronal,'ML (L→R)'); ylabel(axCoronal,'DV');

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
        axis(axSagittal,'image');
        set(axSagittal,'YDir','reverse','XTickLabel',{},'YTickLabel',{});
        xlabel(axSagittal,'AP'); ylabel(axSagittal,'DV');
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
        figPx   = fig.CurrentPoint;      % [x y] in figure pixels (from bottom-left)
        % Pixel coordinate relative to axes inner area
        px = figPx(1) - axPos(1);
        py = figPx(2) - axPos(2);
        xl = axProbe.XLim;
        yl = axProbe.YLim;
        cx = xl(1) + px / axPos(3) * diff(xl);
        cy = yl(1) + py / axPos(4) * diff(yl);
        inside = (px >= 0) && (px <= axPos(3)) && (py >= 0) && (py <= axPos(4));
    end

    function onMouseDown(~,~)
        [cx, cy, inside] = figPtToAxProbe();
        if ~inside, return; end

        xl = axProbe.XLim; yl = axProbe.YLim;
        shanks = unique(D2A.Shank);
        nS  = numel(shanks);
        gap = 0.05;
        sw  = (1 - gap*(nS+1)) / nS;
        snapTol = diff(yl) * 0.03;   % 3% of depth range

        for si = 1:nS
            x0 = gap + (si-1)*(sw+gap);
            x1 = x0 + sw;
            if cx < x0 || cx > x1, continue; end
            bnd = boundaries{si};
            if isempty(bnd.depths), continue; end
            [minDist, bi] = min(abs(bnd.depths - cy));
            if minDist < snapTol
                dragState.active      = true;
                dragState.shankIdx    = si;
                dragState.boundaryIdx = bi;
                fig.Pointer           = 'crosshair';
                return;
            end
        end
    end

    function onMouseMove(~,~)
        [cx, cy, inside] = figPtToAxProbe();

        % Update cursor: show crosshair when hovering near any boundary line
        if ~dragState.active
            shanks  = unique(D2A.Shank);
            nS      = numel(shanks);
            gap     = 0.05;
            sw      = (1 - gap*(nS+1)) / nS;
            yl      = axProbe.YLim;
            snapTol = diff(yl) * 0.03;
            nearBnd = false;
            if inside
                for si = 1:nS
                    x0 = gap + (si-1)*(sw+gap);
                    x1 = x0 + sw;
                    if cx < x0 || cx > x1, continue; end
                    bnd = boundaries{si};
                    if ~isempty(bnd.depths) && min(abs(bnd.depths - cy)) < snapTol
                        nearBnd = true; break;
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
        si = dragState.shankIdx;
        bi = dragState.boundaryIdx;
        dragState.active = false;
        fig.Pointer = 'arrow';

        % Commit: reassign area labels around the moved boundary
        D2A        = applyBoundaryMove(D2A, si, bi, boundaries{si}.depths(bi));
        boundaries = computeBoundaries(D2A);
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
            areas = D2A_in.Area(idx);
            trans = find(~strcmp(areas(1:end-1), areas(2:end)));
            b.depths    = (D2A_in.Depth(idx(trans)) + D2A_in.Depth(idx(trans+1))) / 2;
            b.areaBelow = areas(trans);    % lower depth side
            b.areaAbove = areas(trans+1); % higher depth side
            b.rowBelow  = idx(trans);
            b.rowAbove  = idx(trans+1);
            bnd{si_} = b;
        end
    end

    function D2A_out = applyBoundaryMove(D2A_in, si, bi, newBoundaryDepth)
        % Reassign area labels for rows between the two areas flanking boundary bi.
        D2A_out = D2A_in;
        shanks  = unique(D2A_in.Shank);
        shk     = shanks(si);
        bnd_    = boundaries{si};

        areaBelow = bnd_.areaBelow{bi};
        areaAbove = bnd_.areaAbove{bi};
        colBelow  = getAreaColor(D2A_in, areaBelow);
        colAbove  = getAreaColor(D2A_in, areaAbove);

        % Find all rows on this shank that belong to either neighbouring area
        mask   = D2A_in.Shank == shk & ismember(D2A_in.Area, {areaBelow, areaAbove});
        zIdx   = find(mask);
        zDepths = D2A_in.Depth(zIdx);

        isBelow = zDepths <= newBoundaryDepth;
        D2A_out.Area(zIdx( isBelow)) = {areaBelow};
        D2A_out.Color(zIdx( isBelow)) = {colBelow};
        D2A_out.Area(zIdx(~isBelow)) = {areaAbove};
        D2A_out.Color(zIdx(~isBelow)) = {colAbove};
    end

    function col = getAreaColor(D2A_in, area)
        idx = find(strcmp(D2A_in.Area, area), 1);
        if ~isempty(idx)
            col = D2A_in.Color{idx};
        else
            col = 'aaaaaa';
        end
    end

    %% ── Track overlay voxel helper ───────────────────────────────────────────
    function tc = buildTrackVox(D2A_in, useHist, hV, avSz_)
        tc     = [];
        coords = D2A_in.Coordinates;
        valid  = ~any(isnan(coords), 2);
        coords = coords(valid, :);
        if isempty(coords), return; end
        hexcols = D2A_in.Color(valid);

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
