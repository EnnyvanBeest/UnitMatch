function run_EMD_batch_onMerged()
% Runs Yuan et al.'s own, unmodified EMD pipeline (NT_main / EMD_unit_match,
% from github.com/janelia-TDHarrisLab/Yuan-Neuron_Tracking) on every session
% pair of every merged-dataset group that run_deepunitmatch_batch_onMerged.py
% already runs UMPy/DeepUnitMatch on.
%
% This is step 2 of a 3-step pipeline -- see the docstring at the top of
% UnitMatchPy/PaperAnalyses/run_emd_batch_onMerged.py for the full picture:
%   1. `python run_emd_batch_onMerged.py --stage`   (writes manifest.json +
%      Bombcellgood.mat stubs under BASE_OUTPUT/<X>/EMD/_stage/)
%   2. this script                                  (writes Output.mat under
%      BASE_OUTPUT/<X>/EMD/result_<i>_<j>/)
%   3. `python run_emd_batch_onMerged.py --aggregate` (rebuilds a final_matches
%      matrix from all pairs and computes the same functional-score AUCs used
%      for UMPy/DeepUnitMatch)
%
% Nothing in the cloned Neuron_Tracking repo is modified. create_EMD_input.m
% hardcodes bUseKSlabel=0 / bUseBombcellLabel=1, so per-session units are
% still filtered through a Bombcellgood.mat -- but since the manifest already
% only lists units _prepare_session() would keep, that file's unitType is
% trivially all ones (no additional unit is excluded here).
%
% EDIT these two paths if the merged-data layout ever moves; they must match
% BASE_INPUT / BASE_OUTPUT in run_deepunitmatch_batch_onMerged.py exactly.
BASE_INPUT  = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2_merged\merged_data_v2';
BASE_OUTPUT = '\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026_V3_OnMergedData';

% EDIT these to match your local clones if different.
NEURON_TRACKING_REPO = 'C:\Users\EnnyB\Documents\GitHub\Neuron_Tracking';
NPY_MATLAB_REPO = 'C:\Users\EnnyB\Documents\GitHub\npy-matlab';

addpath(genpath(NEURON_TRACKING_REPO));
addpath(genpath(NPY_MATLAB_REPO));

% A pair is redone when its Output.mat is missing, older than this group's
% own manifest.json (the merged data changed upstream -- e.g.
% generate_merged_dataset.py re-merged after a DUM fix changed which units
% get merged), or older than REDO_FROM_DATE (an EMD-side algorithm fix, e.g.
% the length(f1)/size(f1,1) bug fix in EMD_unit_match.m, landed after that
% point). Set REDO_FROM_DATE to NaT to disable the date-based check (freshness
% vs manifest.json still applies); a far-future date forces redo of everything.
REDO_FROM_DATE = datetime(2026, 7, 22, 19, 0, 0);

[groups, parents] = find_merged_groups(BASE_INPUT);
fprintf('Found %d group(s) under %s\n', numel(groups), BASE_INPUT);

for gi = 1:numel(groups)
    merged_dir = groups{gi};
    % parent is found(i).folder from find_merged_groups -- i.e. dir()'s own,
    % already-correct containing-folder value. Deriving this instead via
    % fileparts(merged_dir) is fragile (observed to misbehave on these deeply
    % nested UNC paths, silently leaving 'DeepUnitMatch' in the "parent" and
    % therefore in emd_dir) so we avoid string surgery on merged_dir entirely.
    parent = parents{gi};
    fprintf('\n[%d/%d] %s\n', gi, numel(groups), merged_dir);

    if ~startsWith(parent, BASE_INPUT)
        fprintf('  WARNING: %s is not under BASE_INPUT, skipping.\n', parent);
        continue
    end
    if strcmp(parent, BASE_INPUT)
        emd_dir = fullfile(BASE_OUTPUT, 'EMD');
    else
        subfolder = extractAfter(parent, [BASE_INPUT filesep]);
        emd_dir = fullfile(BASE_OUTPUT, subfolder, 'EMD');
    end
    manifest_path = fullfile(emd_dir, '_stage', 'manifest.json');

    if ~isfile(manifest_path)
        fprintf('  No manifest at %s -- run "python run_emd_batch_onMerged.py --stage" first. Skipping.\n', manifest_path);
        continue
    end

    try
        run_group(merged_dir, emd_dir, manifest_path, REDO_FROM_DATE);
    catch ME
        fprintf('  Group FAILED: %s\n', ME.message);
        disp(getReport(ME));
    end
end

fprintf('\nAll done.\n');
end


function run_group(merged_dir, emd_dir, manifest_path, REDO_FROM_DATE)
manifest = jsondecode(fileread(manifest_path));
sessions = manifest.sessions;
nSessions = numel(sessions);
if nSessions < 2
    fprintf('  Fewer than 2 sessions in manifest, skipping.\n');
    return
end

manifest_info = dir(manifest_path);
manifest_datetime = datetime(manifest_info.datenum, 'ConvertFrom', 'datenum');

stage_root = fullfile(emd_dir, '_stage');

% ---- load mean waveforms + cluster id lists for every session up front ----
mwf = cell(nSessions, 1);
clusterIds = cell(nSessions, 1);
for s = 1:nSessions
    folder = sessions(s).folder;
    ks_dir = fullfile(merged_dir, folder);
    wf_dir = find_rawwaveforms_dir(ks_dir);
    ids = double(sessions(s).cluster_ids(:));
    nUnit = numel(ids);

    mw = [];
    for u = 1:nUnit
        raw = readNPY(fullfile(wf_dir, sprintf('Unit%d_RawSpikes.npy', ids(u))));
        % raw: (spike_width, nChan, cv) -> mean over cv -> (spike_width, nChan) -> (nChan, spike_width)
        raw2d = squeeze(mean(raw, 3))';
        if isempty(mw)
            mw = zeros(nUnit, size(raw2d, 1), size(raw2d, 2));
        end
        mw(u, :, :) = raw2d;
    end
    mwf{s} = mw;
    clusterIds{s} = ids;

    % trivial Bombcellgood.mat stub: every listed unit is already good.
    sess_stage_dir = fullfile(stage_root, folder);
    if ~isfolder(sess_stage_dir), mkdir(sess_stage_dir); end
    unitType = ones(nUnit, 1); %#ok<NASGU>
    save(fullfile(sess_stage_dir, 'Bombcellgood.mat'), 'unitType');
end

% ---- channel positions + probe geometry, derived from the first session ----
chan_pos_raw = readNPY(fullfile(merged_dir, sessions(1).folder, 'channel_positions.npy'));
if size(chan_pos_raw, 2) >= 3
    chan_pos_emd = double(chan_pos_raw(:, 2:3));
else
    chan_pos_emd = double(chan_pos_raw(:, 1:2));
end
ux = sort(unique(chan_pos_emd(:, 1)));
uz = sort(unique(chan_pos_emd(:, 2)));
xStep = median(diff(ux));
zStep = median(diff(uz));
if isempty(xStep) || xStep == 0, xStep = 32; end
if isempty(zStep) || zStep == 0, zStep = 15; end
fprintf('  Probe geometry: xStep=%.1f, zStep=%.1f (derived from channel_positions.npy)\n', xStep, zStep);

ts = size(mwf{1}, 3);

% ---- every session pair, mirroring the exhaustive r1<r2 loop the Python ----
% ---- batch script uses for UMPy/DeepUnitMatch                          ----
for i = 1:nSessions - 1
    for j = i + 1:nSessions
        folder1 = sessions(i).folder;
        folder2 = sessions(j).folder;
        result_dir = fullfile(emd_dir, sprintf('result_%s_%s', folder1, folder2));
        out_mat = fullfile(result_dir, 'Output.mat');
        if isfile(out_mat)
            out_info = dir(out_mat);
            out_datetime = datetime(out_info.datenum, 'ConvertFrom', 'datenum');
            fresh_vs_manifest = out_datetime >= manifest_datetime;
            fresh_vs_date = isnat(REDO_FROM_DATE) || out_datetime >= REDO_FROM_DATE;
            if fresh_vs_manifest && fresh_vs_date
                fprintf('  Skipping %s vs %s (exists and fresh)\n', folder1, folder2);
                continue
            end
            fprintf('  Re-matching %s vs %s (existing Output.mat is stale)...\n', folder1, folder2);
        else
            fprintf('  Matching %s vs %s ...\n', folder1, folder2);
        end

        clear input
        input.input_path = stage_root;
        input.data_path1 = folder1;
        input.data_path2 = folder2;
        input.chan_pos = chan_pos_emd;
        input.chan_map = (0:size(chan_pos_emd, 1) - 1)'; % unused: nChanPos == nChanMW here
        input.EMD_path = fullfile(result_dir, 'EMD_intermediate');
        input.KSLabel_name = 'cluster_KSLabel.tsv'; % unused (bUseKSlabel=0 in create_EMD_input)
        input.shank = -1;
        input.fs = 30000;
        input.ts = ts;
        input.l2_weights = 1500;
        input.threshold = 10; % z distance threshold for matched units, um
        input.validation = 0;
        input.xStep = xStep;
        input.zStep = zStep;
        input.dim_mask = logical([1, 1, 1, 0, 0, 0, 0, 0, 0, 1]);
        input.dim_mask_physical = logical([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]);
        input.dim_mask_wf = logical([0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
        input.diagDistCalc = false;
        input.result_path = result_dir;
        input.input_name = 'input_pre.mat';
        input.input_name_post = 'input_post.mat';
        input.filename_pre = 'EMD_pre.mat';
        input.filename_post = 'EMD_post.mat';

        try
            NT_main(input, mwf{i}, mwf{j});

            % Stash the per-session cluster-id lists (parallel to mwf{i}/mwf{j}
            % rows) alongside Yuan's own Output.mat, so the Python aggregator
            % can translate f1_labels/f2_labels (positions into mwf{i}/mwf{j})
            % back to actual cluster IDs.
            s_out = load(out_mat, 'output');
            output = s_out.output;
            output.cluster_ids_1 = clusterIds{i};
            output.cluster_ids_2 = clusterIds{j};
            save(out_mat, 'output');
        catch ME
            fprintf('  FAILED %s vs %s: %s\n', folder1, folder2, ME.message);
            disp(getReport(ME));
        end
    end
end
end


function [groups, parents] = find_merged_groups(baseInput)
% Mirrors find_merged_groups() in run_deepunitmatch_batch_onMerged.py: every
% 'DeepUnitMatch' folder under baseInput whose immediate children are
% numbered session folders (0, 1, 2, ...).
%
% This deliberately does NOT use dir(fullfile(baseInput,'**','DeepUnitMatch')):
% when the final path segment is a literal name (no wildcard) that resolves
% to an existing directory, dir() lists that folder's *contents* -- including
% its own '.' self-entry, which passes an [entries.isdir] filter -- rather
% than matching 'DeepUnitMatch' as a name against each parent's listing. That
% silently produced e.g. '...\1\DeepUnitMatch\.' as a "found" path. A manual
% directory-by-directory walk (mirroring Python's os.walk here) avoids this
% class of bug entirely, since every dir() call below is either on a genuine
% wildcard-free concrete folder (to list its children) or filtered by name
% match against those children -- never both conflated at once.
groups = {};
parents = {};
walk_for_groups(baseInput);

    function walk_for_groups(currentDir)
        entries = dir(currentDir);
        entries = entries([entries.isdir] & ~ismember({entries.name}, {'.', '..'}));
        for k = 1:numel(entries)
            childName = entries(k).name;
            childPath = fullfile(currentDir, childName);
            if strcmp(childName, 'DeepUnitMatch')
                sub = dir(childPath);
                sub = sub([sub.isdir] & ~ismember({sub.name}, {'.', '..'}));
                if ~isempty(sub) && any(~isnan(str2double({sub.name})))
                    groups{end + 1} = childPath; %#ok<AGROW>
                    parents{end + 1} = currentDir; %#ok<AGROW>
                end
                % DeepUnitMatch group folders only ever contain numeric
                % session subfolders -- nothing more to find by recursing in.
                continue
            end
            walk_for_groups(childPath);
        end
    end
end


function p = find_rawwaveforms_dir(ks_dir)
% Mirrors the three RawWaveforms locations util.paths_from_KS() checks.
candidates = {
    fullfile(ks_dir, 'RawWaveforms')
    fullfile(ks_dir, 'qMetrics', 'RawWaveforms')
    fullfile(ks_dir, 'bombcell', 'RawWaveforms')
    };
for k = 1:numel(candidates)
    if isfolder(candidates{k})
        p = candidates{k};
        return
    end
end
error('Could not find RawWaveforms folder under %s', ks_dir);
end
