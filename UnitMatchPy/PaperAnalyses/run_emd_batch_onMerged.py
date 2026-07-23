# Bolt-on batch driver that gets Yuan et al. EMD (Earth Mover's Distance,
# github.com/janelia-TDHarrisLab/Yuan-Neuron_Tracking) matches onto the same
# merged-dataset groups/session-pairs that run_deepunitmatch_batch_onMerged.py
# already runs UMPy/DeepUnitMatch on, so EMD can be compared to both using the
# exact same functional-score AUC framework (test.AUC over a boolean
# final_matches matrix), instead of re-deriving UnitMatch.mat-based AUCs the
# way MATLAB/Paper_Figures/ComparetoYuanetal.m does.
#
# EMD itself is MATLAB-only and is not reimplemented here. This script only
# prepares inputs for, and aggregates outputs from, an unmodified run of
# Yuan's own NT_main/EMD_unit_match pipeline. Usage is three steps:
#
#   1. python run_emd_batch_onMerged.py --stage
#      For every merged group (same discovery as find_merged_groups()), works
#      out -- without loading any waveform data -- exactly which good units
#      _prepare_session() would use for each session (same label + same
#      RawSpikes-file-exists checks), and writes:
#        BASE_OUTPUT/<X>/EMD/_stage/manifest.json
#        BASE_OUTPUT/<X>/EMD/_stage/<session>/Bombcellgood.mat
#      manifest.json records, per session and in the same order/positions
#      _prepare_session()'s good_units arrays use, the list of cluster IDs to
#      feed EMD. Bombcellgood.mat is a trivial unitType=1 stub (all units are
#      already pre-filtered to good) -- it exists purely so the unmodified
#      Yuan create_EMD_input.m (which hardcodes bUseBombcellLabel=1) has
#      something to load; no unit is actually excluded by it here.
#
#   2. Run the MATLAB driver (calls Yuan's own NT_main/EMD_unit_match
#      unmodified) for every group and every session pair:
#          matlab -batch "run_EMD_batch_onMerged"
#      See MATLAB/Paper_Figures/EMD_integration/run_EMD_batch_onMerged.m.
#      This reads RawSpikes.npy/channel_positions.npy directly and writes
#      Yuan-format Output.mat (plus the original cluster IDs, for step 3) to
#        BASE_OUTPUT/<X>/EMD/result_<i>_<j>/Output.mat
#
#   3. python run_emd_batch_onMerged.py --aggregate
#      Reloads every pair's Output.mat, translates EMD's thresholded matches
#      (output.results_wth) from mwf-array positions back to cluster IDs and
#      then into the (session_id, good_units) coordinate system UMPy/DUM use,
#      assembles one boolean final_matches matrix per group, computes the
#      same five functional-score AUCs, and saves MatchTable.csv /
#      MatchingOverview.png to BASE_OUTPUT/<X>/EMD/ -- a sibling of
#      DeepUnitMatch/ and UMPy/, ready for the same post-hoc AUC comparisons.

import os
import sys
import json
import datetime
import argparse
import traceback

import numpy as np
import scipy.io as sio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))

import batch_lock
import UnitMatchPy.utils as util
import UnitMatchPy.overlord as ov
import UnitMatchPy.save_utils as su
from DeepUnitMatch.testing import test

from run_deepunitmatch_batch_onMerged import (
    BASE_INPUT,
    BASE_OUTPUT,
    find_merged_groups,
    build_good_units_from_labels,
    _prepare_session,
)

# See batch_lock.sentinel_is_fresh() / run_deepunitmatch_batch_onMerged.py's
# REDO_FROM_DATE for what these do. Split into two separate cutoffs because
# --stage and --aggregate depend on different things going stale:
#   STAGE_REDO_FROM_DATE: manifest.json/Bombcellgood.mat depend on the merged
#     tree's good-unit population (cluster_bc_unitType.tsv, RawSpikes.npy
#     presence), which changes if generate_merged_dataset.py re-merges
#     upstream (e.g. after a DUM fix changes which units get merged).
#   AGGREGATE_REDO_FROM_DATE: the saved AUCs/MatchTable.csv depend on
#     test.py's functional-score functions (ISI_correlations, FR_diff,
#     ISI_CV_diff, refpop_correlations, natim_correlations) and test.AUC
#     itself, independently of whether the underlying EMD matches changed --
#     e.g. the FR_diff/ISI_CV_diff NaN-vs-0 fallback fix.
# None falls back to plain "skip if present"; a far-future date forces redo.
STAGE_REDO_FROM_DATE = datetime.datetime(2026, 7, 22, 19, 0, 0)
AGGREGATE_REDO_FROM_DATE = datetime.datetime(2026, 7, 22, 19, 0, 0)


# ── path helpers ─────────────────────────────────────────────────────────────


def get_emd_dir(merged_dir):
    """Return the EMD output directory for a given merged-data group dir."""
    subfolder = os.path.relpath(os.path.dirname(merged_dir), BASE_INPUT)
    return os.path.join(BASE_OUTPUT, subfolder, "EMD")


def get_stage_dir(merged_dir):
    return os.path.join(get_emd_dir(merged_dir), "_stage")


def emd_results_exist(merged_dir):
    """Return True when the EMD sentinel output file is present and fresh (see AGGREGATE_REDO_FROM_DATE)."""
    sentinel = os.path.join(get_emd_dir(merged_dir), "MatchingOverview.png")
    return batch_lock.sentinel_is_fresh(sentinel, AGGREGATE_REDO_FROM_DATE)


# ── stage phase ──────────────────────────────────────────────────────────────


def stage_group(merged_dir):
    """
    Work out, per session, the exact set/order of good units _prepare_session()
    would use -- same GOOD/NON-SOMA GOOD label check plus "RawSpikes file
    actually exists" check as build_good_units_from_labels() /
    load_waveforms_for_good_units() -- without loading the (large) waveform
    arrays themselves; the MATLAB driver reads RawSpikes.npy directly.

    Writes manifest.json (session order/cluster-ID lists) plus a trivial
    Bombcellgood.mat (unitType=1 for every listed unit) per session.
    """
    session_names = sorted(
        (
            d
            for d in os.listdir(merged_dir)
            if d.isdigit() and os.path.isdir(os.path.join(merged_dir, d))
        ),
        key=int,
    )
    if not session_names:
        print(f"  ERROR: no numbered session folders found under {merged_dir}")
        return None

    ks_dirs = [os.path.join(merged_dir, d) for d in session_names]
    try:
        wave_paths, unit_label_paths, _ = util.paths_from_KS(ks_dirs)
    except Exception as e:
        print(f"  ERROR in paths_from_KS: {e}")
        traceback.print_exc()
        return None

    good_units_per_session = build_good_units_from_labels(unit_label_paths)

    stage_dir = get_stage_dir(merged_dir)
    os.makedirs(stage_dir, exist_ok=True)

    sessions_manifest = []
    for folder, wave_path, g_units in zip(
        session_names, wave_paths, good_units_per_session
    ):
        kept_ids = [
            int(uid)
            for uid in g_units.flatten()
            if os.path.exists(os.path.join(wave_path, f"Unit{int(uid)}_RawSpikes.npy"))
        ]
        if not kept_ids:
            print(f"  Session {folder}: no loadable good units, dropping from manifest.")
            continue

        sess_stage_dir = os.path.join(stage_dir, folder)
        os.makedirs(sess_stage_dir, exist_ok=True)
        sio.savemat(
            os.path.join(sess_stage_dir, "Bombcellgood.mat"),
            {"unitType": np.ones((len(kept_ids), 1), dtype=float)},
        )
        sessions_manifest.append({"folder": folder, "cluster_ids": kept_ids})
        print(f"    [{folder}] {len(kept_ids)} good unit(s)")

    if len(sessions_manifest) < 2:
        print(f"  Fewer than 2 usable sessions, skipping group: {merged_dir}")
        return None

    manifest = {"merged_dir": merged_dir, "sessions": sessions_manifest}
    with open(os.path.join(stage_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Staged {len(sessions_manifest)} session(s) -> {stage_dir}")
    return stage_dir


def run_stage():
    print(f"Scanning for merged-data groups under:\n  {BASE_INPUT}\n")
    groups = find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")
        stage_path = os.path.join(get_stage_dir(merged_dir), "manifest.json")
        if batch_lock.sentinel_is_fresh(stage_path, STAGE_REDO_FROM_DATE):
            print(f"  Skipping (manifest exists and is fresh): {stage_path}")
            continue
        try:
            stage_group(merged_dir)
        except Exception as e:
            print(f"  Staging FAILED: {e}")
            traceback.print_exc()

    print(
        "\nStaging done. Now run the MATLAB driver, e.g.:\n"
        '  matlab -batch "run_EMD_batch_onMerged"\n'
        "(see MATLAB/Paper_Figures/EMD_integration/run_EMD_batch_onMerged.m),\n"
        "then re-run this script with --aggregate."
    )


# ── aggregate phase ──────────────────────────────────────────────────────────


def aggregate_group(merged_dir):
    emd_dir = get_emd_dir(merged_dir)
    manifest_path = os.path.join(get_stage_dir(merged_dir), "manifest.json")
    if not os.path.isfile(manifest_path):
        print(f"  No manifest found (run --stage first): {manifest_path}")
        return
    with open(manifest_path) as f:
        manifest = json.load(f)
    manifest_sessions = manifest["sessions"]

    sess = _prepare_session(merged_dir)
    if sess is None:
        print(f"  _prepare_session failed for {merged_dir}")
        return

    param = sess["param"]
    session_id = sess["session_id"]
    session_switch = sess["session_switch"]
    good_units = sess["good_units"]
    channel_pos = sess["channel_pos"]
    waveform = sess["waveform"]
    n_total = waveform.shape[0]

    if len(manifest_sessions) != len(good_units):
        print(
            f"  WARNING: manifest has {len(manifest_sessions)} session(s) but "
            f"_prepare_session loaded {len(good_units)}; re-run --stage. Skipping."
        )
        return

    # cluster_id -> local index within this session's good_units, per session
    id_to_pos = [
        {int(cid): pos for pos, cid in enumerate(gu.flatten())} for gu in good_units
    ]

    final_matches = np.zeros((n_total, n_total), dtype=bool)
    n_pairs_total = len(manifest_sessions) * (len(manifest_sessions) - 1) // 2
    n_pairs_loaded = 0
    n_rows_skipped = 0

    for r1 in range(len(manifest_sessions)):
        for r2 in range(r1 + 1, len(manifest_sessions)):
            folder1 = manifest_sessions[r1]["folder"]
            folder2 = manifest_sessions[r2]["folder"]
            out_path = os.path.join(
                emd_dir, f"result_{folder1}_{folder2}", "Output.mat"
            )
            if not os.path.isfile(out_path):
                print(f"    Missing EMD result (run the MATLAB step): {out_path}")
                continue

            mat = sio.loadmat(out_path, squeeze_me=True, struct_as_record=False)
            output = mat["output"]
            results = np.atleast_2d(getattr(output, "results_wth", np.empty((0, 7))))
            n_pairs_loaded += 1
            if results.size == 0:
                continue

            cluster_ids_1 = np.atleast_1d(output.cluster_ids_1).astype(int)
            cluster_ids_2 = np.atleast_1d(output.cluster_ids_2).astype(int)

            for row in results:
                pos2_in_all = int(row[1]) - 1  # col2: f2 label (1-based position in mwf2)
                pos1_in_all = int(row[2]) - 1  # col3: f1 label (1-based position in mwf1)
                cid1 = int(cluster_ids_1[pos1_in_all])
                cid2 = int(cluster_ids_2[pos2_in_all])
                if cid1 not in id_to_pos[r1] or cid2 not in id_to_pos[r2]:
                    n_rows_skipped += 1
                    continue
                g1 = session_switch[r1] + id_to_pos[r1][cid1]
                g2 = session_switch[r2] + id_to_pos[r2][cid2]
                final_matches[g1, g2] = True
                final_matches[g2, g1] = True

    if n_rows_skipped:
        print(f"    ({n_rows_skipped} matched row(s) could not be mapped back to good_units)")

    n_matches = int(np.sum(final_matches)) // 2
    print(f"  {n_matches} EMD match(es) assembled ({n_pairs_loaded}/{n_pairs_total} pairs loaded)")

    functional_scores = {}
    try:
        isicorr = test.ISI_correlations(param)
        auc_isi = test.AUC(final_matches, isicorr, session_id)
        print(f"AUC (ISI correlations):            {auc_isi:.3f}")
        functional_scores["ISI_correlations"] = isicorr

        isikl = test.ISI_KL_divergence(param)
        auc_isikl = test.AUC(final_matches, -isikl, session_id)
        print(f"AUC (ISI KL divergence):           {auc_isikl:.3f}")
        functional_scores["ISI_KL_divergence"] = isikl

        isiwass = test.ISI_wasserstein_distance(param)
        auc_isiwass = test.AUC(final_matches, -isiwass, session_id)
        print(f"AUC (ISI Wasserstein distance):    {auc_isiwass:.3f}")
        functional_scores["ISI_wasserstein_distance"] = isiwass

        refpopcorr = test.refpop_correlations(param, matches=final_matches)
        auc_refpop = test.AUC(final_matches, refpopcorr, session_id)
        print(f"AUC (ref. pop. correlation):       {auc_refpop:.3f}")
        functional_scores["refpop_correlations"] = refpopcorr

        frdiff = test.FR_diff(param)
        auc_fr = test.AUC(final_matches, -frdiff, session_id)
        print(f"AUC (firing rate difference):      {auc_fr:.3f}")
        functional_scores["FR_diff"] = frdiff

        cvdiff = test.ISI_CV_diff(param)
        auc_cv = test.AUC(final_matches, -cvdiff, session_id)
        print(f"AUC (ISI CV difference):           {auc_cv:.3f}")
        functional_scores["ISI_CV_diff"] = cvdiff

        try:
            natimcorr = test.natim_correlations(param)
            auc_natim = test.AUC(final_matches, natimcorr, session_id)
            print(f"AUC (nat. image correlations):     {auc_natim:.3f}")
            functional_scores["natim_correlations"] = natimcorr
        except Exception:
            pass  # natim data may not be available for every session
    except Exception as e:
        print(f"  WARNING: functional score computation failed: {e}")
        functional_scores = {}

    clus_info = {
        "good_units": good_units,
        "session_switch": session_switch,
        "session_id": session_id,
        "original_ids": np.concatenate(good_units),
    }
    extracted_wave_properties = ov.extract_parameters(
        waveform, channel_pos, clus_info, param
    )

    os.makedirs(emd_dir, exist_ok=True)
    su.save_to_output(
        emd_dir,
        {},
        np.argwhere(final_matches),
        final_matches.astype(float),
        extracted_wave_properties["avg_centroid"],
        extracted_wave_properties["avg_waveform"],
        extracted_wave_properties["avg_waveform_per_tp"],
        extracted_wave_properties["max_site"],
        final_matches.astype(float),
        final_matches,
        clus_info,
        param,
        UIDs=None,
        matches_curated=None,
        save_match_table=True,
        functional_scores=functional_scores if functional_scores else None,
    )
    su.save_auc_summary(
        emd_dir, test.auc_summary_from_functional_scores(functional_scores, final_matches, session_id)
    )

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(final_matches, cmap="viridis", aspect="auto")
    ax.set_title(f"EMD matches (n={n_matches})")
    ax.set_xlabel("Unit")
    ax.set_ylabel("Unit")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(emd_dir, "MatchingOverview.png"), dpi=150)
    plt.close(fig)

    if functional_scores:
        score_meta = {
            "ISI_correlations": ("ISI correlations", "viridis"),
            "ISI_KL_divergence": ("ISI KL divergence", "magma"),
            "ISI_wasserstein_distance": ("ISI Wasserstein distance", "magma"),
            "refpop_correlations": ("Ref. pop. correlations", "viridis"),
            "FR_diff": ("Firing rate difference", "magma"),
            "ISI_CV_diff": ("ISI CV difference", "magma"),
            "natim_correlations": ("Nat. image correlations", "viridis"),
        }
        keys = [k for k in score_meta if k in functional_scores]
        fig, axes = plt.subplots(1, len(keys), figsize=(5 * len(keys), 5))
        if len(keys) == 1:
            axes = [axes]
        for ax, key in zip(axes, keys):
            title, cmap = score_meta[key]
            im = ax.imshow(functional_scores[key], cmap=cmap, aspect="auto")
            ax.set_title(title)
            fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(emd_dir, "FunctionalScores.png"), dpi=150)
        plt.close(fig)

    print(f"  Results saved to: {emd_dir}")


def run_aggregate():
    print(f"Scanning for merged-data groups under:\n  {BASE_INPUT}\n")
    groups = find_merged_groups()
    if not groups:
        print("No merged-data groups found.")
        return
    print(f"Found {len(groups)} group(s).\n")

    for i, merged_dir in enumerate(groups):
        print(f"\n[{i + 1}/{len(groups)}] {merged_dir}")
        if emd_results_exist(merged_dir):
            print(f"  Skipping (results exist and are fresh): {get_emd_dir(merged_dir)}")
            continue
        try:
            aggregate_group(merged_dir)
        except Exception as e:
            print(f"  Aggregation FAILED: {e}")
            traceback.print_exc()

    print("\nAll done.")


# ── entry point ───────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Stage/aggregate EMD matching around the unmodified MATLAB "
        "Yuan et al. pipeline, for the merged-dataset groups used by "
        "run_deepunitmatch_batch_onMerged.py."
    )
    parser.add_argument(
        "--stage", action="store_true", help="Write manifests + Bombcellgood.mat stubs."
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="Reload MATLAB EMD results and compute functional-score AUCs.",
    )
    args = parser.parse_args()

    if not args.stage and not args.aggregate:
        parser.print_help()
        return

    if args.stage:
        run_stage()
    if args.aggregate:
        run_aggregate()


if __name__ == "__main__":
    main()
