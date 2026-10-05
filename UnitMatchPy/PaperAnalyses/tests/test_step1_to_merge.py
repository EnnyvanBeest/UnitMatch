"""Integration test: step 1 (unified pipeline on non-merged data) -> step 2 (merge) -> merged-data loading.

Run:  python PaperAnalyses/tests/test_step1_to_merge.py        (~5-15 min, mostly network I/O)

Uses one real recording location (AL032/19011111882/1, 2 sessions; raw data and its
UnitMatch.mat are only READ). Everything written goes to a local temp folder (printed at the
start): step-1 output, merged data, logs. Checks:
  1. step-1 session: tagged for the shared pipeline (group, stage, non-merged layout)
  2. DeepUnitMatch + UMPy via the shared pipeline: outputs, MatchTable with unique IDs,
     unique-ID match folders, natural-image scores found in the non-merged layout (if recorded)
  3. log events under stage "step1" with the right location key
  4. generate_merged_dataset merges from that step-1 output
  5. the merged location loads for the merged-data runs (_prepare_session)
Ends with "ALL CHECKS PASSED".
"""
import json
import os
import shutil
import sys
import tempfile

import pandas as pd

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../PaperAnalyses
sys.path.insert(0, PA)

GROUP = "AL032/19011111882/1"


def main():
    import pipeline_config as cfg

    root = os.path.join(tempfile.gettempdir(), "dum_test_step1_to_merge")
    print("test output folder:", root)
    shutil.rmtree(root, ignore_errors=True)
    cfg.UNMERGED_OUTPUT = os.path.join(root, "unmerged")
    cfg.MERGED_DATA = os.path.join(root, "merged")
    cfg.LOG_DIR = os.path.join(root, "logs")

    import pipeline_log as plog
    import run_deepunitmatch_batch as s1
    import run_deepunitmatch_batch_onMerged as onm
    import generate_merged_dataset as gm

    # point the modules' copies of the paths at the temp folder
    s1.BASE_OUTPUT = cfg.UNMERGED_OUTPUT
    gm.DUM_NONMERGED_DATAPATH, gm.MERGED_DATAPATH = cfg.UNMERGED_OUTPUT, cfg.MERGED_DATA
    onm.BASE_INPUT = cfg.MERGED_DATA

    mat = os.path.join(s1.BASE_INPUT, *GROUP.split("/"), "UnitMatch", "UnitMatch.mat")

    # 1. session preparation
    sess = s1._prepare_session(mat)
    assert sess is not None
    assert sess["group"] == GROUP and sess["log_stage"] == "step1" and sess["merged_architecture"] is False
    print(f"1. step-1 session OK: {sess['waveform'].shape[0]} units, {sess['param']['n_sessions']} sessions, "
          f"{sess['n_units_dropped_by_snippets']} dropped by snippets")

    # 2. both methods through the shared pipeline
    model = s1.test.load_trained_model(device=s1.DEVICE)
    s1.run_deep_unit_match(sess, model)
    s1.run_umpy(sess)
    for method, out in (("DeepUnitMatch", s1.get_save_dir(mat)), ("UMPy", s1.get_umpy_save_dir(mat))):
        assert os.path.isfile(os.path.join(out, "MatchingOverview.png")), method
        header = pd.read_csv(os.path.join(out, "MatchTable.csv"), nrows=0).columns
        assert {"UID 1", "UID 2", "UID Cons 1", "UID Cons 2"} <= set(header), method
        for suffix in ("AssignUniqueID", "AssignUniqueID_Conservative"):
            assert os.path.isfile(os.path.join(f"{out}_{suffix}", "AUC_summary.json")), (method, suffix)
        summary = json.load(open(os.path.join(out, "AUC_summary.json")))
        print(f"2. {method}: {summary.get('n_matches_across_sessions')} across-session matches, "
              f"natural-image AUC {'present' if 'natim_correlations' in summary else 'absent (no natural-image data?)'}")

    # 3. logging
    ev = plog.latest_events("step1")
    for method in ("DeepUnitMatch", "UMPy"):
        e = ev.get(("step1", GROUP, method))
        assert e is not None and e["status"] == "done", (method, e)
    print("3. logged under stage step1 for", GROUP)

    # 4. merge from the new step-1 output
    source_dir = s1.get_save_dir(mat)
    gm.run_merging_process([os.path.join(source_dir, "UMparam.pickle")], [source_dir])
    merged_dir = os.path.join(cfg.MERGED_DATA, *GROUP.split("/"), "DeepUnitMatch")
    sessions = sorted(d for d in os.listdir(merged_dir) if d.isdigit())
    assert sessions and all(os.path.isfile(os.path.join(merged_dir, d, gm.MERGE_COMPLETE_MARKER)) for d in sessions)
    merge_ev = plog.latest_events("merge").get(("merge", GROUP, "merge"))
    assert merge_ev and merge_ev["status"] == "done", merge_ev
    print(f"4. merged: {len(sessions)} sessions; {merge_ev['message']}")

    # 5. merged data loads for the merged-data runs
    msess = onm._prepare_session(merged_dir)
    assert msess is not None and msess["param"]["n_sessions"] == sess["param"]["n_sessions"]
    print(f"5. merged location loads: {msess['waveform'].shape[0]} units (step 1 had {sess['waveform'].shape[0]})")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
