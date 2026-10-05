# Completeness check for the paper pipeline: which mice / recording locations
# made it through each stage, and which were silently dropped.
#
# Background: a whole recording location (CB018 location 1: 6 sessions, 303
# good units) never got step-1 output -- run_deepunitmatch_batch.py printed an
# error to the console and moved on, nothing was logged -- so CB018 only
# entered the analysis through a single-session location and was then excluded
# for having no across-session matches. This script makes such losses visible.
#
# Stages checked, per recording location ("group" = mouse/probe/location):
#   raw        UnitMatch.mat under the raw KS tree (run_deepunitmatch_batch.BASE_INPUT):
#              sessions, and sessions/units that have good units
#   step1      DeepUnitMatch/ and UMPy/ output on the non-merged data
#              (run_deepunitmatch_batch.BASE_OUTPUT), via the same
#              MatchingOverview.png sentinel the batch script uses
#   merged     the merged tree written by generate_merged_dataset.py: every
#              session folder present, merge_complete.flag written, good and
#              merged-away (MERGED) units
#   onmerged   DeepUnitMatch/ and UMPy/ output on the merged data
#              (run_deepunitmatch_batch_onMerged.BASE_OUTPUT)
#
# A location with fewer than 2 sessions containing good units can never yield
# across-session matches; it is reported as "single session", not as an error.
#
# Output (in --out, default <merged-analysis output root>/pipeline_reports):
#   completeness_<timestamp>.csv   one row per location, all stages
#   completeness_<timestamp>.txt   per-mouse summary + every problem found
# Exit code 1 if any location is missing from a stage it should have reached.
#
# Read-only with respect to all data/output trees; only writes the report.
#
# Example:
#   python check_pipeline_completeness.py                 # all stages, all mice
#   python check_pipeline_completeness.py --stages raw step1 merged
#   python check_pipeline_completeness.py --mice CB018 --out C:\temp\reports

import argparse
import datetime
import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))

import run_deepunitmatch_batch as step1
import run_deepunitmatch_batch_onMerged as onmerged
import generate_merged_dataset as gm
import pipeline_log as plog

STAGES = ("raw", "step1", "merged", "onmerged")
GOOD_LABELS = {"GOOD", "NON-SOMA GOOD"}
SENTINEL = "MatchingOverview.png"


def find_raw_groups(mice=None):
    """{group key 'mouse/probe/loc': UnitMatch.mat path} under the raw KS tree."""
    groups = {}
    for root, dirs, files in os.walk(step1.BASE_INPUT):
        rel = os.path.relpath(root, step1.BASE_INPUT)
        if mice and rel != "." and rel.split(os.sep)[0] not in mice:
            dirs[:] = []
            continue
        if os.path.basename(root) == "UnitMatch" and "UnitMatch.mat" in files:
            key = os.path.dirname(rel).replace(os.sep, "/")
            groups[key] = os.path.join(root, "UnitMatch.mat")
            dirs[:] = []
    return groups


def raw_info(mat_path):
    ks_dirs, _, recses_all, good_id = step1.load_unitmatchemat(mat_path)
    good_per_session = [int(((recses_all == i + 1) & good_id.astype(bool)).sum()) for i in range(len(ks_dirs))]
    return {
        "raw_sessions": len(ks_dirs),
        "raw_sessions_with_good": sum(n > 0 for n in good_per_session),
        "raw_good_units": sum(good_per_session),
    }


def output_done(root, key, model):
    return os.path.isfile(os.path.join(root, *key.split("/"), model, SENTINEL))


def merged_info(key, n_ks_dirs):
    """Per-session state of the merged tree for one location (layout as written by generate_merged_dataset)."""
    group_dir = os.path.join(onmerged.BASE_INPUT, *key.split("/"), "DeepUnitMatch")
    info = {"merged_sessions": 0, "merged_flags": 0, "merged_good_units": 0, "merged_absorbed": 0}
    if not os.path.isdir(group_dir):
        return info
    for idx in range(n_ks_dirs if n_ks_dirs is not None else 0):
        sess_dir = os.path.join(group_dir, str(idx))
        if not os.path.isdir(sess_dir):
            continue
        info["merged_sessions"] += 1
        info["merged_flags"] += os.path.isfile(os.path.join(sess_dir, gm.MERGE_COMPLETE_MARKER))
        tsv = os.path.join(sess_dir, "cluster_bc_unitType.tsv")
        if os.path.isfile(tsv):
            labels = pd.read_csv(tsv, sep="\t")["bc_unitType"].astype(str).str.upper()
            info["merged_good_units"] += int(labels.isin(GOOD_LABELS).sum())
            info["merged_absorbed"] += int((labels == gm.MERGED_UNIT_LABEL).sum())
    return info


def check(stages, mice=None):
    raw_groups = find_raw_groups(mice)
    step1_groups = set()
    for root, dirs, files in os.walk(step1.BASE_OUTPUT):
        if os.path.basename(root) in ("DeepUnitMatch", "UMPy"):
            key = os.path.relpath(os.path.dirname(root), step1.BASE_OUTPUT).replace(os.sep, "/")
            if not mice or key.split("/")[0] in mice:
                step1_groups.add(key)
            dirs[:] = []
    keys = sorted(set(raw_groups) | step1_groups)
    logged = plog.latest_events()  # {(stage, group, condition): latest event}

    rows = []
    for i, key in enumerate(keys):
        print(f"[{i + 1}/{len(keys)}] {key}", flush=True)
        row = {"group": key, "mouse": key.split("/")[0], "problems": []}
        n_ks_dirs = None
        if key not in raw_groups:
            row["problems"].append("output exists but no raw UnitMatch.mat")
        elif "raw" in stages:
            try:
                row.update(raw_info(raw_groups[key]))
                n_ks_dirs = row["raw_sessions"]
            except Exception as e:
                row["problems"].append(f"raw UnitMatch.mat unreadable: {e}")
        single = row.get("raw_sessions_with_good", 2) < 2
        row["single_session"] = single

        if "step1" in stages:
            for model in ("DeepUnitMatch", "UMPy"):
                done = output_done(step1.BASE_OUTPUT, key, model)
                row[f"step1_{model}"] = done
                if not done:
                    row["problems"].append(f"step1 {model} output missing")

        if "merged" in stages:
            row.update(merged_info(key, n_ks_dirs))
            if n_ks_dirs is not None and row["merged_flags"] < n_ks_dirs:
                row["problems"].append(
                    f"merged tree incomplete ({row['merged_flags']}/{n_ks_dirs} sessions with merge_complete.flag)"
                )

        if "onmerged" in stages:
            for model in ("DeepUnitMatch", "UMPy"):
                done = output_done(onmerged.BASE_OUTPUT, key, model)
                row[f"onmerged_{model}"] = done
                if not done:
                    row["problems"].append(f"merged-data {model} output missing")

        # why: latest logged failure/skip of this location in any stage (pipeline_log.py)
        reasons = [
            f"[{e['stage']}/{e['condition']} {e['status']} {e['time']} on {e['machine']}] {e['message']}"
            for (stage, group, _), e in logged.items()
            if group == key and e["status"] in ("failed", "skipped")
        ]
        row["logged_reason"] = " | ".join(reasons)
        rows.append(row)
    return pd.DataFrame(rows)


def summarise(df):
    lines = []
    per_mouse = df.groupby("mouse").agg(
        locations=("group", "size"),
        single_session=("single_session", "sum"),
        with_problems=("problems", lambda p: sum(len(x) > 0 for x in p)),
    )
    trackable = df[~df["single_session"] & (df["problems"].str.len() == 0)].groupby("mouse").size()
    per_mouse["complete_trackable"] = trackable.reindex(per_mouse.index).fillna(0).astype(int)
    lines.append("Per mouse (complete_trackable = locations with >=2 sessions and no problems):")
    lines.append(per_mouse.to_string())
    lost = per_mouse.index[per_mouse["complete_trackable"] == 0].tolist()
    lines.append(f"\n{len(per_mouse)} mice found; {len(per_mouse) - len(lost)} have at least one complete, trackable location.")
    if lost:
        lines.append(f"MICE WITHOUT ANY COMPLETE TRACKABLE LOCATION: {', '.join(lost)}")
    problems = df[df["problems"].str.len() > 0]
    lines.append(f"\n{len(problems)} location(s) with problems:")
    for _, r in problems.iterrows():
        lines.append(f"  {r['group']}: " + "; ".join(r["problems"]))
        if r.get("logged_reason"):
            lines.append(f"      logged: {r['logged_reason']}")
        else:
            lines.append("      logged: no failure logged (never attempted, or run before logging existed)")
    singles = df[df["single_session"]]["group"].tolist()
    if singles:
        lines.append(f"\nSingle-session locations (cannot be tracked, not an error): {', '.join(singles)}")
    return "\n".join(lines), len(problems)


def main():
    parser = argparse.ArgumentParser(description="Check which mice/locations made it through each pipeline stage.")
    parser.add_argument("--stages", nargs="+", choices=STAGES, default=list(STAGES))
    parser.add_argument("--mice", nargs="+", default=None, help="Only check these mice")
    parser.add_argument("--out", default=os.path.join(onmerged.BASE_OUTPUT, "pipeline_reports"))
    args = parser.parse_args()

    df = check(args.stages, mice=args.mice)
    text, n_problems = summarise(df)
    print("\n" + text)

    os.makedirs(args.out, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.assign(problems=df["problems"].str.join("; ")).to_csv(
        os.path.join(args.out, f"completeness_{stamp}.csv"), index=False
    )
    with open(os.path.join(args.out, f"completeness_{stamp}.txt"), "w") as f:
        f.write(f"stages: {args.stages}\nmice: {args.mice or 'all'}\n\n{text}\n")
    print(f"\nReport written to {args.out}")
    sys.exit(1 if n_problems else 0)


if __name__ == "__main__":
    main()
