# Single source of truth for every path/setting the paper pipeline shares.
#
# All batch, training, comparison and figure scripts read their roots from
# here instead of hard-coding them, so pointing the whole pipeline at a new set
# of folders (e.g. for a full rerun) is a one-line change of RUN_NAME below.
#
# Pipeline order and which root each stage reads/writes:
#   step 1  run_deepunitmatch_batch.py          RAW_KS_BASE      -> UNMERGED_OUTPUT
#           generate_metadata_index.py          UNMERGED_OUTPUT  -> METADATA_INDEX_PATH
#   step 2  generate_merged_dataset.py          UNMERGED_OUTPUT  -> MERGED_DATA
#   step 3  model training                      MERGED_DATA      -> TRAINING_CACHE (+ local checkpoints)
#   step 4  run_*_onMerged*.py, EMD, DANT       MERGED_DATA      -> ANALYSIS_OUTPUT
#   step 6  sql.py / fast_testing.py / figures  ANALYSIS_OUTPUT  -> DATABASE_PATH, RESULTS_DIR
#   any     check_pipeline_completeness.py      all of the above -> REPORTS_DIR
#
# RAW_KS_BASE is input only and is never written to.

import os

# Network share holding all inputs and outputs.
SHARE_ROOT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE"

# Raw Kilosort data + the original UnitMatch.mat per recording location (input only).
RAW_KS_BASE = os.path.join(SHARE_ROOT, "FullAnimal_KSChanMap")

# Name of this run of the pipeline. Every output root below is derived from
# it, so a rerun writes into fresh folders and can never mix with (or be
# skipped because of) outputs from an earlier run.
# Previous run (summer 2026): DeepUM_NatMeth2026V2 (unmerged),
# DeepUM_NatMeth2026V2_merged/merged_data_v2 (merged data),
# DeepUM_NatMeth2026_V3_OnMergedData (analysis output).
RUN_NAME = "DeepUM_Oct2026"

# step 1: DeepUnitMatch/UMPy on the non-merged data (also the input of step 2).
UNMERGED_OUTPUT = os.path.join(SHARE_ROOT, f"{RUN_NAME}_unmerged")
# Recording dates per session, built from UNMERGED_OUTPUT by generate_metadata_index.py.
METADATA_INDEX_PATH = os.path.join(UNMERGED_OUTPUT, "metadata_index.json")

# step 2: merged dataset (KS-style tree, one folder per session).
MERGED_ROOT = os.path.join(SHARE_ROOT, f"{RUN_NAME}_merged")
MERGED_DATA = os.path.join(MERGED_ROOT, "merged_data")

# step 3 (train_paper_models.py): preprocessed training snippets, per-job
# bookkeeping (status, locks, mouse subsets), and the final checkpoint of every
# trained model, copied here from the training machine's local ModelExp so
# any machine can evaluate it.
TRAINING_CACHE = os.path.join(MERGED_ROOT, "training_snippets")
TRAINING_STATE_ROOT = os.path.join(MERGED_ROOT, "training_state")
MODELS_ROOT = os.path.join(MERGED_ROOT, "models")

# step 4 / 6: every matching output on the merged data, plus everything derived from it.
ANALYSIS_OUTPUT = os.path.join(SHARE_ROOT, f"{RUN_NAME}_OnMergedData")
RESULTS_DIR = os.path.join(ANALYSIS_OUTPUT, "results")
FIGURES_DIR = os.path.join(ANALYSIS_OUTPUT, "figures")
REPORTS_DIR = os.path.join(ANALYSIS_OUTPUT, "pipeline_reports")
# The DeepUnitMatch model of every step-4 run on the merged data
# (DeepUnitMatch/, DUM_legacy/, DUM_nodrift/, DUM_maxdist=*/ ...): the main
# paper model trained by train_paper_models.py (--jobs paper); the same model
# is shipped as the package default, DeepUnitMatch/utils/model.
PAPER_MODEL_PATH = os.path.join(MODELS_ROOT, "DUM_paper", "model.pt")
# Step 1 uses the model from before this rerun (the paper model is trained on
# step 1's output), kept in the package as utils/model_PreJuly2026.
STEP1_MODEL_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "DeepUnitMatch", "utils", "model_PreJuly2026")

# step 6: sql.py puts the match tables of a group of models (output folder
# names) into one SQLite database per group, matchtables_<group>.db, and
# fast_testing.py evaluates a group from its database (both: --group).
# One database for every model would be far too large (the previous run's
# per-group databases were 40-130 GB each). Every group starts with
# DeepUnitMatch: fast_testing always reads its natural-image and
# reference-population columns, and its free-N match count sets N for the
# fixed-N comparison.
MODEL_GROUPS = {
    # main comparison (run_deepunitmatch_batch_onMerged.py, EMD, DANT)
    "algo": ["DeepUnitMatch", "UMPy", "DUM_legacy", "EMD",
             "DANT", "DANT_no_functional", "DANT_fixed", "DANT_no_functional_fixed"],
    # paper models (train_paper_models.py --jobs paper; DUM_paper = DeepUnitMatch)
    "models": ["DeepUnitMatch", *[f"DUM_W_ij={w}" for w in (1, 5, 15, 20)],
               *[f"DUM_n_output={n}" for n in (8, 32, 128)],
               "DUM_untrained", "DUM_AE_only", "DUM_finetuned_only"],
    # pipeline conditions
    "nodrift": ["DeepUnitMatch", "UMPy", "DUM_nodrift", "UMPy_nodrift"],
    "spatial_only": ["DeepUnitMatch", "UMPy", "UMPy_spatialonly"],
    # pipeline variants: which parts of the unified pipeline matter for DUM
    # (run_deepunitmatch_batch_onMerged_pipelines.py)
    "pipelines_DUM": ["DeepUnitMatch","DUM_legacy",
                  "DUM_pairCandidates", "DUM_pairBayes", "DUM_pairCandidates_pairBayes", "DUM_weights1to1"],
    "pipelines_DUM2": ["DeepUnitMatch","DUM_pairCandidates_pairBayes_pairPrior", "DUM_pairCandidates_pairBayes_pairPrior_bothDirections"],
    "pipelines_UMPy": ["DeepUnitMatch", "UMPy",
                  "UMPy_pairCandidates", "UMPy_pairBayes", "UMPy_pairCandidates_pairBayes", "UMPy_weights1to1"],
    "pipelines_UMPy2": ["DeepUnitMatch","UMPy_pairCandidates_pairBayes_pairPrior", "UMPy_pairCandidates_pairBayes_pairPrior_bothDirections"],
    "maxdist": ["DeepUnitMatch", "UMPy",
                *[f"{m}_maxdist={d}" for d in ("20", "50", "100", "inf") for m in ("DUM", "UMPy")]],
    # mouse-count cross-validation (the 18-mouse point is DeepUnitMatch)
    "xval_1": ["DeepUnitMatch", *[f"xval_m1_{r}" for r in (1, 2, 3)]],
    "xval_3-6": ["DeepUnitMatch", *[f"xval_m{k}_{r}" for k in (3, 6) for r in (1, 2, 3)]],
    "xval_12": ["DeepUnitMatch", *[f"xval_m12_{r}" for r in (1, 2, 3)]],
}


def database_path(group):
    """SQLite database of one model group (see MODEL_GROUPS)."""
    return os.path.join(ANALYSIS_OUTPUT, f"matchtables_{group}.db")


DATABASE_PATH = database_path("algo")  # the main comparison

# Batch scripts skip work whose output sentinel already exists. With fresh
# output roots per run this stays None ("skip if present"); set a
# datetime.datetime only to force recomputing outputs older than that moment
# within the current run (see batch_lock.sentinel_is_fresh).
REDO_FROM_DATE = None

# Per-location done/failed/skipped events of every stage (pipeline_log.py),
# one JSON-lines file per stage and machine; read by check_pipeline_completeness.py.
LOG_DIR = os.path.join(SHARE_ROOT, f"{RUN_NAME}_pipeline_logs")

# Every root the batch scripts may leave .processing*.lock files in
# (used by clear_stale_locks.py).
LOCK_ROOTS = [UNMERGED_OUTPUT, MERGED_ROOT, ANALYSIS_OUTPUT]
