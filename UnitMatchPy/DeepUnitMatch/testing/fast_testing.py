import os
import pandas as pd
import numpy as np
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from sklearn.neighbors import KernelDensity
from testing.test import remove_conflicts, directional_filter
from utils.helpers import (
    PROJECT_ROOT,
    index_dates_from_loc,
    pick,
    create_um_lookup,
    get_locations_from_sqlite,
    save_final_results,
    save_intermediate_results,
)


# PROJECT_ROOT (defined in utils.helpers) is the directory holding the shared
# data/results and metadata_index.json alongside the sibling repos.
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")


def test_models_optimized(col_names, fixed_n=True, save_names=None):
    """
    Optimized version of test_models with major performance improvements:
    1. Parallel processing using multiprocessing
    2. Batch database operations
    3. Eliminate repeated CSV writes
    4. Cache expensive computations
    5. Early filtering and validation
    """
    results_dir = RESULTS_DIR

    if fixed_n:
        results_dir = os.path.join(results_dir, "N_set_by_UM")
        # Load UM_results once and create a lookup dictionary for fast access
        UM_results = pd.read_csv(os.path.join(results_dir, "MatchProbNew_results.csv"))
        um_lookup = create_um_lookup(UM_results)
    else:
        um_lookup = None

    # Get all valid locations first to avoid repeated file system checks
    valid_locations = get_locations_from_sqlite(db_path=os.path.join(PROJECT_ROOT, "matchtables.db"))
    print(f"Found {len(valid_locations)} valid locations to process")

    if not valid_locations:
        print("No valid locations found!")
        return

    # Process locations in parallel
    max_workers = min(os.cpu_count() - 1, 6)

    results_accumulator = {col_name: [] for col_name in col_names}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_location = {
            executor.submit(
                process_single_location,
                location,
                col_names,
                um_lookup,
                fixed_n=fixed_n,
            ): location
            for location in valid_locations
        }

        # Process completed jobs
        for i, future in enumerate(as_completed(future_to_location)):
            location = future_to_location[future]
            try:
                location_results = future.result()
                if location_results:
                    for col_name, results in location_results.items():
                        if results is not None:
                            results_accumulator[col_name].append(results)

                # Save intermediate results every 20 locations to avoid losing work
                if (i + 1) % 20 == 0:
                    save_intermediate_results(
                        results_accumulator, results_dir, col_names
                    )
                    print(f"Processed {i + 1}/{len(valid_locations)} locations")

            except Exception as e:
                print(f"Error processing {location}: {e}")
                continue

    save_final_results(results_accumulator, results_dir, save_names=save_names)
    print("Processing complete!")


def get_threshold_df(mt: pd.DataFrame, metric: str = "DNNSim"):
    mt = mt.loc[
        (mt["RecSes1"] == mt["RecSes2"]), :
    ]  # Only use within day rows to compute threshold

    # On-diagonal means same neuron. Off-diagonal means different neurons.
    on_diag = mt.loc[(mt["ID1"] == mt["ID2"]), [metric]]
    off_diag = mt.loc[(mt["ID1"] != mt["ID2"]), [metric]]

    if len(on_diag) < 2 or len(off_diag) < 2:
        # can't do KDE in this case, or define a threshold properly.
        return mt[metric].max()

    # sanity check that the categories are being loaded correctly
    assert sum(on_diag.index.isin(off_diag.index)) == 0

    # Kernel density estimation (distributions are more useful than histograms)
    kde_on = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(on_diag.values)
    kde_off = KernelDensity(kernel="gaussian", bandwidth=0.01).fit(off_diag.values)
    x = np.linspace(min(off_diag[metric]), max(on_diag[metric]), 1000).reshape(-1, 1)
    y_on = np.exp(kde_on.score_samples(x))
    y_off = np.exp(kde_off.score_samples(x))

    # Find the threshold where the distributions intersect
    thresh = np.argwhere(np.diff(np.sign(y_off - y_on)))
    if len(thresh) == 0:
        thresh = len(x) - 1
    elif len(thresh) > 1:
        thresh = thresh[-1]
    return x[thresh].item()


def get_matches_1model(mt, metric, fixed_n: int = None):
    """
    Get the matches for one model, given a merged match table (for one session pair).
    If a fixed number of matches is set, we don't do similarity thresholding.
    Always do conflict resolution.

    Spatial filtering is not applied here: the model column already encodes spatial
    information (e.g. the _40Spat / _80Spat / _NoSpat variants).
    """
    # Filter data once before processing to avoid repeated filtering
    across = mt.loc[
        (mt["RecSes1"] != mt["RecSes2"]), [metric, "RecSes1", "RecSes2", "ID1", "ID2"]
    ]

    if not fixed_n:
        if metric == "MatchProb":
            thresh = 0.5
        else:

            thresh = get_threshold_df(mt, metric=metric)

            within = mt.loc[
                (mt["RecSes1"] == mt["RecSes2"]),
                [metric, "ID1", "ID2", "RecSes1", "RecSes2"],
            ]
            # Correct for different median similarities between within- and across-day sets.
            diff = np.median(within[metric]) - np.median(across[metric])
            thresh = thresh - diff

        # Apply thresholds to generate matches
        matches = across.loc[
            across[metric] >= thresh, ["RecSes1", "RecSes2", "ID1", "ID2", metric]
        ]
    else:
        # For fixed_n, sort once and reuse
        if metric == "CentroidDist":
            matches = across.sort_values(by=metric, ascending=True)
        else:
            matches = across.sort_values(by=metric, ascending=False)

    # Only allow a match if it is above threshold when comparing in both directions
    matches = directional_filter(matches)

    if len(matches) != 0:

        # Resolve conflict matches by only keeping the match with highest similarity
        matches, _ = remove_conflicts(matches, metric)

    # Ensure consistent ordering (keep a single direction per session pair)
    matches = matches.loc[matches["RecSes1"] < matches["RecSes2"]]
    matches = matches.sort_values(by=metric, ascending=False)

    if fixed_n is not None:
        matches = matches.head(fixed_n)

    return matches.index.to_list()


def all_results_1model(mt: pd.DataFrame, mouse: str, probe: str, loc: str, model_name: str, fixed_n=None):
    """
    Get all the results for a single model.
    """
    # Early exit for invalid data
    if len(mt) < 20**2:
        return None

    sessions = sorted(mt["RecSes1"].unique())  # Sort once and reuse

    # Per-session recording dates come straight from the committed metadata index
    # (built by build_metadata_index.py) -- no data drive / metadata.json needed.
    metadata_cache = index_dates_from_loc(mouse, probe, loc)

    results_data = []

    if fixed_n is not None:
        session_pairs = fixed_n[["r1", "r2"]].drop_duplicates().values.tolist()
    else:
        session_pairs = [(r1, r2) for r1 in sessions for r2 in sessions if r1 < r2]

    # Process each session pair
    for r1, r2 in session_pairs:
        # Get the subset of the match table
        df = pick(mt, r1, r2)

        # Get matches
        if fixed_n is not None:
            n_value = fixed_n.loc[
                (fixed_n["r1"] == r1) & (fixed_n["r2"] == r2), "N"
            ].values[0]
            match_indices = get_matches_1model(df, metric=model_name, fixed_n=n_value)
        else:
            match_indices = get_matches_1model(df, metric=model_name, fixed_n=None)

        # Skip if no matches found
        if not match_indices:
            continue

        # Calculate metrics
        if "newRPC" in df.columns and "newISI" in df.columns:
            RpcAuc = AUC(df, match_indices, "newRPC")
            IsiAuc = AUC(df, match_indices, "newISI")
        else:
            RpcAuc = AUC(df, match_indices, "refPopCorr")
            IsiAuc = AUC(df, match_indices, "ISICorr")
        N = len(match_indices)

        date1 = metadata_cache.get(r1)
        date2 = metadata_cache.get(r2)
        if date1 is None or date2 is None:
            print(f"Skipping session pair {r1}, {r2} due to missing metadata")
            continue

        # Store results
        results_data.append(
            {
                "mouse": mouse,
                "probe": probe,
                "loc": loc,
                "day1": date1,
                "day2": date2,
                "r1": r1,
                "r2": r2,
                "AUCrpc": RpcAuc,
                "AUCisi": IsiAuc,
                "N": N,
                "delta_days": (date2 - date1).days,
            }
        )

    # Early exit if no results
    if not results_data:
        return None

    # Create result dataframes
    results = pd.DataFrame(results_data)
    results["day1"] = pd.to_datetime(results["day1"])
    results["day2"] = pd.to_datetime(results["day2"])

    # Normalise days
    day0 = min(min(results["day1"]), min(results["day2"]))
    results["day1"] = (results["day1"] - day0).dt.days
    results["day2"] = (results["day2"] - day0).dt.days

    return results


def process_single_location(location_data, col_names, um_lookup, fixed_n=True):
    """Process a single location - designed to be called in parallel"""
    mouse, probe, loc = location_data

    try:
        # Connect to database (each process needs its own connection)
        conn = sqlite3.connect(os.path.join(PROJECT_ROOT, "matchtables.db"))
        cursor = conn.cursor()

        # Check if table exists and get columns
        try:
            cursor.execute(f"PRAGMA table_info('{mouse}_{probe}_{loc}')")
            columns = [column[1] for column in cursor.fetchall()]
            if not columns:
                print(f"Table {mouse}_{probe}_{loc} has no columns.")
                conn.close()
                return None
        except:
            print(f"Table {mouse}_{probe}_{loc} does not exist.")
            conn.close()
            return None

        # Build SQL query based on available columns
        if "newRPC" in columns:
            base_cols = "ID1,ID2,RecSes1,RecSes2,newRPC,newISI"
        else:
            base_cols = "ID1,ID2,RecSes1,RecSes2,refPopCorr,ISICorr"

        # Check which model columns exist
        existing_cols = [col for col in col_names if col in columns]
        if not existing_cols:
            print(f"No requested columns found in {mouse}_{probe}_{loc}.")
            conn.close()
            return None

        query_cols = base_cols + "," + ",".join(existing_cols)

        # Single database query for all needed data
        mt = pd.read_sql_query(f"SELECT {query_cols} FROM {mouse}_{probe}_{loc};", conn)
        conn.close()

        # Get fixed_n data
        if fixed_n:
            lookup_key = (mouse, probe, loc)
            UM_N = um_lookup.get(lookup_key)
        else:
            UM_N = None

        # Process each model
        results = {}
        for col_name in existing_cols:
            model_name = col_name
            result = all_results_1model(
                mt, mouse, probe, loc, model_name=model_name, fixed_n=UM_N
            )
            if result is not None:
                results[col_name] = result

        print(f"Completed {mouse}/{probe}/{loc}")

        return results if results else None

    except Exception as e:
        print(f"Failed to process {mouse}/{probe}/{loc}: {e}")
        return None


def AUC(mt, indices, func_metric):
    """
    Optimized version of AUC calculation function.
    """
    matches_indices = set(indices)
    P = len(matches_indices)
    if P < 1:
        return np.nan

    # Filter only across-session matches
    across = mt.loc[mt["RecSes1"] != mt["RecSes2"]]

    # Get the total number of non-matches
    N_a = len(across) - P

    if N_a <= 0:  # Edge case protection
        return np.nan

    sorted_indices = across.sort_values(by=func_metric, ascending=False).index

    is_match = np.array([idx in matches_indices for idx in sorted_indices], dtype=bool)

    # Calculate cumulative values
    tp_cumsum = np.cumsum(is_match)
    fp_cumsum = np.cumsum(~is_match)

    # Calculate recall and false positive rate arrays
    recall = tp_cumsum / P
    fpr = fp_cumsum / N_a

    auc = np.trapz(recall, fpr)

    return auc


if __name__ == "__main__":
    start = time.time()

    col_names = [
        "NBProb18mice_40Spat",
        "NBProb18mice_80Spat",
        "NBProb18mice_NoSpat",
        "NBProb18mice_10Spat",
        "NBProbuntrained_NoSpat",
    ]

    test_models_optimized(col_names, fixed_n=False)
    test_models_optimized(col_names, fixed_n=True)
    end = time.time()
    print(f"Total time taken: {end - start} seconds")
