import os
import pandas as pd
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed
from func_match import AUC_within_session
from more_testing import all_results_1model
from tqdm import tqdm
from testing.test import pick
import json
import time

# Paths relative to the repo root. PROJECT_ROOT is the repo's parent directory,
# which holds the data/results folders alongside this repo.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = os.path.normpath(os.path.join(REPO_ROOT, os.pardir))
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
ALL_DATA_DIR = os.path.join(PROJECT_ROOT, "ALL_DATA")


def test_models_optimized(col_names, filter50=False, fixed_n=True, save_names=None):
    """
    Optimized version of test_models with major performance improvements:
    1. Parallel processing using multiprocessing
    2. Batch database operations
    3. Eliminate repeated CSV writes
    4. Cache expensive computations
    5. Early filtering and validation
    """
    root = r"/Volumes/UnionSine/DUM_DATA"
    results_dir = RESULTS_DIR

    if filter50:
        results_dir = os.path.join(results_dir, "filter50")

    if fixed_n:
        results_dir = os.path.join(results_dir, "N_set_by_UM")

        # Load UM_results once and create a lookup dictionary for fast access
        UM_results = pd.read_csv(os.path.join(results_dir, "MatchProbNew_results.csv"))
        um_lookup = create_um_lookup(UM_results)

    else:
        um_lookup = None

    # Get all valid locations first to avoid repeated file system checks
    valid_locations = get_valid_locations(root)
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
                filter50,
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


def create_um_lookup(UM_results):
    """Create a fast lookup dictionary for UM_results"""
    lookup = {}
    for _, row in UM_results.iterrows():
        key = (row["mouse"], row["probe"], row["loc"])
        if key not in lookup:
            lookup[key] = []
        lookup[key].append(row.to_dict())
    return {k: pd.DataFrame(v) for k, v in lookup.items()}


def get_valid_locations(root, mice=None):
    """Pre-filter valid locations to avoid repeated file system checks"""
    valid_locations = []

    if mice is None:
        mice = os.listdir(root)

    for mouse in mice:
        mouse_path = os.path.join(root, mouse)
        if not os.path.isdir(mouse_path):
            continue

        for probe in os.listdir(mouse_path):
            probe_path = os.path.join(mouse_path, probe)
            if not os.path.isdir(probe_path):
                continue

            for loc in os.listdir(probe_path):
                loc_path = os.path.join(probe_path, loc)
                if not os.path.isdir(loc_path):
                    continue

                mt_path = os.path.join(loc_path, "merged_mt.csv")
                if os.path.exists(mt_path):
                    valid_locations.append((mouse, probe, loc, mt_path))

    return valid_locations


def get_locations_from_sqlite():
    """Get all locations from the SQLite database"""
    conn = sqlite3.connect("matchtables.db")
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    conn.close()

    locations = []
    for table in tables:
        name = table[0]
        parts = name.split("_")
        if len(parts) == 3:
            mouse, probe, loc = parts
            locations.append((mouse, probe, loc))
    return locations


def process_single_location(
    location_data, col_names, um_lookup, filter50, fixed_n=True
):
    """Process a single location - designed to be called in parallel"""
    mouse, probe, loc, mt_path = location_data

    try:
        # Connect to database (each process needs its own connection)
        conn = sqlite3.connect("matchtables.db")
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
            base_cols = "ID1,ID2,RecSes1,RecSes2,newRPC,newISI,EucledianDistance"
        else:
            base_cols = "ID1,ID2,RecSes1,RecSes2,refPopCorr,ISICorr,EucledianDistance"

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
                mt, mt_path, model_name=model_name, fixed_n=UM_N, filter50=filter50
            )
            if result is not None:
                results[col_name] = result

        print(f"Completed {mouse}/{probe}/{loc}")

        return results if results else None

    except Exception as e:
        print(f"Failed to process {mouse}/{probe}/{loc}: {e}")
        return None


def save_intermediate_results(results_accumulator, results_dir, col_names):
    """Save intermediate results to avoid losing work"""
    for col_name in col_names:
        if results_accumulator[col_name]:
            df = pd.concat(results_accumulator[col_name], ignore_index=True)
            temp_path = os.path.join(results_dir, f"{col_name}_results_temp.csv")
            df.to_csv(temp_path, index=False)


def save_final_results(results_accumulator, results_dir, save_names=None):
    """Save final consolidated results"""
    os.makedirs(results_dir, exist_ok=True)

    for i, col_name in enumerate(results_accumulator):
        if save_names and i < len(save_names):
            save_name = save_names[i]
        else:
            save_name = col_name

        temp_path = os.path.join(results_dir, f"{save_name}_results_temp.csv")
        final_path = os.path.join(results_dir, f"{save_name}_results.csv")

        df_new = pd.concat(results_accumulator[col_name], ignore_index=True)

        if os.path.exists(temp_path):
            df_temp = pd.read_csv(temp_path)
            df_final = pd.concat([df_temp, df_new], ignore_index=True)
            os.remove(temp_path)  # Clean up temp file
        else:
            df_final = df_new

        df_final.drop_duplicates(inplace=True)
        df_final = df_final.reset_index(drop=True)

        df_final.to_csv(final_path, index=False)
        print(f"Saved {len(df_final)} results for {save_name}")


def within_session_AUCs():

    root = ALL_DATA_DIR
    locations = get_valid_locations(root)

    columns1 = "RecSes1,RecSes2,ID1,ID2,refPopCorr,ISICorr"
    columns2 = "RecSes1,RecSes2,ID1,ID2,newRPC,newISI"

    for location in tqdm(locations):
        mouse, probe, loc, mt_path = location
        new_path = os.path.join(root, mouse, probe, loc, "AUC_info.json")

        conn = sqlite3.connect("matchtables.db")
        cursor = conn.cursor()
        cursor.execute(f"PRAGMA table_info('{mouse}_{probe}_{loc}')")
        cols = [column[1] for column in cursor.fetchall()]
        if "newISI" not in cols or "newRPC" not in cols:
            mt = pd.read_sql_query(
                f"SELECT {columns1} from {mouse}_{probe}_{loc}", conn
            )
            new = False
        else:
            mt = pd.read_sql_query(
                f"SELECT {columns2} from {mouse}_{probe}_{loc}", conn
            )
            new = True
        sessions = mt["RecSes1"].unique()
        data = {}

        for ses in sessions:
            df = pick(mt, ses, ses, True)
            matches = df.loc[df["ID1"] == df["ID2"]].index
            rpcauc = AUC_within_session(df, matches, "newRPC" if new else "refPopCorr")
            isiauc = AUC_within_session(df, matches, "newISI" if new else "ISICorr")
            N = len(df["ID1"].unique())
            data[int(ses)] = {"N": N, "RPC_AUC": rpcauc, "ISI_AUC": isiauc}

        with open(new_path, "w") as f:
            json.dump(data, f, indent=4)


if __name__ == "__main__":
    start = time.time()

    # col_names = ["NBProb2", "NBProb3", "NBProbuntrained2", "NBProbuntrained3", "NBProbunfinetuned2", "NBProbunfinetuned3", "NBProbuntrainedAE2", "NBProbuntrainedAE3"]
    # col_names = ['MatchProbNew', 'NBProbuntrained', 'NBProbuntrainedAE', 'NBProbunfinetuned']
    col_names = [
        "NBProb18mice_40Spat",
        "NBProb18mice_80Spat",
        "NBProb18mice_NoSpat",
        "NBProb18mice_10Spat",
        "NBProbuntrained_NoSpat",
    ]

    test_models_optimized(col_names, filter50=False, fixed_n=False)
    test_models_optimized(col_names, filter50=False, fixed_n=True)
    end = time.time()
    print(f"Total time taken: {end - start} seconds")

    # within_session_AUCs()

    # mouse, probe, loc = "AL031", "19011116684", '1'
    # mt_path = r"C:\Users\suyash\UCL\DeepUnitMatch\ALL_DATA\AL031\19011116684\1\merged_mt.csv"
    # process_single_location((mouse, probe, loc, mt_path), col_names=['NBProb18mice'], um_lookup=None, filter50=False, fixed_n=False)
