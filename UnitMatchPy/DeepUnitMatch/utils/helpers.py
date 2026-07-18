import pandas as pd
import numpy as np
import os
import h5py
import json
import sqlite3
import datetime


# PROJECT_ROOT is the directory that holds both sibling repos (DeepMatch + DeepUnitMatch)
# alongside the shared data/results and the committed metadata_index.json.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = os.path.normpath(os.path.join(REPO_ROOT, os.pardir, os.pardir, os.pardir))
METADATA_INDEX_PATH = os.path.join(PROJECT_ROOT, "metadata_index.json")

_metadata_index_cache = None


def create_dataframe(good_units, prob_matrix, session_list=None):
    # Create dataframe for all session pairs

    # Get number of sessions and units per session
    if session_list is None:
        n_sessions = len(good_units)
        session_list = list(range(n_sessions))
    session_switch = np.cumsum(
        [np.asarray(units).squeeze().shape[0] for units in good_units]
    )
    session_switch = np.insert(session_switch, 0, 0)

    d = {
        "RecSes1": [],
        "RecSes2": [],
        "ID1": [],
        "ID2": [],
        "Prob": [],
    }

    # Generate all pairs of sessions (including within-session pairs)
    for i, ses1 in enumerate(session_list):
        for j, ses2 in enumerate(session_list):
            units_ses1 = np.asarray(good_units[i]).squeeze().astype(int)
            units_ses2 = np.asarray(good_units[j]).squeeze().astype(int)
            n_units_ses1 = units_ses1.shape[0]
            n_units_ses2 = units_ses2.shape[0]

            # Extract probability matrix block for this session pair
            prob_block = prob_matrix[
                session_switch[i] : session_switch[i + 1],
                session_switch[j] : session_switch[j + 1],
            ]

            # Add all unit pairs from this session pair
            d["RecSes1"].extend([ses1] * (n_units_ses1 * n_units_ses2))
            d["RecSes2"].extend([ses2] * (n_units_ses1 * n_units_ses2))
            d["ID1"].extend(np.repeat(units_ses1, n_units_ses2).tolist())
            d["ID2"].extend(np.tile(units_ses2, n_units_ses1).tolist())
            d["Prob"].extend(prob_block.ravel().tolist())

    df = pd.DataFrame(d)
    return df


def get_unit_id(filepath: str):
    fp = os.path.basename(filepath)
    if fp[:4] == "Unit" and fp[-14:] == "_RawSpikes.npy":
        fp = fp.replace("Unit", "")
        id = fp.replace("_RawSpikes.npy", "")
        if "+" in id:
            id = id[: id.find("+")]
        if "#" in id:
            id = id.replace("#", "")
        try:
            return int(id)
        except:
            print(id)
            raise ValueError(f"Invalid filepath format for this waveform: {filepath}")
    else:
        raise ValueError(f"Invalid filepath format for this waveform: {filepath}")


def read_pos(path, skip_removed_units=True):
    files = os.listdir(path)
    x = []
    y = []
    filenames = []

    for file in files:
        if skip_removed_units:
            if "+" in file or "#" in file:
                continue
        fp = os.path.join(path, file)
        with h5py.File(fp, "r") as f:
            # waveform = f['waveform'][()]
            MaxSitepos = f["MaxSitepos"][()]
        x.append(MaxSitepos[0])
        y.append(MaxSitepos[1])
        filenames.append(get_unit_id(file))
    output = {"ID": filenames, "x": x, "y": y}
    return pd.DataFrame(output).sort_values("ID").reset_index(drop=True)

def _load_metadata_index():
    """Load and cache the committed metadata_index.json (built once per process)."""
    global _metadata_index_cache
    if _metadata_index_cache is None:
        if not os.path.isfile(METADATA_INDEX_PATH):
            raise FileNotFoundError(
                f"Metadata index not found at {METADATA_INDEX_PATH}; "
                f"run build_metadata_index.py to create it."
            )
        with open(METADATA_INDEX_PATH) as f:
            _metadata_index_cache = json.load(f)
    return _metadata_index_cache

def loc_key_from_mt_path(mt_path: str):
    """Machine-independent 'mouse/probe/loc' key for a match-table path.

    mt_path is always <data_root>/<mouse>/<probe>/<loc>/<mt_name>.csv, so the key
    is derived purely from the trailing path components -- independent of the data
    root prefix and of the match-table filename. Used by both the metadata index
    builder (build_metadata_index.py) and expids_from_index() so the keys agree.
    """
    loc_dir = os.path.dirname(mt_path)
    loc = os.path.basename(loc_dir)
    probe = os.path.basename(os.path.dirname(loc_dir))
    mouse = os.path.basename(os.path.dirname(os.path.dirname(loc_dir)))
    return f"{mouse}/{probe}/{loc}"

def _index_entry(key: str):
    """Look up a 'mouse/probe/loc' key in the committed metadata index."""
    index = _load_metadata_index()
    if key not in index:
        raise KeyError(
            f"'{key}' is not in the metadata index ({METADATA_INDEX_PATH}); "
            f"re-run build_metadata_index.py if this location is new."
        )
    return index[key]


def expids_from_index(mt_path: str):
    """Data-drive-free replacement for expids_from_metadata.

    Returns (exp_ids, metadata) with the exact same shape as expids_from_metadata,
    but reads everything from the committed metadata_index.json instead of scanning
    the loc directory's per-experiment metadata.json files. Only the mt_path (its
    trailing mouse/probe/loc components) is needed -- the data drive does not have
    to be mounted.

    exp_ids maps each RecSes (int) to its experiment folder name; metadata is the
    loc-level dict with "mouse", "probe" and "loc".
    """
    entry = _index_entry(loc_key_from_mt_path(mt_path))
    exp_ids = {int(recses): name for recses, name in entry["exp_ids"].items()}
    metadata = {k: entry[k] for k in ("mouse", "probe", "loc")}
    return exp_ids, metadata


def index_dates_from_loc(mouse, probe, loc):
    """Return {RecSes (int): datetime.date | None} for a location from the index.

    Reads the per-session recording dates baked into metadata_index.json by
    build_metadata_index.py, so no per-experiment metadata.json (and no mounted
    data drive) is needed at test time.
    """
    entry = _index_entry(f"{mouse}/{probe}/{loc}")
    dates = {}
    for recses, date_str in entry.get("dates", {}).items():
        dates[int(recses)] = (
            datetime.date.fromisoformat(date_str) if date_str else None
        )
    return dates

def pick(mt, r1, r2, tight: bool = False):
    if tight:
        return mt.loc[(mt["RecSes1"] == r1) & (mt["RecSes2"] == r2)]
    else:
        return mt.loc[(mt["RecSes1"].isin([r1, r2])) & (mt["RecSes2"].isin([r1, r2]))]


def create_um_lookup(UM_results):
    """Create a fast lookup dictionary for UM_results"""
    lookup = {}
    for _, row in UM_results.iterrows():
        key = (row["mouse"], row["probe"], row["loc"])
        if key not in lookup:
            lookup[key] = []
        lookup[key].append(row.to_dict())
    return {k: pd.DataFrame(v) for k, v in lookup.items()}


def get_locations_from_sqlite(db_path="matchtables.db"):
    """Get all locations from the SQLite database"""
    conn = sqlite3.connect(db_path)
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
