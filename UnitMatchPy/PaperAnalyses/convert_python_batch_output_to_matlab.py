import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io


def _load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def _load_match_table_csv(path):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def _load_waveform_info(path):
    if not os.path.exists(path):
        return None
    data = np.load(path, allow_pickle=True)
    out = {}
    for key in data.files:
        out[key] = data[key]
    return out


def convert_python_output_to_matlab(save_dir, output_mat_path=None):
    """Convert the Python batch output into a MATLAB-compatible UnitMatch.mat file.

    The converter writes a new file alongside the Python outputs and leaves the
    existing files intact. It builds the same core variables MATLAB expects:
    UniqueIDConversion, MatchTable, WaveformInfo and UMparam.
    """
    save_dir = Path(save_dir)
    if output_mat_path is None:
        output_mat_path = save_dir / "UnitMatch.mat"
    else:
        output_mat_path = Path(output_mat_path)

    output_mat_path.parent.mkdir(parents=True, exist_ok=True)
    if output_mat_path.exists():
        raise FileExistsError(
            f"{output_mat_path} already exists; refusing to overwrite"
        )

    umparam = _load_pickle(save_dir / "UMparam.pickle")
    clus_info = _load_pickle(save_dir / "ClusInfo.pickle")

    match_prob = np.load(save_dir / "MatchProb.npy")
    match_table = _load_match_table_csv(save_dir / "MatchTable.csv")
    waveform_info = _load_waveform_info(save_dir / "WaveformInfo.npz")

    original_ids = np.asarray(clus_info["original_ids"]).reshape(-1)
    session_id = np.asarray(clus_info["session_id"]).reshape(-1)
    good_units = clus_info.get("good_units", [])
    if len(good_units) == 0:
        good_id = np.ones(len(original_ids), dtype=bool)
    else:
        good_id = np.zeros(len(original_ids), dtype=bool)
        for sess_idx, units in enumerate(good_units):
            if len(units) == 0:
                continue
            units = np.asarray(units).reshape(-1)
            good_id[np.isin(original_ids, units) & (session_id == sess_idx)] = True

    if match_table is None:
        n_units = len(original_ids)
        xx, yy = np.meshgrid(np.arange(n_units), np.arange(n_units))
        pair_id1 = xx.ravel()
        pair_id2 = yy.ravel()
        recses1 = np.repeat(np.arange(n_units) + 1, n_units)
        recses2 = np.tile(np.arange(n_units) + 1, n_units)
        match_table = pd.DataFrame(
            {
                "ID1": pair_id1,
                "ID2": pair_id2,
                "RecSes1": recses1,
                "RecSes2": recses2,
                "UID1": np.nan,
                "UID2": np.nan,
                "MatchProb": match_prob.ravel(),
                "TotalScore": match_prob.ravel(),
                "EucledianDistance": np.nan,
            }
        )

    # Normalize columns to the MATLAB names expected by downstream code.
    match_table = match_table.copy()
    rename_map = {
        "RecSes 1": "RecSes1",
        "RecSes 2": "RecSes2",
        "UM Probabilities": "MatchProb",
        "Matches": "MatchProb",
        "Matches Currated": "MatchProb",
        "UID orig 1": "UID1",
        "UID orig 2": "UID2",
        "UID Lib 1": "UID1",
        "UID Lib 2": "UID2",
        "UID 1": "UID1",
        "UID 2": "UID2",
        "UID Cons 1": "UID1",
        "UID Cons 2": "UID2",
    }
    for old_name, new_name in rename_map.items():
        if old_name in match_table.columns and new_name not in match_table.columns:
            match_table[new_name] = match_table[old_name]

    if "MatchProb" not in match_table.columns:
        match_table["MatchProb"] = np.asarray(match_prob).ravel()
    if "TotalScore" not in match_table.columns:
        match_table["TotalScore"] = np.asarray(match_prob).ravel()
    if "EucledianDistance" not in match_table.columns:
        match_table["EucledianDistance"] = np.nan

    # MATLAB-compatible column order matching UnitMatch.m.
    matlab_cols = [
        "ID1",
        "ID2",
        "RecSes1",
        "RecSes2",
        "UID1",
        "UID2",
        "MatchProb",
        "TotalScore",
        "EucledianDistance",
    ]
    for col in matlab_cols:
        if col not in match_table.columns:
            match_table[col] = np.nan
    match_table = match_table[
        matlab_cols + [c for c in match_table.columns if c not in matlab_cols]
    ]

    # Reconstruct a MATLAB-compatible table as a struct array with table-like
    # Properties.VariableNames metadata so MATLAB code can access it via
    # MatchTable.Properties.VariableNames.
    table_struct = {
        "Properties": {
            "VariableNames": np.array(list(match_table.columns), dtype=object),
        }
    }
    for col in match_table.columns:
        values = np.asarray(match_table[col])
        if values.dtype.kind in {"U", "S", "O"}:
            values = values.astype(object)
        table_struct[col] = values

    # Build UniqueIDConversion with the same field names MATLAB code expects.
    uid_conversion = {
        "UniqueID": np.asarray(original_ids).reshape(1, -1),
        "OriginalClusID": np.asarray(original_ids).reshape(1, -1),
        "recsesAll": np.asarray(session_id + 1).reshape(1, -1),
        "GoodID": np.asarray(good_id, dtype=int).reshape(1, -1),
    }

    clusinfo = {
        "RecSesID": np.asarray(session_id + 1).reshape(1, -1),
        "Good_ID": np.asarray(good_id, dtype=int).reshape(1, -1),
    }

    # Preserve a lightweight WaveformInfo structure.
    wf_info = {}
    if waveform_info is not None:
        for key, value in waveform_info.items():
            wf_info[key] = np.asarray(value)
    if "ProjectedLocation" not in wf_info:
        wf_info["ProjectedLocation"] = np.zeros((3, len(original_ids), 1))
    if "ProjectedWaveform" not in wf_info:
        wf_info["ProjectedWaveform"] = np.zeros((1, 1))
    if "MaxChannel" not in wf_info:
        wf_info["MaxChannel"] = np.zeros(len(original_ids), dtype=int)

    # Save as a MATLAB v7.3 file.
    scipy.io.savemat(
        str(output_mat_path),
        {
            "clusinfo": clusinfo,
            "UniqueIDConversion": uid_conversion,
            "MatchTable": table_struct,
            "WaveformInfo": wf_info,
            "UMparam": umparam,
        },
        format="5",
    )
    return str(output_mat_path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        raise SystemExit(
            "Usage: convert_python_batch_output_to_matlab.py <save_dir> [output_mat_path]"
        )
    save_dir = sys.argv[1]
    output_mat_path = sys.argv[2] if len(sys.argv) > 2 else None
    print(convert_python_output_to_matlab(save_dir, output_mat_path=output_mat_path))
