import json
import os
import pickle
from collections import OrderedDict
from pathlib import Path


def parse_date_from_path(path: str):
    """Extract the 4th folder segment from a server path like /server/Subjects/subject/date."""
    parts = [part for part in str(path).replace('\\', '/').split('/') if part]
    if len(parts) >= 4:
        return parts[3]
    return None


def infer_subject_probe_loc(umparam_path: str, data_root: str):
    umparam_path_obj = os.path.normpath(umparam_path)
    data_root_obj = os.path.normpath(data_root)
    rel_parts = os.path.relpath(umparam_path_obj, data_root_obj).split(os.sep)
    if len(rel_parts) < 4:
        return None

    # Expected layout: Subject/Probe/Loc/DeepUnitMatch/UMparam.pickle
    subject_name = rel_parts[0]
    probe_name = rel_parts[1]
    loc_name = rel_parts[2]
    return subject_name, probe_name, loc_name


def find_umparam_files(data_root: str):
    root = Path(data_root)
    return sorted(str(p) for p in root.glob('**/DeepUnitMatch/UMparam.pickle') if p.is_file())


def build_metadata_index(data_dir: str, output_name: str = 'metadata_index.json'):
    data_root = os.path.abspath(data_dir)
    print(data_root)
    if not os.path.exists(data_root):
        raise FileNotFoundError(f'Data directory does not exist: {data_root}')

    index = OrderedDict()
    seen = {}

    for umparam_path in find_umparam_files(data_root):
        try:
            with open(umparam_path, 'rb') as f:
                umparam = pickle.load(f)
        except Exception as exc:
            print(f'Skipping {umparam_path}: {exc}')
            continue

        parsed = infer_subject_probe_loc(umparam_path, data_root)
        print(parsed)
        if parsed is None:
            continue

        subject_name, probe_name, loc_name = parsed
        key = (subject_name, probe_name, loc_name)

        ks_dirs = umparam.get('KS_dirs', {})

        dates = OrderedDict()
        exp_paths = OrderedDict()
        exp_ids = OrderedDict()
        for exp_idx, exp_path in enumerate(ks_dirs):
            exp_key = str(exp_idx+1) # index 1??
            exp_paths[exp_key] = str(exp_path)
            exp_ids[exp_key] = str(exp_idx) # index 0??
            date_value = parse_date_from_path(str(exp_path))
            dates[exp_key] = date_value if date_value is not None else ''

        entry = OrderedDict([
            ('mouse', subject_name),
            ('probe', probe_name),
            ('loc', loc_name),
            ('dates', dates),
            ('exp_paths', exp_paths),
            ('exp_ids', exp_ids),
        ])

        seen[key] = {'path': umparam_path, 'entry': entry}

    for key, item in seen.items():
        index[f'{key[0]}/{key[1]}/{key[2]}'] = item['entry']

    output_path = os.path.join(data_root, output_name)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(index, f, indent=4)

    print(f'Wrote {len(index)} entries to {output_path}')
    return output_path


if __name__ == '__main__':
    import sys
    data_dir = sys.argv[1] if len(sys.argv) > 1 else r'\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2'
    output_path = build_metadata_index(data_dir)
    print(output_path)
