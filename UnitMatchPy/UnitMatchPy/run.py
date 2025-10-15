# This is an example usage script. It will need to be modified to run on your data.
# The script processes data from lots of mice concurrently so it is fast.
# The end result of the script is to write MatchProb values to a SQL database with a match table for each mouse/probe/location combination.
# It also assumes that the waveform snippets are stored in a specific directory structure: root/mouse/probe/location/experiment_id/processed_waveforms

import os, sys
sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.path.dirname(os.getcwd()), 'UnitMatchPy'))
import bayes_functions as bf
import utils as util
import overlord as ov
import numpy as np
import default_params as default_params
from concurrent.futures import ProcessPoolExecutor
import pandas as pd


def process_single_location(location_data, path_lookup, param):

    try:
        mouse, probe, loc, mt_path = location_data

        print(f"Processing {mouse}, {probe}, {loc}")

        # Check if already processed (with retry logic)
        table_name = f'{mouse}_{probe}_{loc}'

        # Create robust database connection with retry logic
        conn, cursor = util.create_robust_db_connection(r"C:\Users\suyash\UCL\DeepUnitMatch\DeepMatch\matchtables.db")
        
        def _check_if_processed():
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [column[1] for column in cursor.fetchall()]
            if "MatchProbNew" in columns:
                # additionally check if the column contains NULLs
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE MatchProbNew IS NULL")
                count_null = cursor.fetchone()[0]
                if count_null == 0:
                    return True
            return False
        
        if util.retry_db_operation(_check_if_processed):
            print(f"Skipping {mouse}, {probe}, {loc} as it is already processed.")
            conn.close()
            return None

        KS_dirs = path_lookup.loc[(path_lookup['mouse'] == mouse) & (path_lookup['probe'] == probe) & (path_lookup['loc'] == loc), 'recordings'].to_list()[0]
        param['KS_dirs'] = KS_dirs
        wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(KS_dirs)
        param = util.get_probe_geometry(channel_pos[0], param)

        waveform, session_id, session_switch, within_session, good_units, param = util.load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only = False)
        waveform, session_switch, within_session, session_id, good_units, param = util.filter_good_units_and_merge(mt_path, mouse, KS_dirs, waveform, session_switch, param)
        # waveform, session_switch, within_session, session_id, good_units, param = util.filter_good_units(mouse, probe, loc, conn, waveform, session_switch, param)

        mt = pd.read_sql_query(f"SELECT RecSes1,RecSes2,ID1,ID2 FROM {table_name}", conn)
        ids = {}
        for session in mt['RecSes1'].unique():
            df = mt.loc[mt['RecSes1'] == session]
            ids[session] = df['ID1'].unique()

        clus_info = {'session_switch' : session_switch, 'session_id' : session_id}

        # STEP 1
        # Extract parameters from waveform
        extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)

        # STEP 2, 3, 4
        # Extract metric scores
        total_score, candidate_pairs, scores_to_include, predictors  = ov.extract_metric_scores(extracted_wave_properties, session_switch, within_session, param, niter  = 2)

        prior_match = 1 - (param['n_expected_matches'] / param['n_units']**2 )
        priors = np.array((prior_match, 1-prior_match))

        labels = candidate_pairs.astype(int)
        cond = np.unique(labels)
        parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one = 1)
        probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)

        util.add_col_to_sql(conn, cursor, mouse, probe, loc, "MatchProbNew", probability[:, 1], overwrite=True)
        conn.close()

        print(f"Finished processing {mouse}, {probe}, {loc}")
        return probability[:, 1]
    
    except Exception as e:
        print(f"Error processing {mt_path}: {e}")
        try:
            if 'conn' in locals():
                conn.close()
        except:
            pass
        return None


if __name__ == "__main__":
    param = default_params.get_default_param()
    root = r"path/to/data/root"  # Root directory containing all data folders - data should be organised under this as mouse/probe/location

    paths = pd.DataFrame(util.read_datapaths(os.listdir(root), r'/server/data/base'))
    valid_locs = util.get_valid_locations(root)
    # paths is a pandas dataframe with columns mouse, probe, loc, recordings, where recordings is a list of KS_dirs for that mouse/probe/location
    # valid_locs is a list of tuples (mouse, probe, loc, mt_path) for each location with a match table
    
    max_workers = min(os.cpu_count() - 1, 12)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_location = {
            executor.submit(process_single_location, location, paths, param): location
            for location in valid_locs
        }

