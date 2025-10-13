import os, sys
sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.path.dirname(os.getcwd()), 'UnitMatchPy'))
import bayes_functions as bf
import utils as util
import overlord as ov
import numpy as np
import matplotlib.pyplot as plt
import save_utils as su
import GUI as gui
import assign_unique_id as aid
import default_params as default_params
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import sqlite3


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
        
        # if util.retry_db_operation(_check_if_processed):
        #     print(f"Skipping {mouse}, {probe}, {loc} as it is already processed.")
        #     conn.close()
        #     return None

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

        prior_match = 1 - (param['n_expected_matches'] / param['n_units']**2 ) # freedom of choose in prior prob
        priors = np.array((prior_match, 1-prior_match))

        # Construct distributions (kernels) for Naive Bayes Classifier
        labels = candidate_pairs.astype(int)
        cond = np.unique(labels)
        parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one = 1)

        # Get probability of each pair of being a match
        probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)
        
        # Add column to SQL with retry logic
        util.add_col_to_sql(conn, cursor, mouse, probe, loc, "MatchProbNew", probability[:, 1], overwrite=True)
        
        # Close connection properly
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

def check_ids(mouse, probe, loc):

    directory = os.path.join(r"C:\Users\suyash\UCL\DeepUnitMatch\ALL_DATA", mouse, probe, loc)
    conn = sqlite3.connect(r"C:\Users\suyash\UCL\DeepUnitMatch\DeepMatch\matchtables.db")
    mt = pd.read_sql_query(f"SELECT RecSes1,RecSes2,ID1,ID2 FROM {mouse}_{probe}_{loc}", conn)


if __name__ == "__main__":
    param = default_params.get_default_param()
    root = r"C:\Users\suyash\UCL\DeepUnitMatch\ALL_DATA"

    # paths = pd.DataFrame(util.read_datapaths(os.listdir(root)))
    # valid_locs = util.get_valid_locations(root)
    # # max_workers = min(os.cpu_count() - 1, 12)
    # max_workers = 2

    # with ProcessPoolExecutor(max_workers=max_workers) as executor:
    #     # Submit all jobs
    #     future_to_location = {
    #         executor.submit(process_single_location, location, paths, param): location
    #         for location in valid_locs
    #     }

    mouse, probe, loc = "JF084", "Probe1", "IMRO_10"
    paths = pd.DataFrame(util.read_datapaths([mouse]))

    prob = process_single_location((mouse, probe, loc, os.path.join(root, mouse, probe, loc, "merged_mt.csv")), paths, param)






