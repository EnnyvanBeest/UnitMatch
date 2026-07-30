import sqlite3
import pandas as pd
import os
import numpy as np

# Paths relative to the repo root. PROJECT_ROOT is the repo's parent directory,
# which holds the data/results folders alongside this repo.
# PROJECT_ROOT = r'\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2'
PROJECT_ROOT = r'\\znas\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026_V3_OnMergedData'


def pandas_to_sqlite_type(dtype):
    """Maps pandas dtype to a suitable SQLite data type string."""
    if pd.api.types.is_integer_dtype(dtype):
        return "INTEGER"
    elif pd.api.types.is_float_dtype(dtype):
        return "REAL"
    elif pd.api.types.is_bool_dtype(dtype):
        return "INTEGER"  # Store bools as 0 or 1
    elif pd.api.types.is_datetime64_any_dtype(dtype):
        return "TEXT"  # Store datetimes as ISO format strings
    # Add other specific types if needed (e.g., timedelta)
    else:
        # Default to TEXT for object, string, category, etc.
        return "TEXT"


def merge_match_tables(df_1, df_2, sub1, sub2):
    """Merge two match-table dataframes by columns.

    The predefined shared columns are checked first: if they are present in both
    frames and equal within 1e-13, one copy is kept; if they differ more than
    that, an error is raised. All other columns are treated as unique to each
    input and are suffixed.

    COULD JUST USE df.merge() INSTEAD?
    """

    shared_columns_to_keep = ["RecSes 1", "RecSes 2", "ID1", "ID2", 
                              "ISI_correlations", "ISI_KL_divergence", "ISI_wasserstein_distance",
                              "FR_diff", "ISI_CV_diff",
                              "UID orig 1", "UID orig 2",
                              #"natim_correlations"
                              ]
                              
    tol = 1e-12

    suff_sub1 = f"_{sub1}" if sub1 is not None else ""
    suff_sub2 = f"_{sub2}" if sub2 is not None else ""

    df_1 = df_1.copy()
    df_2 = df_2.copy()

    for col in shared_columns_to_keep:
        if col in df_1.columns and col in df_2.columns:
            if not np.allclose(
                df_1[col].to_numpy(),
                df_2[col].to_numpy(),
                equal_nan=True,
                rtol=0.0,
                atol=tol,
            ):
                raise ValueError(f"Shared column '{col}' differs between inputs.")
            df_2.drop(columns=[col], inplace=True)
        elif col in df_1.columns or col in df_2.columns:
            print(f"Expected shared column '{col}' is missing from one input frame.")

    for col in list(df_1.columns):
        if col not in shared_columns_to_keep:
            df_1.rename(columns={col: f"{col}{suff_sub1}"}, inplace=True)

    for col in list(df_2.columns):
        if col not in shared_columns_to_keep:
            df_2.rename(columns={col: f"{col}{suff_sub2}"}, inplace=True)

    return pd.concat([df_1, df_2], axis=1)


def import_csv_to_sqlite(mt_paths, models, m, p, l):

    db_file = os.path.join(PROJECT_ROOT, "matchtables.db")

    # Define the table name in SQLite
    table_name = f"{m}_{p}_{l}"  # Change this

    print(f"Reading data from the provided paths...")
    df_merged = pd.read_csv(mt_paths[0])
    if len(mt_paths) > 1:
        df_2 = pd.read_csv(mt_paths[1])
        print(f"Merging models {models[0]} and {models[1]}...")
        df_merged = merge_match_tables(df_merged, df_2, models[0], models[1])
        for idx, path in enumerate(mt_paths[2:]):
            if os.path.exists(path):
                df_tomerge = pd.read_csv(path)
                print(f"Merging with model {models[idx + 2]}...")
                df_merged = merge_match_tables(df_merged, df_tomerge, None, models[idx + 2])
            else:
                print(f"Warning: File not found at {path}. Skipping this file.")
    df_merged.rename(columns={'RecSes 1': 'RecSes1', 'RecSes 2': 'RecSes2'}, inplace=True)

    column_defs = []

    print("Generating SQL column definitions from merged DataFrame...")
    for col_name in df_merged.columns:
        dtype = df_merged[col_name].dtype
        sql_type = pandas_to_sqlite_type(dtype)
        # Quote column names to handle spaces or reserved keywords
        safe_col_name = f'"{col_name}"'
        column_defs.append(f"{safe_col_name} {sql_type}")
        print(f"  - Column: '{col_name}', Pandas dtype: {dtype}, SQL type: {sql_type}")

    conn = None
    try:
        # Connect to the SQLite database (creates the file if it doesn't exist)
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        print(f"Successfully connected to database '{db_file}'")

        create_table_sql = (
            f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join(column_defs)});"
        )

        cursor.execute(create_table_sql)
        print(f"Table '{table_name}' checked/created successfully.")

        print(f"Importing data into table '{table_name}'...")
        # Use pandas.DataFrame.to_sql for easy import
        # if_exists='replace': Drops the table first and recreates it. Useful for initial load.
        # if_exists='append': Adds data to the existing table. Use if running multiple times.
        # if_exists='fail': Raises an error if the table already exists.
        df_merged.to_sql(table_name, conn, if_exists="replace", index=False)
        # index=False prevents pandas index from being written as a column

        print("Data imported successfully.")

        # Commit the changes (important!)
        conn.commit()
        print("Changes committed.")

    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
        if conn:
            conn.rollback()  # Roll back changes if anything went wrong
    except pd.errors.EmptyDataError:
        print(
            f"Error: CSV file '{mt_paths[0]}' or '{mt_paths[1]}' or another one is empty or not found."
        )
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        if conn:
            conn.rollback()
    finally:
        # Close the connection
        if conn:
            conn.close()
            print("Database connection closed.")


if __name__ == "__main__":
    data_root = PROJECT_ROOT
    # models = ["DeepUnitMatch", "UMPy"]
    models = ["DeepUnitMatch",
              "UMPy", 
              "DUM_maxdist=20", "DUM_maxdist=50", "DUM_maxdist=100", "DUM_maxdist=inf", 
              "UMPy_maxdist=20", "UMPy_maxdist=50", "UMPy_maxdist=100", "UMPy_maxdist=inf",
              "DUM_W_ij=1","DUM_W_ij=5","DUM_W_ij=10","DUM_W_ij=15","DUM_W_ij=20", 
              "n_output=8_after_ae_and_finetune", "n_output=32_after_ae_and_finetune", "n_output=128_after_ae_and_finetune", "n_output=256_after_ae_and_finetune", 
              "EMD", "DANT", "DANT_no_functional"
              ]

    for mouse in os.listdir(data_root):
        if os.path.isdir(os.path.join(data_root, mouse)):
            print(f"Processing mouse: {mouse}")
            for probe in os.listdir(os.path.join(data_root, mouse)):
                print(f"  Processing probe: {probe}")
                for loc in os.listdir(os.path.join(data_root, mouse, probe)):
                    print(f"    Processing location: {loc}")
                    mt_paths = []
                    models_found = []
                    for model in models:
                        mt_path = os.path.join(
                            data_root, mouse, probe, loc, model, "MatchTable.csv"
                        )
                        if not os.path.exists(mt_path):
                            print(f"      Warning: MatchTable.csv not found for {model} at {mt_path}")
                        else:
                            mt_paths.append(mt_path)
                            models_found.append(model)

                    if len(mt_paths) > 0:
                        import_csv_to_sqlite(mt_paths, models_found, mouse, probe, loc)
                    else:
                        print('No valid MatchTable.csv found.')
