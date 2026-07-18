import sqlite3
import pandas as pd
import os

# Paths relative to the repo root. PROJECT_ROOT is the repo's parent directory,
# which holds the data/results folders alongside this repo.
PROJECT_ROOT = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2"


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


def merge_match_tables(df_um, df_dum):
    """Merge two match-table dataframes by columns.

    Columns that exist in both frames and have identical values are kept once.
    Columns that exist in both frames but differ are renamed to <name>_um and
    <name>_dum so both versions are preserved. Columns unique to one frame are
    kept with their original names.
    """
    df_um = df_um.copy()
    df_dum = df_dum.copy()

    shared_columns = sorted(set(df_um.columns) & set(df_dum.columns))
    for col in shared_columns:
        if df_um[col].equals(df_dum[col]):
            df_dum.drop(columns=[col], inplace=True)
        else:
            df_um.rename(columns={col: f"{col}_um"}, inplace=True)
            df_dum.rename(columns={col: f"{col}_dum"}, inplace=True)

    return pd.concat([df_um, df_dum], axis=1)


def import_csv_to_sqlite(mt_path_um, mt_path_dum, m, p, l):

    db_file = os.path.join(PROJECT_ROOT, "matchtables.db")

    # Define the table name in SQLite
    table_name = f"{m}_{p}_{l}"  # Change this

    print(f"Reading data from '{mt_path_um}' and '{mt_path_dum}'...")
    df_um = pd.read_csv(mt_path_um)
    df_dum = pd.read_csv(mt_path_dum)

    df_merged = merge_match_tables(df_um, df_dum)

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
            f"Error: CSV file '{mt_path_um}' or '{mt_path_dum}' is empty or not found."
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

    for mouse in os.listdir(data_root):
        if os.path.isdir(os.path.join(data_root, mouse)):
            print(f"Processing mouse: {mouse}")
            for probe in os.listdir(os.path.join(data_root, mouse)):
                print(f"  Processing probe: {probe}")
                for loc in os.listdir(os.path.join(data_root, mouse, probe)):
                    print(f"    Processing location: {loc}")
                    mt_path_um = os.path.join(
                        data_root, mouse, probe, loc, "DeepUnitMatch", "MatchTable.csv"
                    )
                    mt_path_dum = os.path.join(
                        data_root, mouse, probe, loc, "UMPy", "MatchTable.csv"
                    )
                    if os.path.exists(mt_path_um) & os.path.exists(mt_path_dum):
                        import_csv_to_sqlite(mt_path_um, mt_path_dum, mouse, probe, loc)
