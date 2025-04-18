import pandas as pd
import glob
import sys

def merge_transposed_dataframes(file_paths, columns_to_keep):
    """
    Merges a list of CSV files, ensuring that:
    - Files with only two columns and multiple rows are transposed, making the first column the header
      and the second column the only row.
    - Only specified columns are retained for merging.

    Parameters:
    - file_paths (list): List of file paths to CSV files.
    - columns_to_keep (list): List of column names to retain after transposition.

    Returns:
    - pd.DataFrame: A merged DataFrame containing the specified columns from all valid files.
    """
    dataframes = []

    for file in file_paths:
        df = pd.read_csv(file)

        # Transpose if the dataframe has only two columns and multiple rows
        if df.shape[1] == 2 and df.shape[0] > 1:
            df = df.set_index(df.columns[0]).transpose().reset_index(drop=True)

        # Keep only the desired columns if they exist
        if all(col in df.columns for col in columns_to_keep):
            dataframes.append(df[columns_to_keep])

    # Merge non-empty DataFrames
    return pd.concat(dataframes, axis=0).reset_index(drop=True)

def main(file_paths, output):
    columns_to_keep = ['query', 'p_value', 'source', 'term_name', 'highlighted', 'term_size', 'term_id']
    merged_df = merge_transposed_dataframes(file_paths, columns_to_keep)
    merged_df.to_csv(output, index=False)
    print(f"Merged data saved to {output}")

if __name__ == "__main__":
    file_paths = glob.glob(sys.argv[1]+"/*.csv")  # Update this path as needed
    output = sys.argv[2]  # Define your output file
    main(file_paths, output)
