import pandas as pd
import os
import glob

def merge_counts(input_dir, output_file):
    # Get all CSV files in the input directory
    all_files = glob.glob(os.path.join(input_dir, "*.csv"))

    # List to store dataframes
    df_list = []

    # Read each file and extract the required columns (2nd and 3rd)
    for file in all_files:
        print(f"Processing {file}...")
        df = pd.read_csv(file, header=0)  # Read the file with headers

        # Extract the necessary columns (gene_id and sample_id)
        df_filtered = df.iloc[:, [1, 2]]

        # Get the sample name by removing '_counts' from the filename
        sample_name = os.path.basename(file).split('.')[0].replace('_counts', '')

        # Rename columns for clarity
        df_filtered.columns = ['gene_id', sample_name]

        # Append dataframe to the list
        df_list.append(df_filtered)

    # Merge the dataframes on 'Gene_ID' (join by 'Gene_ID' and keep all genes from all files)
    merged_df = df_list[0]  # Start with the first dataframe
    for df in df_list[1:]:
        merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')  # Merge all by 'Gene_ID'

    # Save the merged counts to the output file
    merged_df.to_csv(output_file, index=False)

    print(f"Merged counts saved to {output_file}")

# If this script is called directly, provide the required paths
if __name__ == "__main__":
    import sys
    input_dir = sys.argv[1]  # First argument is input directory
    output_file = sys.argv[2]  # Second argument is output file
    merge_counts(input_dir, output_file)
