# for breast cancer RNA-seq ("GSE281571", "GSE281303", "GSE273890", "GSE62944", "GSE42568")
import os
import pandas as pd
from Bio import Entrez
import urllib.request
import gzip

# Set up Entrez for GEO search (replace with your actual email)
Entrez.email = "dwcy212@163.com"

# Retrieve GEO Data Function
def retrieve_geo_data(gse_id, output_file):
    """
    Function to retrieve GEO data using the Entrez API.
    Args:
        gse_id (str): The GEO dataset ID to retrieve.
        output_file (str): The file path to save the fetched data.
    Returns:
        DataFrame: A Pandas DataFrame containing the retrieved data.
    """
    # Construct the FTP URL for downloading the dataset
    ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    
    try:
        # Download the gzipped file using urllib
        gz_output_file = output_file + ".gz"
        urllib.request.urlretrieve(ftp_url, gz_output_file)
        print(f"Downloaded {gse_id} data to {gz_output_file}")
        
        # Unzip the gzipped file and load it into a DataFrame
        with gzip.open(gz_output_file, 'rt') as f:
            df = pd.read_csv(f, sep='\t', comment='!', low_memory=False)  # Use '!' to skip metadata lines

        # Save the uncompressed DataFrame to a CSV file
        df.to_csv(output_file, index=False)
        print(f"Raw data saved to {output_file}")

        return df
    except Exception as e:
        raise Exception(f"Failed to retrieve data from {ftp_url}. Error: {e}")


# Clean Data Function
def clean_data(df):
    """
    Function to clean the GEO data.
    Args:
        df (DataFrame): The DataFrame to clean.
    Returns:
        DataFrame: A cleaned DataFrame.
    """
    # Remove rows with missing values
    df = df.dropna()
    # Remove duplicate entries
    df = df.drop_duplicates()
    return df

# Normalize Data Function
def normalize_data(df):
    """
    Function to normalize the GEO data.
    Args:
        df (DataFrame): The DataFrame to normalize.
    Returns:
        DataFrame: A normalized DataFrame.
    """
    # Example normalization: Min-Max scaling on numeric columns
    numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns
    df[numeric_cols] = (df[numeric_cols] - df[numeric_cols].min()) / (df[numeric_cols].max() - df[numeric_cols].min())
    return df

# Save Data Function
def save_data(df, output_file):
    """
    Function to save the cleaned and normalized data to a CSV file.
    Args:
        df (DataFrame): The DataFrame to save.
        output_file (str): The file path to save the data.
    """
    df.to_csv(output_file, index=False)

# Main Code
if __name__ == "__main__":
    # Define paths and dataset IDs
    gse_ids = ["GSE281571", "GSE281303", "GSE273890", "GSE62944", "GSE42568"]
    output_dir = os.path.expanduser("~/code/sample1/statistical_analysis/GEO")
    os.makedirs(output_dir, exist_ok=True)
    
    combined_df = pd.DataFrame()
    
    # Retrieve, clean, and normalize each dataset
    for gse_id in gse_ids:
        raw_output_file = os.path.join(output_dir, f"{gse_id}_raw_data.csv")
        cleaned_output_file = os.path.join(output_dir, f"{gse_id}_cleaned_data.csv")
        normalized_output_file = os.path.join(output_dir, f"{gse_id}_normalized_data.csv")
        
        # Retrieve data
        df = retrieve_geo_data(gse_id, raw_output_file)
        
        # Clean data
        df_cleaned = clean_data(df)
        save_data(df_cleaned, cleaned_output_file)
        
        # Normalize data
        df_normalized = normalize_data(df_cleaned)
        save_data(df_normalized, normalized_output_file)
        
        # Combine datasets
        combined_df = pd.concat([combined_df, df_normalized], axis=0)
    
    # Save the combined dataset
    combined_output_file = os.path.join(output_dir, "BRCA_geo_data.csv")
    save_data(combined_df, combined_output_file)
    
    print(f"Data preprocessing complete. Combined data saved to {combined_output_file}")