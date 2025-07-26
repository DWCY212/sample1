import os
import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer  # Imported KNNImputer

# Function to load data from each dataset
def load_data(file_path, prefix):
    df = pd.read_csv(file_path, index_col=0)
    df.columns = [f"{prefix}_{col}" for col in df.columns]
    return df

# Normalize Data Function (TPM)
def normalize_data_tpm(df, gene_lengths):
    """
    Function to normalize the data using TPM.
    Args:
        df (DataFrame): The DataFrame to normalize.
        gene_lengths (Series): A Series containing gene lengths in kilobases.
    Returns:
        DataFrame: A normalized DataFrame with TPM values.
    """
    # Calculate RPK (Reads Per Kilobase)
    rpk = df.div(gene_lengths, axis=0)
    
    # Calculate scaling factor per sample (sum of all RPKs)
    scaling_factors = rpk.sum(axis=0) / 1e6
    scaling_factors = scaling_factors.replace(0, 1e-6).infer_objects(copy=False)  # Prevent division by zero
    
    # Calculate TPM
    tpm = rpk.div(scaling_factors, axis=1)
    
    return tpm

def main():
    # Define paths for preprocessed TCGA and GEO datasets
    tcga_file = os.path.expanduser("~/code/sample1/statistical_analysis/TCGA/tcga_normalized_data.csv")
    geo_file = os.path.expanduser("~/code/sample1/statistical_analysis/GEO/BRCA_geo_data.csv")
    
    # Load datasets with unique column names
    tcga_df = load_data(tcga_file, "TCGA")
    geo_df = load_data(geo_file, "GEO")
    
    # Remove duplicate rows if any exist
    tcga_df = tcga_df.loc[~tcga_df.index.duplicated(keep='first')]
    geo_df = geo_df.loc[~geo_df.index.duplicated(keep='first')]
    
    # Load gene lengths for TPM normalization
    gene_lengths_file = os.path.expanduser("~/code/sample1/integration/gene_lengths.csv")
    gene_lengths = pd.read_csv(gene_lengths_file, index_col=0).squeeze()
    
    # Normalize datasets using TPM
    tcga_df = normalize_data_tpm(tcga_df, gene_lengths)
    geo_df = normalize_data_tpm(geo_df, gene_lengths)
    
    # Combine datasets
    combined_df = pd.concat([tcga_df, geo_df], axis=1, join='inner')
    
    # Replace infinite values with NaN (if any arise during normalization)
    combined_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # **Impute missing values using KNNImputer**
    imputer = KNNImputer(n_neighbors=5)
    imputed_data = imputer.fit_transform(combined_df.T)  # Transpose to impute samples
    combined_df_imputed = pd.DataFrame(imputed_data.T, index=combined_df.index, columns=combined_df.columns)
    
    # Save the combined data for the next script
    output_dir = os.path.expanduser("~/code/sample1/integration")  # Updated output directory
    os.makedirs(output_dir, exist_ok=True)
    combined_output_file = os.path.join(output_dir, "combined_expression_data.csv")
    combined_df_imputed.to_csv(combined_output_file)
    
    print(f"Data preprocessing complete. Combined data saved to {combined_output_file}")

if __name__ == "__main__":
    main()
