import os
import pandas as pd
import numpy as np
import rpy2
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer  # Imported KNNImputer

# Activate pandas2ri to enable conversions between pandas DataFrames and R data frames
pandas2ri.activate()

# Import the sva R package
sva = importr('sva')

# Function to perform batch effect correction using ComBat (via rpy2)
def combat_batch_correction(expression_df, batch_labels, mean_only=True):
    # Ensure the number of samples matches the number of batch labels
    assert expression_df.shape[1] == len(batch_labels), "Mismatch between number of samples and batch labels"
    
    # Convert batch labels to R factor
    r_expression_df = pandas2ri.py2rpy(expression_df)
    r_batch_labels = ro.conversion.py2rpy(batch_labels.astype('category'))
    
    # Debugging output
    print("Number of samples in expression_df:", expression_df.shape[1])
    print("Number of genes in expression_df:", expression_df.shape[0])
    print("Number of batch labels:", len(batch_labels))
    print("Batch labels and counts:")
    print(batch_labels.value_counts())
    print("Variance of first 5 genes:")
    print(expression_df.var(axis=1).head())
    
    # Run ComBat with mean.only and mod set to None
    corrected_matrix = sva.ComBat(dat=r_expression_df, batch=r_batch_labels, mod=ro.NULL, mean_only=mean_only)
    
    corrected_df = pandas2ri.rpy2py(corrected_matrix)
    corrected_df.columns = expression_df.columns
    corrected_df.index = expression_df.index

    return corrected_df

def main():
    # Load the combined data from the previous script
    combined_data_file = os.path.expanduser("~/code/sample1/integration/combined_expression_data.csv")
    combined_df = pd.read_csv(combined_data_file, index_col=0)
    
    # Create batch labels for each dataset
    batch_labels = pd.Series(
        ["TCGA"] * len([col for col in combined_df.columns if col.startswith("TCGA")]) +
        ["GEO"] * len([col for col in combined_df.columns if col.startswith("GEO")]),
        index=combined_df.columns
    )
    
    # Filter out batches with fewer than two samples
    valid_batches = batch_labels.value_counts()[batch_labels.value_counts() > 1].index
    valid_samples = batch_labels[batch_labels.isin(valid_batches)].index
    
    # Filter the combined dataframe and batch labels
    filtered_combined_df = combined_df[valid_samples]
    filtered_batch_labels = batch_labels[valid_samples]
    
    # Ensure all data is numeric
    filtered_combined_df = filtered_combined_df.apply(pd.to_numeric, errors='coerce')
    
    # Replace infinite values with NaN (if any)
    filtered_combined_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # **Impute missing values using KNNImputer (if any missing values remain)**
    if filtered_combined_df.isnull().values.any():
        imputer = KNNImputer(n_neighbors=5)
        imputed_data = imputer.fit_transform(filtered_combined_df.T)  # Transpose to impute samples
        filtered_combined_df = pd.DataFrame(imputed_data.T, index=filtered_combined_df.index, columns=filtered_combined_df.columns)
    
    # Apply log transformation to ensure variability
    filtered_combined_df = np.log1p(filtered_combined_df)
    
    # Add a small constant to ensure numerical stability
    filtered_combined_df += 1e-6
    
    # Ensure each batch has more than three samples
    valid_batches = filtered_batch_labels.value_counts()[filtered_batch_labels.value_counts() > 3].index
    filtered_combined_df = filtered_combined_df.loc[:, filtered_batch_labels.isin(valid_batches)]
    filtered_batch_labels = filtered_batch_labels[filtered_batch_labels.isin(valid_batches)]
    
    # Ensure filtered_combined_df and filtered_batch_labels have matching columns and labels
    filtered_combined_df = filtered_combined_df.loc[:, filtered_batch_labels.index]
    filtered_batch_labels = filtered_batch_labels.loc[filtered_combined_df.columns]
    
    # Additional Debugging Information
    print("Checking for all-zero rows and columns before ComBat...")
    zero_variance_rows = filtered_combined_df.var(axis=1) == 0
    zero_variance_columns = filtered_combined_df.var(axis=0) == 0
    
    if zero_variance_rows.any():
        print(f"Warning: Found {zero_variance_rows.sum()} rows with zero variance. These rows will be removed.")
        filtered_combined_df = filtered_combined_df.loc[~zero_variance_rows]
    
    if zero_variance_columns.any():
        print(f"Warning: Found {zero_variance_columns.sum()} columns with zero variance. These columns will be removed.")
        filtered_combined_df = filtered_combined_df.loc[:, ~zero_variance_columns]
    
    # Ensure filtered_combined_df and filtered_batch_labels have matching columns and labels
    filtered_combined_df = filtered_combined_df.loc[:, filtered_batch_labels.index]
    filtered_batch_labels = filtered_batch_labels.loc[filtered_combined_df.columns]
    
    # Apply ComBat for batch effect correction
    try:
        # Debugging the shape of data and batch labels before correction
        print("Filtered combined dataframe shape (genes x samples):", filtered_combined_df.shape)
        print("Filtered batch labels length:", len(filtered_batch_labels))
        assert filtered_combined_df.shape[1] == len(filtered_batch_labels), "Mismatch between number of samples and batch labels"
        corrected_combined_df = combat_batch_correction(filtered_combined_df, filtered_batch_labels, mean_only=True)
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        print("RRuntimeError encountered during ComBat batch correction. Possible causes:")
        print("1. Singular design matrix (e.g., insufficient variability in the data).")
        print("2. Too few samples in one or more batches.")
        print("3. Data preprocessing issues leading to near-zero variance.")
        print("Batch labels and their counts:")
        print(filtered_batch_labels.value_counts())
        print("Error details:", str(e))
        # Attempt to retry with modified parameters
        print("Attempting retry with modified parameters (mean_only=False)...")
        try:
            # Rechecking batch labels and filtered_combined_df
            print("Filtered combined dataframe shape (genes x samples) before retry:", filtered_combined_df.shape)
            print("Filtered batch labels length before retry:", len(filtered_batch_labels))
            corrected_combined_df = combat_batch_correction(filtered_combined_df, filtered_batch_labels, mean_only=False)
        except rpy2.rinterface_lib.embedded.RRuntimeError as e:
            print("Retry failed. Consider adjusting batch sizes or further preprocessing to ensure sufficient variability.")
            corrected_combined_df = None
    
    # Save the final integrated expression matrix
    output_dir = os.path.expanduser("~/code/sample1/integration")  # Updated output directory
    os.makedirs(output_dir, exist_ok=True)
    corrected_output_file = os.path.join(output_dir, "batch_corrected_expression_matrix.csv")
    if corrected_combined_df is not None:
        corrected_combined_df.to_csv(corrected_output_file)
        print(f"Batch effect correction complete. Corrected data saved to {corrected_output_file}")
    else:
        print("Batch correction failed. No output to save.")
    
    # PCA for Visualization
    if corrected_combined_df is not None:
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(corrected_combined_df.T)
    
        # Plot the first two principal components
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(
            principal_components[:, 0],
            principal_components[:, 1],
            c=filtered_batch_labels.factorize()[0],
            cmap='viridis',
            alpha=0.5
        )
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA of Integrated Batch-Corrected Data')
        plt.legend(handles=scatter.legend_elements()[0], labels=filtered_batch_labels.unique())
        plt.show()
    
        # Boxplots for Visualization
        corrected_combined_df.boxplot(figsize=(20, 6), grid=False)
        plt.title('Boxplot of Normalized Expression Data Across Samples')
        plt.xlabel('Samples')
        plt.ylabel('Normalized Expression')
        plt.show()
    else:
        print("Batch correction failed. No PCA or boxplot will be generated.")

if __name__ == "__main__":
    main()
