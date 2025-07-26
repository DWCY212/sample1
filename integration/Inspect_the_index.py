# Load a small sample of your data to check gene identifiers
import pandas as pd

tcga_file = "~/code/sample1/statistical_analysis/TCGA/tcga_normalized_data.csv"
geo_file = "~/code/sample1/statistical_analysis/GEO/BRCA_geo_data.csv"

tcga_df = pd.read_csv(tcga_file, index_col=0, nrows=5)
geo_df = pd.read_csv(geo_file, index_col=0, nrows=5)

print("TCGA gene identifiers:")
print(tcga_df.index.tolist())

print("\nGEO gene identifiers:")
print(geo_df.index.tolist())
