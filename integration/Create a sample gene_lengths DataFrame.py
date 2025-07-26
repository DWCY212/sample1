import os
import pandas as pd

# Create a sample gene_lengths DataFrame
gene_lengths_data = {
    "Gene": ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"],
    "Length_kb": [2.5, 1.8, 3.2, 4.1, 2.0]
}

gene_lengths_df = pd.DataFrame(gene_lengths_data).set_index("Gene")

# Define output directory and file path
output_dir = os.path.expanduser("~/code/sample1/integration")
os.makedirs(output_dir, exist_ok=True)
gene_lengths_file = os.path.join(output_dir, "gene_lengths.csv")

# Save the DataFrame to CSV
gene_lengths_df.to_csv(gene_lengths_file)

print(f"gene_lengths.csv has been recreated at {gene_lengths_file}")
