# example: search for breast cancer datasets in GEO
from Bio import Entrez

# Set up Entrez with your email
Entrez.email = "dwcy212@163.com"

def search_geo_breast_cancer():
    """
    Function to search GEO for breast cancer datasets that can potentially be integrated with TCGA.
    """
    # Define the search term for breast cancer
    search_term = "breast cancer AND RNA-seq[All Fields] OR microarray[All Fields]"

    # Search GEO for datasets matching the term
    handle = Entrez.esearch(db="gds", term=search_term, retmax=20)  # Retrieving the top 20 datasets, can adjust this as needed
    record = Entrez.read(handle)
    id_list = record['IdList']

    if len(id_list) > 0:
        for gse_id in id_list:
            handle = Entrez.esummary(db="gds", id=gse_id)
            gse_summary = Entrez.read(handle)
            print(f"ID: {gse_id}, Title: {gse_summary[0]['title']}, Samples: {gse_summary[0]['n_samples']}")
    else:
        print("No GEO dataset found for the given search term.")

# Execute the search
search_geo_breast_cancer()