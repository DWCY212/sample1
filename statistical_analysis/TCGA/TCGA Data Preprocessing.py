# example line27 TCGA-BRCA for breast cancer
# import Libraries
import os # os: provides functions to interact with the operating system, such as handling file paths.
import pandas as pd
import requests # requests: send HTTP requests, making it useful to interact with web APIs.
from io import StringIO # StringIO: work with strings as if they were file-like objects, which is useful when handling data directly in memory.

# retrive TCGA Data Function
def retrieve_tcga_data(api_url, output_file):
    """
    Function to retrieve TCGA data using an API.
    Args:
        api_url: The API URL from which TCGA data will be fetched.
        api_url (str): The API endpoint to fetch data from.
        output_file (str): The file path to save the fetched data.
    Returns:
        DataFrame: A Pandas DataFrame containing the retrieved data.
    """
    headers = {
        "Content-Type": "application/json"
    }
    params = {
        "filters": {
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TCGA-BRCA"]  # Replace with desired project, e.g., TCGA-BRCA for breast cancer
            }
        },
        "format": "TSV",  # You can change to "JSON" if needed
        "fields": "file_id,file_name,cases.submitter_id,cases.samples.sample_type",  # Specify the fields you want
        "size": "100"  # Limit the number of results, adjust as needed
    }
    response = requests.post(api_url, headers=headers, json=params) # Sends a GET request to the specified API URL.
    if response.status_code == 200: # Checks if the response from the API is successful (200 is the HTTP status code for success).
        data = StringIO(response.text) # Converts the response text into a file-like object.
        df = pd.read_csv(data, sep='\t') # Reads the data using Pandas, assuming it’s tab-separated (\t).
        if df.empty:
            raise ValueError("The retrieved data is empty. Check the API or request parameters.")
        df.to_csv(output_file, index=False) # Saves the retrieved data to the specified CSV file.
        print(f"Raw data saved to {output_file}")
        return df
    else:
        raise Exception(f"Failed to retrieve data from {api_url}, Statue Code: {response.status_code}")

# clean Data Function
def clean_data(df):
    """
    Function to clean the TCGA data.
    Args:
        df (DataFrame): The DataFrame to clean.
    Returns:
        DataFrame: A cleaned DataFrame.
    """
    # Remove rows with missing values
    df = df.dropna() # Removes rows that contain missing values.
    # Remove duplicate entries
    df = df.drop_duplicates() # Removes duplicate rows.
    # Further cleaning steps ...
    # ...
    if df.empty:
        raise ValueError("The cleaned data is empty. Please check the cleaning steps.")
    return df

# normalize Data Function
def normalize_data(df):
    """
    Function to normalize the TCGA data.
    Args:
        df (DataFrame): The DataFrame to normalize.
    Returns:
        DataFrame: A normalized DataFrame.
    """
    # Example normalization: Min-Max scaling on numeric columns. Standardizing numeric columns so that values are scaled between 0 and 1.
    numeric_cols = df.select_dtypes(include=['float64', 'int64']).columns # select_dtypes(include=['float64', 'int64']): Selects numeric columns (either float or integer).
    df[numeric_cols] = (df[numeric_cols] - df[numeric_cols].min()) / (df[numeric_cols].max() - df[numeric_cols].min()) # Min-Max scaling: This technique scales the data so that the minimum value becomes 0 and the maximum becomes 1. This is often used to bring all data to a similar range for easier comparison.
    return df

# save Data Function
def save_data(df, output_file):
    """
    Function to save the cleaned and normalized data to a CSV file.
    Args:
        df (DataFrame): The DataFrame to save.
        output_file (str): The file path to save the data.
    """
    df.to_csv(output_file, index=False) # to_csv(): Saves the provided DataFrame to the specified file path in CSV format; index=False: Prevents Pandas from saving the row indices as a separate column in the CSV file.
    print(f"Data saved to {output_file}")

# main code
if __name__ == "__main__":
    # define paths and API.
    api_url = "https://api.gdc.cancer.gov/files" # replace with actual API
    output_dir = os.path.expanduser("~/code/sample1/statistical_analysis/TCGA")
    os.makedirs(output_dir, exist_ok=True)
    raw_output_file = os.path.join(output_dir, "tcga_raw_data.csv")
    cleaned_output_file = os.path.join(output_dir, "tcga_cleaned_data.csv")
    normalized_output_file = os.path.join(output_dir, "tcga_normalized_data.csv")

    # retieve data
    df = retrieve_tcga_data(api_url, raw_output_file)

    # clean data
    df_cleaned = clean_data(df)
    save_data(df_cleaned, cleaned_output_file)

    # normalize data
    df_normalized = normalize_data(df_cleaned)
    save_data(df_normalized, normalized_output_file)

    print(f"Data preprocessing complete. Saved cleaned data to {cleaned_output_file} and normalized data to {normalized_output_file}")