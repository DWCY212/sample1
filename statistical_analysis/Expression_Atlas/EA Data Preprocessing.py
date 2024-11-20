import os
import pandas as pd
import requests

# Base URL for Expression Atlas API
base_url = "https://www.ebi.ac.uk/gxa/json/experiments"

# Retrieve Expression Atlas Data Function
def retrieve_expression_atlas_data(experiment_id, output_file):
    """
    Function to retrieve Expression Atlas data using the API.
    Args:
        experiment_id (str): The experiment ID to retrieve.
        output_file (str): The file path to save the fetched data.
    Returns:
        DataFrame: A Pandas DataFrame containing the retrieved data.
    """
    # Construct the URL for the experiment
    api_url = f"{base_url}/{experiment_id}"
    
    try:
        # Send GET request to the API
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            
            # Print the response keys to understand the structure
            print("Response keys:", data.keys())

            # Extract profiles data
            if 'profiles' in data:
                profiles = data['profiles']
                df = pd.json_normalize(profiles)
                
                # Save raw data to CSV
                df.to_csv(output_file, index=False)
                print(f"Expression Atlas data saved to {output_file}")
                return df
            else:
                raise KeyError(f"'profiles' key not found in response: {data.keys()}")
        
        else:
            raise Exception(f"Failed to retrieve data from {api_url}, Status Code: {response.status_code}")

    except Exception as e:
        raise Exception(f"Error occurred: {e}")

# Clean Data Function
def clean_data(df):
    """
    Function to clean the Expression Atlas data.
    Args:
        df (DataFrame): The DataFrame to clean.
    Returns:
        DataFrame: A cleaned DataFrame.
    """
    # Flatten columns that contain lists by converting them to strings
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, list)).any():
            df[col] = df[col].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)

    # Remove rows with missing values
    df = df.dropna()

    # Remove duplicate entries
    df = df.drop_duplicates()

    return df

# Normalize Data Function
def normalize_data(df):
    """
    Function to normalize the Expression Atlas data.
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
    # Define experiment ID for breast cancer RNA-seq
    experiment_ids = ["E-MTAB-4801"]  # Only using E-MTAB-4801 for breast cancer cell lines
    output_dir = os.path.expanduser("~/code/sample1/statistical_analysis/Expression_Atlas")
    os.makedirs(output_dir, exist_ok=True)
    
    combined_df = pd.DataFrame()
    
    # Retrieve, clean, and normalize the dataset
    for experiment_id in experiment_ids:
        raw_output_file = os.path.join(output_dir, f"{experiment_id}_raw_data.csv")
        cleaned_output_file = os.path.join(output_dir, f"{experiment_id}_cleaned_data.csv")
        normalized_output_file = os.path.join(output_dir, f"{experiment_id}_normalized_data.csv")
        
        # Retrieve data
        df = retrieve_expression_atlas_data(experiment_id, raw_output_file)
        
        # Clean data
        df_cleaned = clean_data(df)
        save_data(df_cleaned, cleaned_output_file)
        
        # Normalize data
        df_normalized = normalize_data(df_cleaned)
        save_data(df_normalized, normalized_output_file)
        
        # Combine datasets (in this case, only one dataset is being used)
        combined_df = pd.concat([combined_df, df_normalized], axis=0)
    
    # Save the combined dataset
    combined_output_file = os.path.join(output_dir, "BRCA_expression_atlas_data.csv")
    save_data(combined_df, combined_output_file)
    
    print(f"Data preprocessing complete. Combined data saved to {combined_output_file}")