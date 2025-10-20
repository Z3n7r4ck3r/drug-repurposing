import pandas as pd
import requests
import io

def get_signor_interactions():
    """
    Fetches the complete SIGNOR dataset for Homo sapiens.
    """
    url = "https://signor.uniroma2.it/getData.php?organism=9606&format=tsv"
    print("Fetching SIGNOR interactions...")
    response = requests.get(url)
    response.raise_for_status()

    # Use io.StringIO to treat the string content as a file
    data = io.StringIO(response.text)

    df = pd.read_csv(data, sep='	')
    print(f"Fetched {len(df)} interactions from SIGNOR.")
    return df

if __name__ == "__main__":
    interactions_df = get_signor_interactions()
    print(interactions_df.head())
