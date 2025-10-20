import pandas as pd
import requests
import gzip
import io

def get_string_interactions():
    """
    Fetches the STRING protein-protein interaction data for Homo sapiens.
    """
    # The taxonomic ID for Homo sapiens is 9606
    url = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    print("Fetching STRING interactions...")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    # Use gzip to decompress the file in memory
    with gzip.open(io.BytesIO(response.content), 'rt') as f:
        df = pd.read_csv(f, sep=' ')

    print(f"Fetched {len(df)} interactions from STRING.")
    return df

if __name__ == "__main__":
    interactions_df = get_string_interactions()
    print(interactions_df.head())
