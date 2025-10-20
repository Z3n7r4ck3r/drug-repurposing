import pandas as pd
import requests
import gzip
import io

def get_drugcentral_interactions():
    """
    Fetches the DrugCentral drug-target interaction data.
    """
    url = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
    print("Fetching DrugCentral drug-target interaction data...")
    response = requests.get(url)
    response.raise_for_status()

    # The file is a gzipped tsv file
    with gzip.open(io.BytesIO(response.content), 'rt') as f:
        df = pd.read_csv(f, sep='	')

    print(f"Fetched {len(df)} interactions from DrugCentral.")
    return df

if __name__ == "__main__":
    interactions_df = get_drugcentral_interactions()
    print(interactions_df.head())
