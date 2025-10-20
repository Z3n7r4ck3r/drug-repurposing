import pandas as pd
import requests
import io

def get_iuphar_interactions():
    """
    Fetches the IUPHAR/BPS all interaction data.
    """
    url = "https://www.guidetopharmacology.org/DATA/interactions.tsv"
    print("Fetching IUPHAR/BPS interaction data...")
    response = requests.get(url)
    response.raise_for_status()

    # The file is a tsv file
    df = pd.read_csv(io.StringIO(response.text), sep='	')

    print(f"Fetched {len(df)} interactions from IUPHAR/BPS.")
    return df

if __name__ == "__main__":
    interactions_df = get_iuphar_interactions()
    print(interactions_df.head())
