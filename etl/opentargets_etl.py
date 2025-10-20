import pandas as pd
import requests
from pathlib import Path
import pyarrow.parquet as pq
from ftplib import FTP

def get_opentargets_associations():
    """
    Fetches the Open Targets gene-disease associations.
    """

    # Creating a directory to store the data
    data_dir = Path("opentargets_data")
    data_dir.mkdir(exist_ok=True)

    # Get the list of files from the FTP server
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd('pub/databases/opentargets/genetics/latest/d2v2g_scored')
    filenames = ftp.nlst()
    ftp.quit()

    all_associations = []

    for filename in filenames:
        if filename.endswith(".parquet"):
            url = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/d2v2g_scored/{filename}"

            print(f"Downloading {filename}...")
            response = requests.get(url)
            response.raise_for_status()

            filepath = data_dir / filename
            with open(filepath, 'wb') as f:
                f.write(response.content)

            print(f"Reading {filename}...")
            df = pd.read_parquet(filepath)
            all_associations.append(df)

    # Concatenate all dataframes
    associations_df = pd.concat(all_associations, ignore_index=True)

    print(f"Fetched {len(associations_df)} associations from Open Targets.")
    return associations_df

if __name__ == "__main__":
    associations_df = get_opentargets_associations()
    print(associations_df.head())
