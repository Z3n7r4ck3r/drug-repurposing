import pandas as pd
import requests
import zipfile
import io

def get_hpa_protein_atlas_data():
    """
    Fetches the Human Protein Atlas data.
    """
    url = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
    print("Fetching HPA proteinatlas.tsv data...")
    response = requests.get(url)
    response.raise_for_status()

    # The file is a zip archive, so we need to extract it in memory
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        # Find the .tsv file in the zip archive
        tsv_file_name = [name for name in z.namelist() if name.endswith('.tsv')][0]
        with z.open(tsv_file_name) as f:
            df = pd.read_csv(f, sep='	')

    print(f"Fetched {len(df)} records from the Human Protein Atlas.")
    return df

def parse_and_filter_ocular_tissues(df):
    """
    Parses the 'RNA tissue specific nTPM' column and filters for ocular tissues.
    """
    ocular_tissues = ['retina', 'choroid', 'cornea', 'optic nerve', 'uvea', 'lacrimal gland']

    # The 'RNA tissue specific nTPM' column is a string of 'tissue:value;tissue:value;...'
    # We need to parse this string to extract the expression values for each tissue.

    # Create a new DataFrame to store the parsed data
    expression_data = []

    for index, row in df.iterrows():
        gene = row['Gene']
        expression_string = row['RNA tissue specific nTPM']

        if pd.isna(expression_string):
            continue

        tissues = expression_string.split(';')
        for tissue in tissues:
            try:
                tissue_name, value = tissue.split(':')
                if tissue_name.lower() in ocular_tissues:
                    expression_data.append({
                        'Gene': gene,
                        'Tissue': tissue_name,
                        'nTPM': float(value)
                    })
            except ValueError:
                # Handle cases where the split doesn't produce two values
                continue

    ocular_df = pd.DataFrame(expression_data)
    print(f"Found {len(ocular_df)} records for ocular tissues.")
    return ocular_df

if __name__ == "__main__":
    hpa_df = get_hpa_protein_atlas_data()
    ocular_expression_df = parse_and_filter_ocular_tissues(hpa_df)
    print(ocular_expression_df.head())
