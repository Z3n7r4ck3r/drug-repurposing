from pandasgwas.get_associations import get_associations
import pandas as pd

def get_gwas_catalog_associations_for_disease(efo_id):
    """
    Fetches associations from the GWAS Catalog for a specific disease.
    """
    print(f"Fetching GWAS Catalog associations for {efo_id}...")
    associations = get_associations(efo_id=efo_id)

    # The 'associations' attribute is a list of dictionaries
    df = pd.DataFrame(associations.associations)

    print(f"Fetched {len(df)} associations from the GWAS Catalog.")
    return df

if __name__ == "__main__":
    # EFO ID for age-related macular degeneration is EFO_0001365
    associations_df = get_gwas_catalog_associations_for_disease("EFO_0001365")
    print(associations_df.head())
