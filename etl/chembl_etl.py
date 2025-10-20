from chembl_webresource_client.new_client import new_client
import pandas as pd

def get_chembl_bioactivity_for_target(target_chembl_id):
    """
    Fetches bioactivity data from ChEMBL for a specific target.
    """
    print(f"Fetching ChEMBL bioactivity data for target {target_chembl_id}...")
    activity = new_client.activity

    # Search for activities related to the target
    res = activity.filter(target_chembl_id=target_chembl_id).only(['molecule_chembl_id', 'activity_type', 'standard_value', 'standard_units', 'pchembl_value'])

    df = pd.DataFrame(res)

    print(f"Fetched {len(df)} bioactivity records from ChEMBL.")
    return df

def get_target_chembl_id(gene_name):
    """
    Gets the ChEMBL ID for a given gene name.
    """
    print(f"Fetching ChEMBL ID for gene {gene_name}...")
    target = new_client.target
    res = target.filter(target_synonym__icontains=gene_name, target_type='SINGLE PROTEIN').only(['target_chembl_id'])

    if res:
        return res[0]['target_chembl_id']
    else:
        return None

if __name__ == "__main__":
    # Example: VEGFA
    gene_name = "VEGFA"
    target_chembl_id = get_target_chembl_id(gene_name)

    if target_chembl_id:
        bioactivity_df = get_chembl_bioactivity_for_target(target_chembl_id)
        print(bioactivity_df.head())
    else:
        print(f"Could not find ChEMBL ID for gene {gene_name}")
