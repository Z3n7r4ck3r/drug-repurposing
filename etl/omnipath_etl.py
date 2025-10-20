import omnipath as op
import pandas as pd

def get_omnipath_interactions():
    """
    Fetches signed and directed interactions from OmniPath.
    """
    print("Fetching OmniPath interactions...")
    interactions = op.interactions.OmniPath().get(
        directed=True,
        signed=True,
    )
    print(f"Fetched {len(interactions)} interactions from OmniPath.")
    return interactions

if __name__ == "__main__":
    interactions_df = get_omnipath_interactions()
    print(interactions_df.head())
