import pandas as pd
import requests

def get_openfda_adverse_events(search_term):
    """
    Fetches adverse events from openFDA for a given search term.
    """
    # Base URL for the openFDA drug event API
    base_url = "https://api.fda.gov/drug/event.json"

    # Parameters for the API request
    params = {
        'search': f'patient.reaction.reactionmeddrapt:"{search_term}"',
        'limit': 1000
    }

    print("Fetching openFDA adverse events...")
    response = requests.get(base_url, params=params)
    response.raise_for_status()

    results = response.json().get('results', [])

    # Extract relevant data from the results
    adverse_events = []
    for report in results:
        adverse_events.append({
            'safetyreportid': report.get('safetyreportid'),
            'reactions': [reaction['reactionmeddrapt'] for reaction in report['patient']['reaction']],
            'drugs': [drug.get('medicinalproduct') for drug in report['patient']['drug']]
        })

    df = pd.DataFrame(adverse_events)
    print(f"Fetched {len(df)} adverse events from openFDA.")
    return df

if __name__ == "__main__":
    adverse_events_df = get_openfda_adverse_events("EYE")
    print(adverse_events_df.head())
