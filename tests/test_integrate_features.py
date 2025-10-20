import pandas as pd

from model.integrate_features import integrate_features


def test_integrate_features_combines_tables():
    target_scores = pd.DataFrame({"target_id": ["T1"], "disease_id": ["D1"], "score": [0.8]})
    drug_targets = pd.DataFrame({"drug_id": ["DRUG1"], "target_id": ["T1"], "source": ["DrugCentral"]})
    safety = pd.DataFrame({"drug": ["DRUG1"], "reports": [2]})
    trials = pd.DataFrame({"condition": ["D1"], "active_trials": [3]})

    features = integrate_features(target_scores, drug_targets, safety, trials)
    assert "drug_count" in features.columns
    assert features["drug_count"].iloc[0] == 1
    assert features["active_trials"].iloc[0] == 3
