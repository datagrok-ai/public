import pickle
import pandas as pd


def PyKNNApply(
    model: bytes,
    df: pd.DataFrame,
    nameskeys: str,
    namesvalues: str
) -> pd.DataFrame:
    """Apply a trained KNN model to a dataset and return predictions."""
    keys = [k for k in nameskeys.split(",") if k]
    values = [v for v in namesvalues.split(",") if v]

    if keys and values and len(keys) == len(values):
        df = df.rename(columns=dict(zip(values, keys)))

    trained_model = pickle.loads(model)
    preds = trained_model.predict(df.to_numpy())

    return pd.DataFrame({"pred": preds})
