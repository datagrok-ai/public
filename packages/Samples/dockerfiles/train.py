import pickle
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier


def PyKNNTrain(
    df: pd.DataFrame,
    predictColumn: str,
    n_neighbors: int,
    weights: str,
    leaf_size: int,
    p: int,
    metric: str,
    algorithm: str
) -> bytes:
    """Train a KNN classifier and return the serialized model."""
    X = df.drop(columns=[predictColumn])
    y = df[predictColumn].to_numpy()

    model = KNeighborsClassifier(
        n_neighbors=n_neighbors,
        weights=weights,
        leaf_size=leaf_size,
        p=p,
        metric=metric,
        algorithm=algorithm
    )
    model.fit(X, y)

    return pickle.dumps(model)
