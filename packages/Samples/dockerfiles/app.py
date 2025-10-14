import logging
import numpy as np
import pandas as pd
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings, get_logger

from train import PyKNNTrain
from predict import PyKNNApply


logging_level = logging.DEBUG
logging.basicConfig(level=logging_level)

settings = Settings(log_level=logging_level)
app = Celery(settings.celery_name, broker=settings.broker_url)
logger = get_logger()

#name: PyKNNTrain
#meta.mlname: PyKNN
#meta.mlrole: train
#description: Custom Python train function for KNN
#language: python
#input: dataframe df
#input: string predictColumn
#input: int n_neighbors=5 {category: FirstParm}
#input: string weights=uniform {category: Parameters; choices: ["uniform", "distance"]}
#input: int leaf_size=30 {category: Parameters}
#input: int p=1 {category: Parameters; range:1-2}
#input: string metric=minkowski {category: Parameters; choices: ["euclidean", "manhattan", "chebyshev", "minkowski"]}
#input: string algorithm=auto {category: Parameters; choices: ["auto","ball_tree", "kd_tree", "brute"]}
#output: blob model
@app.task(name='knn_train', bind=True, base=DatagrokTask)
def knn_train(
    self,
    df: pd.DataFrame,
    predictColumn: str,
    n_neighbors: int,
    weights: str,
    leaf_size: int,
    p: int,
    metric: str,
    algorithm: str
) -> bytes:
    return PyKNNTrain(df, predictColumn, n_neighbors, weights, leaf_size, p, metric, algorithm)


#name: PyKNNApply
#meta.mlname: PyKNN
#meta.mlrole: apply
#description: Custom Python apply function for KNN
#language: python
#input: blob model
#input: dataframe df
#input: string nameskeys [Original features' names]
#input: string namesvalues [New features' names]
#output: dataframe data_out
@app.task(name='knn_apply', bind=True, base=DatagrokTask)
def knn_apply(self, model: bytes, df: pd.DataFrame, nameskeys: str, namesvalues: str) -> pd.DataFrame:
    return PyKNNApply(model, df, nameskeys, namesvalues)



#name: PyKNNIsApplicable
#meta.mlname: PyKNN
#meta.mlrole: isApplicable
#description: Custom Python isApplicable function for KNN
#language: python
#input: dataframe df
#input: string predictColumn
#output: bool result
@app.task(name='knn_is_applicable', bind=True, base=DatagrokTask)
def knn_is_applicable(self, df: pd.DataFrame, predictColumn: str) -> bool:
    numeric = df.select_dtypes(include=np.number).columns.tolist()
    return len(numeric) == df.shape[1]
