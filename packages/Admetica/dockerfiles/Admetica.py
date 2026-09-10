import os
# torch>=2.6 defaults torch.load to weights_only=True, which rejects the numpy globals
# in chemprop checkpoints; ours are baked into the image, so full loading is safe.
os.environ.setdefault('TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD', '1')
import logging
import tempfile
import json
from time import time
from typing import List, Optional
from numpy.linalg import norm
from concurrent.futures import ThreadPoolExecutor

import torch
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from lightning import pytorch as pl
from lightning.pytorch.accelerators import find_usable_cuda_devices
from io import StringIO
from celery import Celery
from datagrok_celery_task import DatagrokTask, Settings, get_logger

from constants import mean_vectors
from chemprop import data, featurizers, models


logging_level = logging.DEBUG
logging.basicConfig(level=logging_level)

settings = Settings(log_level=logging_level)
app = Celery(settings.celery_name, broker=settings.broker_url)
logger = get_logger()

# Global flag to control whether exceptions should be raised or logged
raise_ex_flag = False  # Default is to log exceptions, not raise them

def is_malformed(smiles: str) -> bool:
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp_file:
        stderr_fd = 2  # file descriptor for stderr
        stderr_backup = os.dup(stderr_fd)

        try:
            os.dup2(tmp_file.fileno(), stderr_fd)
            mol = Chem.MolFromSmiles(smiles)
        finally:
            os.dup2(stderr_backup, stderr_fd)
            os.close(stderr_backup)

        tmp_file.seek(0)
        warning_msg = tmp_file.read().strip()
    
    logger.debug(f"Checking SMILES: {smiles}, Warning: {warning_msg}")

    if mol is None or warning_msg:
        print(f"Invalid SMILES detected: {smiles}. Warning: {warning_msg}")
        if raise_ex_flag:
            raise ValueError(f"Invalid SMILES string: {smiles}. Warning: {warning_msg}")
        return True

    return False

def convert_to_smiles(molecule: str) -> Optional[str]:
  if "M  END" in molecule:
    try:
      mol = Chem.MolFromMolBlock(molecule)
      return Chem.MolToSmiles(mol) if mol else ''
    except Exception as e:
      logger.error(f"Error converting molblock to SMILES: {str(e)}")
      if raise_ex_flag:
        raise ValueError("Error converting molblock to SMILES") from e
      return None
  return molecule

def parallel_process_smiles(smis: List[str]) -> List[Optional[str]]:
  with ThreadPoolExecutor(max_workers=8) as executor:
    valid_smiles = list(executor.map(convert_to_smiles, smis))
  return valid_smiles

def make_chemprop_predictions(mpnn, trainer, smis: List[str], batch_size: int = 512) -> np.ndarray:
  valid_indices = [i for i, smi in enumerate(smis) if not is_malformed(smi) and smi != '']
  valid_smiles = [smis[i] for i in valid_indices]
  invalid_indices = [i for i in range(len(smis)) if i not in valid_indices]

  if not valid_smiles:
    return np.full(len(smis), np.nan, dtype=float)

  test_data = [data.MoleculeDatapoint.from_smi(smi) for smi in valid_smiles]
  featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
  test_dset = data.MoleculeDataset(test_data, featurizer=featurizer)

  test_loader = data.build_dataloader(
    test_dset,
    shuffle=False,
    batch_size=batch_size
  )

  with torch.inference_mode():
    test_preds = trainer.predict(mpnn, test_loader)

  test_preds = [pred.item() for batch in test_preds for pred in batch]
  for index in invalid_indices:
    test_preds.insert(index, np.nan)

  return np.array(test_preds, dtype=object)

def find_model(model: str) -> Optional[str]:
  target = f'{model.lower()}.ckpt'
  return next((file for file in os.listdir() if file.lower() == target), None)

def predict_for_model(model: str, smis: List[str], chunk_size: int = 1000) -> pd.DataFrame:
  model_name = find_model(model)
  if not model_name:
    raise ValueError(f"No matching model extension found for model '{model}'")

  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
  mpnn = models.MPNN.load_from_checkpoint(model_name, map_location=device)

  logger.debug(f"CUDA available: {torch.cuda.is_available()}")
  if torch.cuda.is_available():
    logger.debug(f"Usable CUDA devices: {find_usable_cuda_devices(1)}")
  logger.debug(f'Model device: {next(mpnn.parameters()).device}')

  trainer = pl.Trainer(
    logger=True,
    enable_progress_bar=True,
    accelerator="gpu" if torch.cuda.is_available() else "cpu",
    devices=1,
    precision=16 if torch.cuda.is_available() else 32
  )

  start = time()
  chunks = [make_chemprop_predictions(mpnn, trainer, smis[j:j + chunk_size])
            for j in range(0, len(smis), chunk_size)]
  predictions = np.concatenate(chunks) if chunks else np.empty(0, dtype=object)
  logger.debug(f'Chemprop prediction for {model} took {time() - start}')
  return pd.DataFrame(predictions, columns=[model])

def predict(molecules: pd.Series, models: str):
  models_res = models.split(",")
  result_dfs = []

  smis = parallel_process_smiles(molecules.fillna('').tolist())
  for model in models_res:
    result_dfs.append(predict_for_model(model, smis))

  final_df = (
    pd.concat(result_dfs, axis=1)
      .loc[:, lambda df: ~df.columns.duplicated()]
      .apply(pd.to_numeric, errors='coerce')
  )
  return final_df


#name: runAdmetica
#meta.cache: all
#meta.cache.invalidateOn: 0 0 1 * *
#input: string csv
#input: string models
#input: bool raiseException = false
#output: dataframe result
@app.task(name='run_admetica', bind=True, base=DatagrokTask)
def run_admetica(self, csv: str, models: str, raiseException: bool=False) -> pd.DataFrame:
  global raise_ex_flag
  raise_ex_flag = raiseException
  
  df = pd.read_csv(
    StringIO(csv),
    skip_blank_lines=False,
    keep_default_na=False,
    na_values=['']
  ).fillna('')
  
  molecules = df.iloc[:, 0]
  return predict(molecules, models)

#name: checkHealth
#output: string result
@app.task(name='check_health', base=DatagrokTask)
def check_health() -> str:
  return json.dumps({"status": "ok"})