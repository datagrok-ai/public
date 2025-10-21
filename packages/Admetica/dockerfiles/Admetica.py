import os
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

def make_chemprop_predictions(smis: List[str], checkpoint_path: str, batch_size: int = 512) -> np.ndarray:
  # Check for GPU availability
  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

  # Load the model onto the appropriate device
  mpnn = models.MPNN.load_from_checkpoint(checkpoint_path, map_location=device)

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

  logger.debug('Check GPU availability')
  logger.debug(f"CUDA available: {torch.cuda.is_available()}")

  if torch.cuda.is_available():
    logger.debug(f"Usable CUDA devices: {find_usable_cuda_devices(1)}")
  else:
    logger.debug("No usable CUDA devices found. Using CPU.")

  logger.debug(f'Model device: {next(mpnn.parameters()).device}')

  with torch.inference_mode():
    trainer = pl.Trainer(
      logger=True,
      enable_progress_bar=True,
      accelerator="gpu" if torch.cuda.is_available() else "cpu",  # Use GPU if available, else CPU
      devices=1,
      precision=16 if torch.cuda.is_available() else 32   # Enable mixed precision if using GPU, otherwise default
    )
    test_preds = trainer.predict(mpnn, test_loader)

  test_preds = [pred.item() for batch in test_preds for pred in batch]
  for index in invalid_indices:
    test_preds.insert(index, np.nan)

  return np.array(test_preds, dtype=object)

def find_model(model: str) -> Optional[str]:
  files_in_dir = os.listdir()
  model_name = next(
    (file for file in files_in_dir if model.lower() in file.lower() and file.lower().endswith('.ckpt')),
    None
  )
  return model_name

def predict_for_model(model: str, smis: List[str]) -> pd.DataFrame:
  model_name = find_model(model)
  if not model_name:
    raise ValueError(f"No matching model extension found for model '{model}'")

  start = time()
  predictions = make_chemprop_predictions(smis, model_name)
  logger.debug(f'Chemprop prediction for {model} took {time() - start}')
  return pd.DataFrame(predictions, columns=[model])

def predict(molecules: pd.Series, models: str, batch_size: int = 1000):
  models_res = models.split(",")
  result_dfs = []
  
  smis = parallel_process_smiles(molecules.fillna('').tolist())
  for model in models_res:
    model_results = []
    for j in range(0, len(smis), batch_size):
      batch_smiles = smis[j:j + batch_size]
      result_df = predict_for_model(model, batch_smiles)
      model_results.append(result_df)
    model_result_df = pd.concat(model_results, axis=0, ignore_index=True)
    result_dfs.append(model_result_df)
    
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