import os
import logging
import tempfile
from time import time
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

def is_malformed(smiles):
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

def convert_to_smiles(molecule):
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

def parallel_process_smiles(smis):
  with ThreadPoolExecutor(max_workers=8) as executor:
    valid_smiles = list(executor.map(convert_to_smiles, smis))
  return valid_smiles

def make_chemprop_predictions(df_test, checkpoint_path, batch_size=512):
  # Check for GPU availability
  device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

  # Load the model onto the appropriate device
  mpnn = models.MPNN.load_from_checkpoint(checkpoint_path, map_location=device)

  smis = df_test.iloc[:, 0].tolist()
  valid_indices = [i for i, smi in enumerate(smis) if not is_malformed(smi) and smi != '']
  valid_smiles = [smis[i] for i in valid_indices]
  invalid_indices = [i for i in range(len(smis)) if i not in valid_indices]

  if not valid_smiles:
    return np.array([""] * len(smis), dtype=object)

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
    test_preds.insert(index, "")

  return np.array(test_preds, dtype=object)

def generate_fingerprint(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol:
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return list(map(int, fingerprint))
  return None

def calculate_similarity(mean_vector, fingerprint):
  return np.dot(mean_vector, fingerprint) / (norm(mean_vector) * norm(fingerprint))

def calculate_chemprop_probability(smiles_list, model):
  mean_vector = mean_vectors[model]
  return [
    calculate_similarity(mean_vector, generate_fingerprint(smiles))
    if generate_fingerprint(smiles) else 0.0
    for smiles in smiles_list
  ]

def find_model(model):
  files_in_dir = os.listdir()
  model_name = next(
    (file for file in files_in_dir if model.lower() in file.lower() and file.lower().endswith('.ckpt')),
    None
  )
  return model_name

def predict_for_model(model, df_test, add_probability):
  model_name = find_model(model)
  if not model_name:
    raise ValueError(f"No matching model extension found for model '{model}'")

  start = time()
  predictions = make_chemprop_predictions(df_test, model_name)
  logger.debug(f'Chemprop prediction for {model} took {time() - start}')
  df = pd.DataFrame(predictions, columns=[model])
  if add_probability:
    probabilities = calculate_chemprop_probability(df_test.iloc[:, 0].tolist(), model)
    df[f'Y_{model}_probability'] = probabilities
  return df

def predict(data: str, models: str, add_probability: bool, batch_size=1000):
  models_res = models.split(",")
  result_dfs = []

  df_test = pd.read_csv(
    StringIO(data),
    skip_blank_lines=False,
    keep_default_na=False,
    na_values=['']
  )
  df_test = df_test.fillna('')
  df_test.iloc[:, 0] = parallel_process_smiles(df_test.iloc[:, 0].tolist())

  for model in models_res:
    model_results = []
    for j in range(0, len(df_test), batch_size):
      batch_df = df_test.iloc[j:j + batch_size]
      result_df = predict_for_model(model, batch_df, add_probability)
      model_results.append(result_df)
    model_result_df = pd.concat(model_results, axis=0, ignore_index=True)
    result_dfs.append(model_result_df)

  final_df = pd.concat(result_dfs, axis=1).loc[:, ~pd.concat(result_dfs, axis=1).columns.duplicated()]
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
  return predict(csv, models, False)
