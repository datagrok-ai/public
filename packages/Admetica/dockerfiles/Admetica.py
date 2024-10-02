import os
import csv
import subprocess
import joblib
import torch
import chemprop
import numpy as np
import logging
from time import time
from numpy.linalg import norm
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ThreadPoolExecutor
from flask import Flask, request
from constants import mean_vectors
from chemprop import data, featurizers, models
from lightning import pytorch as pl
from lightning.pytorch.accelerators import find_usable_cuda_devices
from io import StringIO
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

logging_level = logging.DEBUG
logging.basicConfig(level=logging_level)

models_extensions = [
  "Pgp-Inhibitor.pkl", "Caco2.ckpt", "Lipophilicity.ckpt",
  "Solubility.ckpt", "PPBR.ckpt", "VDss.ckpt", "CYP1A2-Inhibitor.ckpt", "CYP1A2-Substrate.ckpt",
  "CYP3A4-Inhibitor.pkl", "CYP3A4-Substrate.pkl", "CYP2C19-Inhibitor.ckpt", "CYP2C19-Substrate.ckpt",
  "CYP2C9-Inhibitor.ckpt", "CYP2C9-Substrate.ckpt", "CYP2D6-Substrate.ckpt",
  "CYP2D6-Inhibitor.ckpt", "CL-Hepa.ckpt", "CL-Micro.ckpt", "Half-Life.ckpt", "hERG.ckpt", "LD50.ckpt"
]

def is_malformed(smiles):
  try:
    mol = Chem.MolFromSmiles(smiles)
    return mol is None
  except Exception as e:
    logging.error(f"Error validating SMILES '{smiles}': {str(e)}")
    return True

def convert_to_smiles(molecule):
  if "M  END" in molecule:
    try:
      mol = Chem.MolFromMolBlock(molecule)
      return Chem.MolToSmiles(mol) if mol else ''
    except Exception as e:
      logging.error(f"Error converting molblock to SMILES: {str(e)}")
      return None
  return molecule

def parallel_process_smiles(smis):
  with ThreadPoolExecutor(max_workers=8) as executor:
    valid_smiles = list(executor.map(convert_to_smiles, smis))
  return valid_smiles

def make_chemprop_predictions(df_test, checkpoint_path, batch_size=512):
  # Load the model onto GPU
  mpnn = models.MPNN.load_from_checkpoint(checkpoint_path, map_location=torch.device('cuda'))
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

  logging.debug('Check GPU availability')
  logging.debug(torch.cuda.is_available())
  logging.debug(find_usable_cuda_devices(1))
  logging.debug(next(mpnn.parameters()).device)

  with torch.inference_mode():
    trainer = pl.Trainer(
      logger=True,
      enable_progress_bar=True,
      accelerator="gpu",
      devices=find_usable_cuda_devices(1),  # Use one available GPU
      precision=16                          # Enable mixed precision for speedup
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

def handle_model(model, df_test, add_probability):
  test_model_name = next((ext for ext in models_extensions if model in ext), None)
  if not test_model_name:
    raise ValueError(f"No matching model extension found for model '{model}'")

  if 'ckpt' in test_model_name:
    start = time()
    predictions = make_chemprop_predictions(df_test, test_model_name)
    logging.debug(f'Chemprop prediction for {model} took {time() - start}')
    df = pd.DataFrame(predictions, columns=[model])
    if add_probability:
      probabilities = calculate_chemprop_probability(df_test.iloc[:, 0].tolist(), model)
      df[f'Y_{model}_probability'] = probabilities
  else:
    predictions, probabilities = make_euclia_predictions(df_test, test_model_name, add_probability)
    df = pd.DataFrame({f'{model}': predictions})
    if add_probability:
      df[f'{model}_probability'] = probabilities

  return df

def handle_uploaded_data(data, models, add_probability, batch_size=1000):
  models_res = models.split(",")
  result_dfs = []

  df_test = pd.read_csv(
    StringIO(data.decode('utf-8')),
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
      result_df = handle_model(model, batch_df, add_probability)
      model_results.append(result_df)
    model_result_df = pd.concat(model_results, axis=0, ignore_index=True)
    result_dfs.append(model_result_df)

  final_df = pd.concat(result_dfs, axis=1).loc[:, ~pd.concat(result_dfs, axis=1).columns.duplicated()]
  return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
  raw_data = request.data

  models = request.args.get('models')
  add_probability = request.args.get('probability', 'false') == 'true'
  start = time()
  response = handle_uploaded_data(raw_data, models, add_probability)
  logging.debug(f'Time required: {time() - start}')
  return response

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=6678, debug=False)