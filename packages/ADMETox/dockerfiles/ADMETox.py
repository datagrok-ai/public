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
import jaqpotpy
from jaqpotpy import jaqpot
from jaqpotpy.models import MolecularModel
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ThreadPoolExecutor
from flask import Flask, request
from constants import mean_vectors

app = Flask(__name__)

logging_level = logging.DEBUG
logging.basicConfig(level=logging_level)

models_extensions = [
  "Pgp-Inhibitor.pkl", "Pgp-Substrate.pt", "Caco2.pt", "Lipophilicity.pt", 
  "Solubility.pt", "PPBR.pt", "VDss.pt", "CYP1A2-Inhibitor.pt", 
  "CYP3A4-Inhibitor.pkl", "CYP3A4-Substrate.pkl", "CYP2C19-Inhibitor.pt", 
  "CYP2C9-Inhibitor.pt", "CYP2C9-Substrate.pt", "CYP2D6-Substrate.pt", 
  "CYP2D6-Inhibitor.pt", "CL-Hepa.pt", "CL-Micro.pt", "Half-Life.pt",
]

def is_malformed(smiles):
  try:
    mol = Chem.MolFromSmiles(smiles)
    return mol is None
  except Exception as e:
    print(f"Error validating SMILES '{smiles}': {str(e)}")
    return True

def read_csv_and_return_string(file_path):
  with open(file_path, 'r') as file:
    csv_reader = csv.reader(file)
    csv_string = '\n'.join([','.join(row) for row in csv_reader])

  return csv_string

def make_chemprop_predictions(test_path, checkpoint_path, preds_path):
  current_path = os.path.dirname(os.path.realpath(__file__))
  print(f'Cuda available check: {torch.cuda.is_available()}')

  command = [
    "chemprop_predict",
    "--test_path", test_path,
    "--checkpoint_path", checkpoint_path,
    "--preds_path", preds_path,
    "--batch_size", "1024",
    "--gpu", "0",
    "--empty_cache",
    "--num_workers", "0"
  ]

  subprocess.call(command, cwd=current_path)

  with open(preds_path, 'r') as preds_file:
    preds_data = preds_file.read()
  preds_data = preds_data.replace('Invalid SMILES', '')

  with open(preds_path, 'w') as preds_file:
    preds_file.write(preds_data)

def generate_fingerprint(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is not None:
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return [int(bit) for bit in fingerprint]
  else:
    return None

def calculate_similarity(mean_vector, fingerprint):
  return np.dot(mean_vector, fingerprint) / (norm(mean_vector) * norm(fingerprint))

def calculate_chemprop_probability(smiles_list, model):
  mean_vector = mean_vectors[model]
  probabilities = [
    calculate_similarity(mean_vector, generate_fingerprint(smiles))
    if generate_fingerprint(smiles) is not None
    else 0.0
    for smiles in smiles_list
  ]
  return probabilities

def make_euclia_predictions(test_path, checkpoint_path, add_probability):
  current_path = os.path.dirname(os.path.realpath(__file__))
  model = joblib.load(checkpoint_path)
  smiles_list = pd.read_csv(test_path, header=None).iloc[1:, 0].tolist()
  valid_smiles = [smile for smile in smiles_list if Chem.MolFromSmiles(smile) is not None]
  molecules = [Chem.MolFromSmiles(smile) for smile in valid_smiles]
  model(molecules)
  predictions = model.prediction   
  probabilities = model.probability if add_probability else None
  probabilities = [sublist[0] for sublist in probabilities]

  if len(valid_smiles) < len(smiles_list):
    malformed_indexes = [i for i, smile in enumerate(smiles_list) if smile not in valid_smiles]    
    for index in malformed_indexes:
      predictions.insert(index, " ")
      if add_probability:
        probabilities.insert(index, None)
    
  return predictions, probabilities

def handle_model(model, test_data_path, add_probability):
  test_model_name = next((model_ext for model_ext in models_extensions if model in model_ext), None)
  result_path = f'predictions-{model}.csv'
    
  if 'pt' in test_model_name:
    start = time()
    make_chemprop_predictions(test_data_path, test_model_name, result_path)
    end = time()
    logging.debug(f'Chemprop prediction for {test_model_name} took {end - start}')
    current_df = pd.read_csv(result_path)
    smiles_list = current_df.iloc[:, 0].tolist()
    current_df = current_df.drop(current_df.columns[0], axis=1)
    os.remove(result_path)
    new_columns = {col: f"{model}" for col in current_df.columns}
    current_df = current_df.rename(columns=new_columns)
    if add_probability == 'true':
      probabilities = calculate_chemprop_probability(smiles_list, model)
      current_df[f'Y_{model}_probability'] = probabilities
  else:
    predictions, probabilities = make_euclia_predictions(test_data_path, test_model_name, add_probability)
    current_df = pd.DataFrame({f'{model}': predictions})
    if add_probability == 'true':
      current_df[f'{model}_probability'] = probabilities
    
  return current_df

def handle_uploaded_file(test_data_path, models, add_probability, batch_size=1000):
  models_res = models.split(",")
  result_dfs = []
  with open(test_data_path, 'r') as file:
    lines = file.readlines()
    header = lines[0]
    for model in models_res:
      model_results = []
      for j in range(1, len(lines), batch_size):
        batch_lines = [header] + lines[j:j+batch_size]
        batch_data = ''.join(batch_lines)
        batch_data_file = f'batch_{j}.csv'
        with open(batch_data_file, 'w') as batch_file:
          batch_file.write(batch_data)
        result_df = handle_model(model, batch_data_file, add_probability)
        os.remove(batch_data_file)
        model_results.append(result_df)
      model_result_df = pd.concat(model_results, axis=0, ignore_index=True)
      result_dfs.append(model_result_df)

  final_df = pd.concat(result_dfs, axis=1).loc[:, ~pd.concat(result_dfs, axis=1).columns.duplicated()]
  return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
  raw_data = request.data
  test_data_path = 'smiles_data.csv'
    
  with open(test_data_path, 'wb') as file:
    file.write(raw_data)
    
  models = request.args.get('models')
  add_probability = request.args.get('probability', 'false')
  start = time()
  response = handle_uploaded_file(test_data_path, models, add_probability)
  end = time()
  logging.debug(f'Time required: {end - start}')
  os.remove(test_data_path)
  return response

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=8000, debug=False)