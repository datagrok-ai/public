import os
import csv
import subprocess
import joblib
import torch
import chemprop
import numpy as np
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

models_extensions = [
    "Pgp-Inhibitor.pkl", "Pgp-Substrate.pt", "Caco2.pt", "Lipophilicity.pt", 
    "Solubility.pt", "PPBR.pt", "VDss.pt", "CYP1A2-Inhibitor.pt", 
    "CYP3A4-Inhibitor.pkl", "CYP3A4-Substrate.pkl", "CYP2C19-Inhibitor.pt", 
    "CYP2C9-Inhibitor.pt", "CYP2C9-Substrate.pt", "CYP2D6-Substrate.pt", 
    "CYP2D6-Inhibitor.pt", "CL-Hepa.pt", "CL-Micro.pt", "Half-Life.pt",
]

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
        "--gpu", "0",
        "--batch_size", "1024",
        "--empty_cache",
        "--num_workers", "0"
    ]

    subprocess.call(command, cwd=current_path)

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
    print(f'Smiles list: {str(smiles_list)}')
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
    smiles = pd.read_csv(test_path, header=None).iloc[:, 0].tolist()
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles if Chem.MolFromSmiles(smile) is not None]
    model(molecules)
    predictions = model.prediction   
    probabilities = model.probability if add_probability else None
    probabilities = [sublist[0] for sublist in probabilities]
    return predictions, probabilities

def handle_model(model, test_data_path, add_probability):
    test_model_name = next((model_ext for model_ext in models_extensions if model in model_ext), None)
    result_path = f'predictions-{model}.csv'
    
    if 'pt' in test_model_name:
        make_chemprop_predictions(test_data_path, test_model_name, result_path)
        current_df = pd.read_csv(result_path)
        smiles_list = current_df.iloc[:, 0].tolist()
        current_df = current_df.drop(current_df.columns[0], axis=1)
        os.remove(result_path)
        new_columns = {col: f"{model}" for col in current_df.columns}
        current_df = current_df.rename(columns=new_columns)
        print(f'addProbability: {add_probability}')
        if add_probability == 'true':
            probabilities = calculate_chemprop_probability(smiles_list, model)
            current_df[f'Y_{model}_probability'] = probabilities
    else:
        predictions, probabilities = make_euclia_predictions(test_data_path, test_model_name, add_probability)
        current_df = pd.DataFrame({f'{model}': predictions})
        print(f'addProbability: {add_probability}')
        print(type(add_probability))
        if add_probability == 'true':
            current_df[f'{model}_probability'] = probabilities
    
    return current_df

def handle_uploaded_file(test_data_path, models, add_probability, batch_size=3):
    models_res = models.split(",")
    result_dfs = []
    dfs = []
    for i in range(0, len(models_res)):
        print(f'model: {models_res[i]}')
        dfs.append(handle_model(models_res[i], test_data_path, add_probability))

    # Process models in batches
    #for i in range(0, len(models_res), batch_size):
        #batch_models = models_res[i:i + batch_size]
        
        #with ThreadPoolExecutor() as executor:
            # Process each model in the current batch in parallel
            #futures = [executor.submit(handle_model, model, test_data_path, add_probability) for model in batch_models]

            # Collect results as they become available
            #dfs = [future.result() for future in futures]

        # Concatenate results of the current batch
        batch_result_df = pd.concat(dfs, axis=1).loc[:, ~pd.concat(dfs, axis=1).columns.duplicated()]
        result_dfs.append(batch_result_df)

    # Concatenate results of all batches
    final_df = pd.concat(result_dfs, axis=1).loc[:, ~pd.concat(result_dfs, axis=1).columns.duplicated()]
    final_df.to_csv(index=False)
    return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
    raw_data = request.data
    test_data_path = 'smiles_data.csv'
    
    with open(test_data_path, 'wb') as file:
        file.write(raw_data)
    
    models = request.args.get('models')
    add_probability = request.args.get('probability', 'false')
    response = handle_uploaded_file(test_data_path, models, add_probability)
    os.remove(test_data_path)
    return response

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True, debug=False)