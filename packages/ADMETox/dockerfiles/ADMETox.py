import os
import csv
import subprocess
import joblib
import torch
import chemprop
import jaqpotpy
from jaqpotpy import jaqpot
from jaqpotpy.models import MolecularModel
import pandas as pd
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor
from flask import Flask, request

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
    command = [
        "chemprop_predict",
        "--test_path", test_path,
        "--checkpoint_path", checkpoint_path,
        "--preds_path", preds_path
    ]

    subprocess.call(command, cwd=current_path)
    
    #with open(preds_path, 'r') as file:
        #content = file.read()
    
    #return content

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
        os.remove(result_path)
        new_columns = {col: f"{col}_{model}" for col in current_df.columns[1:]}
        current_df = current_df.rename(columns=new_columns)
    else:
        predictions, probabilities = make_euclia_predictions(test_data_path, test_model_name, add_probability)
        current_df = pd.DataFrame({f'Y_{model}': predictions})
        if add_probability:
            current_df[f'Y_{model}_probability'] = probabilities
    
    return current_df

def handle_uploaded_file(test_data_path, models, add_probability):
    models_res = models.split(",")
    
    with ThreadPoolExecutor() as executor:
        # Process each model in parallel
        futures = [executor.submit(handle_model, model, test_data_path, add_probability) for model in models_res]

        # Collect results as they become available
        dfs = [future.result() for future in futures]

    final_df = pd.concat(dfs, axis=1).loc[:, ~pd.concat(dfs, axis=1).columns.duplicated()]
    return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
    raw_data = request.data
    test_data_path = 'smiles_data.csv'
    
    with open(test_data_path, 'wb') as file:
        file.write(raw_data)
    
    models = request.args.get('models')
    add_probability = request.args.get('probability', False)
    response = handle_uploaded_file(test_data_path, models, add_probability)
    os.remove(test_data_path)
    return response

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True, debug=False)