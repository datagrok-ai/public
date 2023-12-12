import os
import io
import csv
import subprocess
import jaqpotpy
from jaqpotpy import jaqpot
from jaqpotpy.models import MolecularModel
import joblib
import rdkit
from rdkit import Chem
import pandas as pd
from multiprocessing import Pool
from os import path as osp
from flask import Flask, request, jsonify
import torch

app = Flask(__name__)

smiles_global = []
models_extensions = [
    "CYP3A4-Inhibitor.pkl",
    "CL-Hepa.pt",
    "Half-Life.pt",
    "CYP2C9-Inhibitor.pt",
    "CYP2C19-Inhibitor.pt",
    "CYP2D6-Substrate.pt"
]

class ProcessHandler:
    def __init__(self):
        self.pool = None
        self.results = []
    def start_process(self, target, args=(), kwargs={}):
        if self.pool is not None:
            raise RuntimeError("A process is already running.")
        self.results = []
        self.pool = multiprocessing.Pool(processes=1)
        self.results.append(self.pool.apply_async(target, args=args, kwds=kwargs))
    def stop_process(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.join()
            self.pool = None
        else:
            raise RuntimeError("No running process to stop.")
    def get_results(self):
        results = [result.get() for result in self.results]
        return results

def read_csv_and_return_string(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        rows = [row for row in csv_reader]
        csv_string = '\n'.join([','.join(row) for row in rows])

    return csv_string

def make_chemprop_predictions(test_path, checkpoint_path, preds_path):
    current_path = os.path.split(os.path.realpath(__file__))[0]
    command = [
        "chemprop_predict",
        "--test_path", test_path,
        "--checkpoint_path", checkpoint_path,
        "--preds_path", preds_path
    ]

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=current_path)
    output, error = process.communicate()
    with open(preds_path, 'r') as file:
        content = file.read()
    return content

def make_euclia_predictions(test_path, checkpoint_path):
    current_path = os.path.split(os.path.realpath(__file__))[0]
    model = joblib.load(checkpoint_path)
    smiles = pd.read_csv(test_path, header=None).iloc[:, 0].tolist()
    print('smiles')
    print(smiles)
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles if Chem.MolFromSmiles(smile) is not None]
    model(molecules)
    print(model.prediction)
    return model.prediction

def handle_uploaded_file(test_data_path, models):
    models_res = [model for model in models.split(",")]
    dfs = []
    for model in models_res:
        #test_model_name = '{}.pt'.format(model)
        test_model_name = next((model_ext for model_ext in models_extensions if model in model_ext), None)
        print(test_model_name)
        result_path = 'predictions-{}.csv'.format(model)
        if 'pt' in test_model_name:
            make_chemprop_predictions(test_data_path, test_model_name, result_path)
            current_df = pd.read_csv(result_path)
            new_columns = {col: f"{col}_{model}" for col in current_df.columns[1:]}
            current_df = current_df.rename(columns=new_columns)
            dfs.append(current_df)
        else:
            results = make_euclia_predictions(test_data_path, test_model_name)
            current_df = pd.DataFrame({'Y_{}'.format(model): results})
            dfs.append(current_df)
    
    final_df = pd.concat(dfs, axis=1)
    final_df = final_df.loc[:, ~final_df.columns.duplicated()]
    return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
    raw_data = request.data
    test_data_path = 'smiles_data.csv'
    with open(test_data_path, 'wb') as file:
        file.write(raw_data)
    models = request.args.get('models')
    print('models')
    print(models)
    response = handle_uploaded_file(test_data_path, models)
    return str(torch.cuda.is_available())

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True)