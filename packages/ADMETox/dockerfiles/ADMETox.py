import os
import io
import csv
import subprocess
import multiprocessing
from flask import Flask, request
from jaqpotpy import jaqpot
from jaqpotpy.models import MolecularModel
import joblib
from rdkit import Chem
import pandas as pd
import torch

app = Flask(__name__)

models_extensions = [
    "CYP3A4-Inhibitor.pkl",
    "CL-Hepa.pt",
    "Half-Life.pt",
    "Caco2.pt",
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

def get_current_path():
    return os.path.split(os.path.realpath(__file__))[0]

def get_files_in_current_directory():
    current_directory = os.getcwd()
    files = [f for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]
    files_string = "\n".join(files)
    return files

def read_csv_and_return_string(file_path):
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        csv_string = '\n'.join([','.join(row) for row in csv_reader])

    return csv_string

def make_chemprop_predictions(test_path, checkpoint_path, preds_path):
    current_path = get_current_path()
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

def make_euclia_predictions(batch_df, checkpoint_path):
    #current_path = get_current_path()
    model = joblib.load(checkpoint_path)
    #smiles = pd.read_csv(test_path, header=None).iloc[:, 0].tolist()
    molecules = [Chem.MolFromSmiles(smile) for smile in batch_df if Chem.MolFromSmiles(smile)]
    model(molecules)
    return model.prediction

def handle_uploaded_file(test_data_path, models):
    models_res = [model for model in models.split(",")]
    dfs = []
    all_smiles_df = pd.read_csv(test_data_path, header=None).iloc[:, 0].tolist()
    all_smiles_df = all_smiles_df[1:]
    batch_size = 100
    num_batches = len(all_smiles_df) // batch_size + 1

    dfs = []
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = (batch_num + 1) * batch_size
        batch_df = all_smiles_df[start_idx:end_idx]
        # Initialize an empty list to store DataFrames for the current batch
        batch_dfs = []
        # Loop through models
        for model in models_res:
            test_model_name = next((model_ext for model_ext in models_extensions if model in model_ext), None)
            result_path = f'predictions-{model}_batch_{batch_num}.csv'
            if 'pt' in test_model_name:
                batch_df_path = f'smiles_data_{batch_num}.csv'
                with open(batch_df_path, 'w') as file:
                    file.write('smiles\n')
                    for smiles in batch_df:
                        file.write(f"{smiles}\n")
                make_chemprop_predictions(batch_df_path, test_model_name, result_path)
                current_df = pd.read_csv(result_path)
                new_columns = {col: f"{col}_{model}" for col in current_df.columns[1:]}
                current_df = current_df.rename(columns=new_columns)
                batch_dfs.append(current_df)
            else:
                results = make_euclia_predictions(batch_df, test_model_name)
                current_df = pd.DataFrame({f'Y_{model}': results})
                batch_dfs.append(current_df)
        batch_result_df = pd.concat(batch_dfs, axis=1)
        dfs.append(batch_result_df)
    final_df = pd.concat(dfs, axis=0)
    final_df = final_df.loc[:, ~final_df.columns.duplicated()]
    return final_df.to_csv(index=False)

@app.route('/df_upload', methods=['POST'])
def df_upload():
    raw_data = request.data
    test_data_path = 'smiles_data.csv'
    with open(test_data_path, 'wb') as file:
        file.write(raw_data)
    models = request.args.get('models')
    response = handle_uploaded_file(test_data_path, models)
    return str(torch.cuda.is_available())

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True, debug=False)