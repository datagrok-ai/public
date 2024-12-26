from flask import Flask, request, jsonify
import pandas as pd
import subprocess
import os
import sys
import logging
import json
import shutil

app = Flask(__name__)

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

@app.route('/run_reinvent', methods=['POST'])
def run_reinvent():
    """
    Endpoint to run the REINVENT process and return processed results.
    """
    try:
        request_data = request.data
        json_data = json.loads(request_data)

        smiles_list = json_data.get('smiles', [])
        logging.debug(f'smiles list: {smiles_list}')
        config_string = json_data.get('config', " ")
        receptor_file_content = json_data.get('receptor', " ")
        adme_models = json_data.get('admeModels', [])

        output_csv = "stage1_1.csv"

        if not smiles_list or not config_string or not receptor_file_content:
            return {"error": "Missing required inputs: 'smiles', 'config', or 'receptor'"}, 400
        
        # Write SMILES to .smi file
        smiles_filename = "smiles.smi"
        with open(smiles_filename, "w") as smiles_file:
            smiles_file.write("\n".join(smiles_list))
        logging.info("SMILES file written: %s", smiles_filename)

        # Write configuration string to .json file
        config_filename = "kras_docking.json"
        with open(config_filename, "w") as config_file:
            config_file.write(config_string)
        logging.info("Config file written: %s", config_filename)

        # Write receptor content to .pdbqt file
        receptor_filename = "kras_receptor.pdbqt"
        with open(receptor_filename, "w") as receptor_file:
            receptor_file.write(receptor_file_content)
        logging.info("Receptor file written: %s", receptor_filename)


        logging.info("Starting the REINVENT process...")
        current_directory = os.getcwd()
        
        files = [f for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]
        pt_files = [f for f in files if f.endswith('.pt')]
        for pt_file in pt_files:
            folder_name = os.path.splitext(pt_file)[0]
            folder_path = os.path.join(current_directory, folder_name)
            os.makedirs(folder_path, exist_ok=True)
            shutil.move(os.path.join(current_directory, pt_file), os.path.join(folder_path, pt_file))
        
        python_path = sys.executable

        # Define configurations
        global_parameters = """
        run_type = "staged_learning"
        tb_logdir = "tb_stage1"
        json_out_config = "_stage1.json"
        device = "cpu"
        """

        prior_filename = f"/{current_directory}/REINVENT4/priors/mol2mol_similarity.prior"
        parameters = f"""
        [parameters]
        prior_file = "{prior_filename}"
        agent_file = "{prior_filename}"
        summary_csv_prefix = "stage1"
        sample_strategy = "multinomial"
        smiles_file = "/{current_directory}/{smiles_filename}"

        batch_size = 1

        use_checkpoint = false
        use_cuda = false
        unique_sequences = true
        randomize_smiles = true
        """

        learning_strategy = """
        [learning_strategy]
        type = "dap"
        sigma = 128
        rate = 0.0001
        """

        stages_template = fstages = f"""
        [[stage]]
        
        max_score = 1.0
        max_steps = 1
        chkpt_file = 'stage1.chkpt'
        
        [stage.scoring]
        type = "arithmetic_mean"
        
        [[stage.scoring.component]]
        [stage.scoring.component.DockStream]
        [[stage.scoring.component.DockStream.endpoint]]
        name = "Docking with Dockstream"
        weight = 1
        params.configuration_path = "/{current_directory}/kras_docking.json"
        params.docker_script_path = "/{current_directory}/DockStream/docker.py"
        params.docker_python_path =   "/opt/conda/envs/DockStream/bin/python"
        transform.type = "reverse_sigmoid"
        transform.high = 0
        transform.low = -10
        transform.k = 0.5
        
        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/hERG/"
        params.rdkit_2d_normalized = true
        params.target_column = "Class"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/BBB/"
        params.rdkit_2d_normalized = true
        params.target_column = "Class"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/CYP1A2-Inhibitor/"
        params.rdkit_2d_normalized = true
        params.target_column = "Activity"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/CYP2C19-Inhibitor/"
        params.rdkit_2d_normalized = true
        params.target_column = "Activity"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/CYP2C9-Inhibitor/"
        params.rdkit_2d_normalized = true
        params.target_column = "Activity"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/CYP2D6-Inhibitor/"
        params.rdkit_2d_normalized = true
        params.target_column = "Activity"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5

        [[stage.scoring.component]]
        [stage.scoring.component.ChemProp]
        
        [[stage.scoring.component.ChemProp.endpoint]]
        name = "ChemProp"
        weight = 0.5
        
        params.checkpoint_dir = "/{current_directory}/CYP3A4-Inhibitor/"
        params.rdkit_2d_normalized = true
        params.target_column = "Activity"
        
        transform.type = "reverse_sigmoid"
        transform.high = 1.0
        transform.low = 0.0
        transform.k = 0.5
        
        """

        # Combine configurations
        config = global_parameters + parameters + learning_strategy + stages_template

        logging.debug("Generated configuration: %s", config)

        # Write configuration to a TOML file
        toml_config_filename = "stage1.toml"
        with open(toml_config_filename, "w") as tf:
            tf.write(config)
        logging.info("Configuration written to %s", toml_config_filename)

        reinvent_help_command = "reinvent --help"
        logging.info("Running REINVENT help command: %s", reinvent_help_command)
        help_result = subprocess.run(reinvent_help_command, shell=True, capture_output=True, text=True)
        
        # Log the output and error (if any) from the help command
        logging.debug("REINVENT help command output: %s", help_result.stdout)
        if help_result.stderr:
            logging.error("REINVENT help command error: %s", help_result.stderr)

        # Run the REINVENT command
        reinvent_command = "reinvent -l stage1.log stage1.toml"
        logging.info("Running REINVENT command: %s", reinvent_command)
        result = subprocess.run(reinvent_command, shell=True, check=True, capture_output=True, text=True)
        logging.debug("REINVENT command output: %s", result.stdout)
        logging.debug("REINVENT command error (if any): %s", result.stderr)

        # Read and process the output CSV
        logging.info("Reading output CSV: %s", output_csv)
        df = pd.read_csv(output_csv)
        #df_cleaned = df[(df["SMILES_state"] != 0)].drop_duplicates(subset=["SMILES"])
        #logging.info("Processed DataFrame with %d records", len(df_cleaned))

        # Return the cleaned DataFrame as JSON
        #return df_cleaned.to_json(orient='records')
        return df.to_json(orient='records')

    except subprocess.CalledProcessError as e:
        logging.error("Subprocess failed: %s", e)
        return {"error": "Subprocess failed", "details": str(e)}, 500
    except FileNotFoundError as e:
        logging.error("File not found: %s", e)
        return {"error": "File not found", "details": str(e)}, 404
    except pd.errors.EmptyDataError as e:
        logging.error("CSV file is empty or invalid: %s", e)
        return {"error": "CSV processing failed", "details": str(e)}, 400
    except Exception as e:
        logging.error("An unexpected error occurred: %s", e)
        return {"error": "Unexpected error", "details": str(e)}, 500

if __name__ == '__main__':
    logging.info("Starting Flask app...")
    app.run(host='0.0.0.0', port=6666, threaded=True)