from flask import Flask, request, jsonify
from flask_cors import CORS
import subprocess
import os
import logging
import tempfile
import re

app = Flask(__name__)
CORS(app)

logging.basicConfig(
  level=logging.DEBUG,
  format="%(asctime)s - %(levelname)s - %(message)s",
  handlers=[
    logging.StreamHandler(),
    logging.FileHandler("app.log")
  ]
)

@app.route('/predict', methods=['POST'])
def predict():
    """
    Endpoint to run the Boltz prediction command and capture real-time logs.
    """
    try:
        yaml_content = request.json.get('yaml')
        msa_content = request.json.get('msa')

        if not yaml_content:
            return jsonify({"error": "YAML content is required."}), 400

        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".yaml", prefix="input_") as yaml_file:
            yaml_file.write(yaml_content)
            yaml_file_path = yaml_file.name

        msa_file_path = None
        if msa_content:
            with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".a3m", prefix="input_") as msa_file:
                msa_file.write(msa_content)
                msa_file_path = msa_file.name

            with open(yaml_file_path, 'r') as file:
                yaml_content = file.read()

            yaml_content = re.sub(r"msa:\s*[\w\.-]+", f"msa: {os.path.basename(msa_file_path)}", yaml_content)

            with open(yaml_file_path, 'w') as file:
                file.write(yaml_content)

        command = [
            "boltz",
            "predict",
            yaml_file_path,
            "--accelerator", "cpu",
            "--output_format", "pdb",
            "--num_workers", "0",
            "--override",
        ]
        
        if not msa_file_path:
            command.append("--use_msa_server")

        logging.info(f"Running command: {' '.join(command)}")

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        for line in process.stdout:
            logging.info(line.strip())

        process.wait()

        os.remove(yaml_file_path)
        if msa_file_path:
            os.remove(msa_file_path)

        if process.returncode == 0:
            yaml_file_name = os.path.splitext(os.path.basename(yaml_file_path))[0]
            # Locate the results folder containing the predictions
            for root, dirs, files in os.walk("."):
                if "results" in root:
                    predictions_folder = os.path.join(root, "predictions", yaml_file_name)
                    if os.path.exists(predictions_folder):
                        pdb_files = [f for f in os.listdir(predictions_folder) if f.endswith(".pdb")]
                        if pdb_files:
                            pdb_file_path = os.path.join(predictions_folder, pdb_files[0])
                            with open(pdb_file_path, 'r') as pdb_file:
                                pdb_content = pdb_file.read()

                            return jsonify({
                                "success": True,
                                "message": "Prediction completed successfully.",
                                "pdb": pdb_content
                            }), 200

            return jsonify({"success": False, "error": "PDB file not found in predictions folder."}), 500
        else:
            return jsonify({"success": False, "error": "Prediction failed."}), 500

    except Exception as e:
        logging.error(f"Unexpected error: {str(e)}")
        return jsonify({"success": False, "error": "An unexpected error occurred.", "details": str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8001)