import subprocess
import os
import logging
import threading

from flask import Flask, request, jsonify
from flask_cors import CORS

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

process_lock = threading.Lock()
current_process = None

@app.route('/aizynthfind', methods=['POST'])
def aizynthfind():
  global current_process

  try:
    config_content = request.json.get('config')
    smiles = request.json.get('smiles')

    output_dir = os.getcwd()
    config_path = os.path.join(output_dir, "aizynthcli_data/config.yml")
    smiles_file_path = os.path.join(output_dir, "smiles.txt")

    if config_content:
      with open(config_path, "w") as config_file:
        config_file.write(config_content)

    with open(smiles_file_path, mode='w') as smiles_file:
      smiles_file.write(smiles)
    
    command = [
      "aizynthcli",
      "--config", config_path,
      "--smiles", smiles_file_path,
      "--output", os.path.join(output_dir, "trees.json")
    ]

    logging.info(f"Running command: {' '.join(command)}")

    with process_lock:
      if current_process and current_process.poll() is None:
        logging.info(f"Terminating previous process (PID: {current_process.pid})")
        current_process.terminate()
        current_process.wait()
      
      current_process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
      )

    stdout, stderr = current_process.communicate()

    if current_process.returncode == 0:
      logging.info(f"Process completed successfully. Return code: {current_process.returncode}")
      logging.info(stdout)

      result_file_path = os.path.join(output_dir, "trees.json")
      with open(result_file_path, 'r') as result_file:
        result_content = result_file.read()
      return jsonify({"result": result_content, "success": True}), 200
    else:
      logging.error(f"Process failed. Return code: {current_process.returncode}")
      logging.error(stderr)
      return jsonify({"success": False, "error": stderr}), 500

  except Exception as e:
    logging.exception("Unexpected error occurred")
    return jsonify({"success": False, "error": str(e)}), 500

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=8000)