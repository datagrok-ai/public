import subprocess
import os
import logging

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

@app.route('/aizynthfind', methods=['POST'])
def aizynthfind():
  try:
    config_content = request.json.get('config')
    smiles = request.json.get('smiles')

    output_dir = os.getcwd()
    config_path = f"{output_dir}/aizynthcli_data/config.yml"
    smiles_file_path = f"{output_dir}/smiles.txt"

    if (config_content):
      with open(config_path, "w") as config_file:
        config_file.write(config_content)

    with open(smiles_file_path, mode='w') as smiles_file:
      smiles_file.write(smiles)
    
    command = [
      "aizynthcli",
      "--config", config_path,
      "--smiles", smiles_file_path,
      "--output", f"{output_dir}/trees.json"
    ]

    logging.info(f"Running command: {' '.join(command)}")

    process = subprocess.Popen(
      command,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True
    )
    
    stdout, stderr = process.communicate()

    if (process.returncode == 0):
      logging.info(f"return code {process.returncode}")
      logging.info(stdout)
      result_file_path = os.path.join(output_dir, 'trees.json')
      with open(result_file_path, 'r') as result_file:
        result_content = result_file.read()
      return jsonify({
        "result": result_content,
        "success": True
      }, 200)
    else:
      logging.info(f"return code {process.returncode}")
      logging.info(stdout)
      logging.info(stderr)
      return jsonify({
        "success": False,
        "error": stderr
      }, 500)
  except Exception as e:
    logging.info(f"Unexpected error: {str(e)}")
    return jsonify({"success": False, "error": str(e)}), 500

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=8000)