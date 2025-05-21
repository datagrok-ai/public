import os
import signal
import logging
import subprocess
from tempfile import NamedTemporaryFile
from werkzeug.utils import secure_filename
from datagrok_api import DatagrokClient

from multiprocessing import Manager
from flask import Flask, request, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

manager = Manager()
shared_state = manager.dict()

logging.basicConfig(
  level=logging.DEBUG,
  format="%(asctime)s - %(levelname)s - %(message)s",
  handlers=[logging.StreamHandler(), logging.FileHandler("app.log")]
)

def is_process_alive(pid):
  try:
    os.kill(pid, 0)
    return True
  except OSError:
    return False

def terminate_previous_process(user_id):
  process_info = shared_state.get(user_id)
  if not process_info or "pid" not in process_info:
    return False

  pid = process_info["pid"]
  if is_process_alive(pid):
    logging.info(f"Terminating process {pid} for user {user_id}.")
    try:
      os.kill(pid, signal.SIGKILL)
      logging.info(f"Process {pid} terminated.")
    except Exception as e:
      logging.error(f"Failed to terminate process {pid}: {e}")
      return False

  shared_state.pop(user_id, None)
  return True

def run_aizynthfind(config_path, smiles, user_id):
  try:
    with NamedTemporaryFile(delete=False, suffix=".txt") as smiles_file:
      smiles_file.write(smiles.encode())
      smiles_file.flush()

      res_output_dir = os.path.join(os.getcwd(), "aizynthcli_data", user_id)
      
      # Create output directory if it doesn't exist
      os.makedirs(res_output_dir, exist_ok=True)

      command = [
        "aizynthcli",
        "--config", config_path,
        "--smiles", smiles_file.name,
        "--output", f"{res_output_dir}/trees.json"
      ]

    logging.info(f"Executing: {' '.join(command)}")

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    logging.info(f"Started aizynthcli with PID: {process.pid}")

    shared_state[user_id] = {"status": "running", "pid": process.pid}

    stdout, stderr = process.communicate()
    return_code = process.returncode

    if shared_state.get(user_id, {}).get("status") == "terminated":
      logging.warning(f"Process {process.pid} was terminated by a new request.")
      return {"success": False, "error": "Process terminated"}

    shared_state[user_id].update({"status": "completed", "return_code": return_code})

    if return_code == 0:
      with open(f"{res_output_dir}/trees.json", "r") as result_file:
        result_content = result_file.read()
      logging.info("aizynthcli completed successfully")
      return {"result": result_content, "success": True}
    else:
      logging.error(f"aizynthcli failed: {stderr}, error code: {return_code}, stdout: {stdout}")
      return {"success": False, "error": stderr}

  except Exception as e:
    logging.exception(f"Error running aizynthfind for user {user_id}")
    shared_state.pop(user_id, None)
    return {"success": False, "error": str(e)}

@app.route('/aizynthfind', methods=['POST'])
def aizynthfind():
  try:
    output_dir = os.getcwd()
    data = request.json
    user_id = data.get("user_id")
    smiles = data.get("smiles")
    config_dir = data.get("config_dir")

    if not smiles or not user_id:
      return jsonify({"success": False, "error": "Missing 'smiles' in the request"}), 400
    
    config_path = os.path.join(output_dir, "aizynthcli_data", f"{config_dir}", "config.yml") if (config_dir and config_dir != '') else os.path.join(output_dir, "aizynthcli_data", "config.yml")
    if not os.path.exists(config_path):
      return jsonify({"success": False, "error": f"File {config_path} not found"}), 404

    logging.info(f"****************** Reading config")
    with open(config_path, "r") as file:
      content = file.read()
    logging.info(content)


    terminate_previous_process(user_id)
    result = run_aizynthfind(config_path, smiles, user_id)

    if result["success"]:
      return jsonify({"success": True, "result": result["result"]}), 200
    else:
      return jsonify({"success": False, "error": result["error"]}), 500

  except Exception as e:
    logging.exception(f"Error handling request for user {user_id}")
    return jsonify({"success": False, "error": str(e)}), 500


@app.route('/save_config_files', methods=['POST'])
def save_config_files():
  try:
      token = request.json.get('token')
      files = request.json.get('files', [])
      folder_name = request.json.get('folder')
      
      if not files:
          return jsonify({"success": False, "error": "No files provided"}), 400
      
      if not folder_name:
          return jsonify({"success": False, "error": "No config folder name provided"}), 400

      if not token:
          return jsonify({ "success": False, "error": "No token provided"}), 400

      api = DatagrokClient(token, 'http://host.docker.internal:8082')
      
      # Create config directory if it doesn't exist
      config_dir = os.path.join(os.getcwd(), "aizynthcli_data", "config")
      if not os.path.exists(config_dir):
          os.makedirs(config_dir)

      # Create config folder if it doesn't exist and clear if exists
      custom_config_dir = os.path.join(config_dir, folder_name)
      if not os.path.exists(custom_config_dir):
          os.makedirs(custom_config_dir)
      else:
          for file in os.listdir(custom_config_dir):
              file_path = os.path.join(custom_config_dir, file)
              try:
                  if os.path.isfile(file_path):
                      os.unlink(file_path)
              except Exception as e:
                  logging.error(f"Error deleting file {file_path}: {str(e)}")

      saved_files = []
      for file_path in files:
          try:
              file_content = api.download_file('System:AppData', file_path)
              
              # Get filename from path
              filename = os.path.basename(file_path)
              
              save_path = os.path.join(custom_config_dir, filename)
              with open(save_path, 'wb') as f:
                  f.write(file_content)
              
              saved_files.append(filename)
              logging.info(f"Successfully saved {filename} to {save_path}")
              
          except Exception as e:
              logging.error(f"Error processing file {file_path}: {str(e)}")
              continue

      if not saved_files:
          return jsonify({"success": False, "error": "Failed to save any files"}), 500

      return jsonify({
          "success": True, "result": saved_files}), 200

  except Exception as e:
      logging.error(f"Error in save_config_files: {str(e)}")
      return jsonify({"success": False, "error": str(e)}), 500

@app.route('/get_configs', methods=['GET'])
def get_configs():
  try:
      config_dir = os.path.join(os.getcwd(), "aizynthcli_data", "config")
      
      if not os.path.exists(config_dir):
          return jsonify({
              "success": True,
              "result": []
          }), 200

      folders = [f for f in os.listdir(config_dir) 
                if os.path.isdir(os.path.join(config_dir, f))]
      
      return jsonify({
          "success": True,
          "result": folders
      }), 200

  except Exception as e:
      logging.error(f"Error in get_user_configs: {str(e)}")
      return jsonify({
          "success": False,
          "error": str(e)
      }), 500


@app.route('/get_user_configs', methods=['POST'])
def get_user_configs():
  try:
    user_id = request.json.get("user_id")
    if not user_id:
      return jsonify({"success": False, "error": "User ID is required"}), 400

    output_dir = os.getcwd()
    user_config_dir = os.path.join(output_dir, "aizynthcli_data", user_id)
      
    # Create directory if it doesn't exist
    os.makedirs(user_config_dir, exist_ok=True)
      
    # Get all files in the directory
    config_files = []
    for file in os.listdir(user_config_dir):
      if os.path.isfile(os.path.join(user_config_dir, file)):
        config_files.append(file)

    return jsonify({"success": True, "configs": config_files}), 200

  except Exception as e:
    logging.error(f"Unexpected error: {str(e)}")
    return jsonify({"success": False, "error": str(e)}), 500


@app.route('/add_user_custom_config', methods=['POST'])
def add_user_custom_config():
  try:  
    config_content = request.json.get('config')
    config_name = request.json.get('config_name')
    user_id = request.json.get("user_id")
    
    if not user_id:
      return jsonify({"success": False, "error": "User ID is required"}), 400
        
    if not config_name or config_name == '':
      return jsonify({"success": False, "error": "Config name is required"}), 400

    if not config_content:
      return jsonify({"success": False, "error": "Config is required"}), 400
        
    output_dir = os.getcwd()
    user_config_path = os.path.join(output_dir, "aizynthcli_data", user_id)
    os.makedirs(user_config_path, exist_ok=True)
    user_config_file = os.path.join(user_config_path, f"{config_name}")

    with open(user_config_file, "w") as config_file:
      config_file.write(config_content)
    
    logging.info(f"Config {config_name} saved successfully for user {user_id}")
    return jsonify({"success": True, "message": f"Config {config_name} uploaded successfully"}), 200
      
  except Exception as e:
      logging.error(f"Error uploading config: {str(e)}")
      return jsonify({"success": False, "error": str(e)}), 500


@app.route('/sync_dir', methods=['POST'])
def add_user_config():
  try:  
    to_dir_name = request.json.get('to_dir_name')
    from_dir_name = request.json.get('from_dir_name')
    token = request.json.get('token')

    if not token:
      return jsonify({ "success": False, "error": "No token provided"}), 400
    
    if not from_dir_name:
      return jsonify({"success": False, "error": "Directory to sync was not provided"}), 400
        
    output_dir = os.getcwd()
    config_path = os.path.join(output_dir, "aizynthcli_data", to_dir_name)

    api = DatagrokClient(token, 'http://host.docker.internal:8082')
    api.sync_dir('System:AppData', from_dir_name, config_path, False)
    
    logging.info(f"Directory {from_dir_name} synchronized successfully with {config_path}")
    return jsonify({"success": True, "message": f"Directory {from_dir_name} synchronized successfully with {config_path}"}), 200
      
  except Exception as e:
      logging.error(f"Error uploading config: {str(e)}")
      return jsonify({"success": False, "error": str(e)}), 500


@app.route('/health_check', methods=['GET'])
def health_check():
  return jsonify({"success": True}), 200

if __name__ == "__main__":
  app.run(host="0.0.0.0", port=8000)