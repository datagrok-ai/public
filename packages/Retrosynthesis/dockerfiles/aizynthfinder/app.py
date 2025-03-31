import os
import signal
import logging
import subprocess
from tempfile import NamedTemporaryFile

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

      command = [
        "aizynthcli",
        "--config", config_path,
        "--smiles", smiles_file.name,
        "--output", "trees.json"
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
      with open("trees.json", "r") as result_file:
        result_content = result_file.read()
      logging.info("aizynthcli completed successfully")
      return {"result": result_content, "success": True}
    else:
      logging.error(f"aizynthcli failed: {stderr}")
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
    user_id = data.get("id")
    smiles = data.get("smiles")
    config_path = os.path.join(output_dir, "aizynthcli_data", "config.yml")

    if not os.path.exists(config_path):
      return jsonify({"success": False, "error": f"Config file not found at {config_path}"}), 400

    if not smiles:
      return jsonify({"success": False, "error": "Missing 'smiles' in the request"}), 400

    terminate_previous_process(user_id)
    result = run_aizynthfind(config_path, smiles, user_id)

    if result["success"]:
      return jsonify({"success": True, "result": result["result"]}), 200
    else:
      return jsonify({"success": False, "error": result["error"]}), 500

  except Exception as e:
    logging.exception(f"Error handling request for user {user_id}")
    return jsonify({"success": False, "error": str(e)}), 500

if __name__ == "__main__":
  app.run(host="0.0.0.0", port=8000)