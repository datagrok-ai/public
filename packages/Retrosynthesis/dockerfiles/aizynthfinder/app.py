import subprocess
import os
import logging
from flask import Flask, request, jsonify
from flask_cors import CORS
import signal
import time
from multiprocessing import Manager

app = Flask(__name__)
CORS(app)

# Initialize multiprocessing manager
manager = Manager()
shared_state = manager.dict()
DEBOUNCE_WINDOW = 2.0  # Seconds between allowed requests


def terminate_previous_process(user_id):
    """Terminate any previous running process for this user with robust error handling"""
    if user_id not in shared_state:
        logging.info(f"No existing process found for user {user_id}")
        return False

    process_info = shared_state[user_id]
    if 'pid' not in process_info:
        logging.info(f"No PID found in process info for user {user_id}")
        return False

    pid = process_info['pid']
    
    # First check if process exists
    try:
        # Cross-platform process existence check
        if os.name == 'posix':  # Unix/Linux/Mac
            if not os.path.exists(f"/proc/{pid}"):
                logging.info(f"Process {pid} already terminated naturally")
                shared_state[user_id]['status'] = 'completed'
                return True
        else:  # Windows
            import ctypes
            PROCESS_QUERY_INFORMATION = 0x0400
            process = ctypes.windll.kernel32.OpenProcess(PROCESS_QUERY_INFORMATION, False, pid)
            if not process:
                logging.info(f"Process {pid} already terminated naturally")
                shared_state[user_id]['status'] = 'completed'
                return True
            ctypes.windll.kernel32.CloseHandle(process)
    except Exception as e:
        logging.warning(f"Error checking process {pid} status: {str(e)}")
        return False

    # Try graceful termination
    try:
        logging.info(f"Sending SIGTERM to process {pid}")
        os.kill(pid, signal.SIGTERM)
        
        # Wait for process to terminate (with timeout)
        for _ in range(5):  # 5 attempts with 0.1s delay
            if os.name == 'posix':
                if not os.path.exists(f"/proc/{pid}"):
                    break
            else:
                if not ctypes.windll.kernel32.OpenProcess(PROCESS_QUERY_INFORMATION, False, pid):
                    break
            time.sleep(0.1)
        else:
            logging.warning(f"Process {pid} didn't terminate, sending SIGKILL")
            os.kill(pid, signal.SIGKILL)
            
        shared_state[user_id].update({
            'status': 'terminated',
            'termination_time': time.time(),
            'termination_method': 'SIGKILL' if _ == 4 else 'SIGTERM'
        })
        logging.info(f"Successfully terminated process {pid}")
        return True
        
    except ProcessLookupError:
        logging.info(f"Process {pid} already terminated")
        shared_state[user_id]['status'] = 'completed'
        return True
    except Exception as e:
        logging.error(f"Error terminating process {pid}: {str(e)}")
        return False

def run_aizynthfind(config_content, smiles, user_id):
    try:
        output_dir = os.getcwd()
        config_path = os.path.join(output_dir, "aizynthcli_data/config.yml")
        smiles_file_path = os.path.join(output_dir, "smiles.txt")
        
        # Write input files
        if config_content:
            with open(config_path, "w") as config_file:
                config_file.write(config_content)
                logging.info(f"Config written to {config_path}")
        
        with open(smiles_file_path, mode='w') as smiles_file:
            smiles_file.write(smiles)
            logging.info(f"SMILES written to {smiles_file_path}")

        command = [
            "aizynthcli",
            "--config", config_path,
            "--smiles", smiles_file_path,
            "--output", os.path.join(output_dir, "trees.json")
        ]
        logging.info(f"Executing command: {' '.join(command)}")

        process = subprocess.Popen(
            command, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True
        )
        logging.info(f"Process started with PID: {process.pid}")

        # Update shared state
        shared_state[user_id] = {
            'status': 'running',
            'pid': process.pid,
            'command': ' '.join(command),
            'start_time': time.time(),
            'last_request_time': time.time()
        }
        logging.info(f"Updated shared state for user {user_id}")

        stdout, stderr = process.communicate()
        return_code = process.returncode
        logging.info(f"Process {process.pid} completed with return code {return_code}")

        # Check if terminated
        if shared_state.get(user_id, {}).get('status') == 'terminated':
            termination_info = {
                'success': False,
                'error': 'Process terminated by new request',
                'terminated': True,
                'pid': process.pid,
                'runtime': time.time() - shared_state[user_id]['start_time']
            }
            logging.warning(f"Process {process.pid} was terminated prematurely")
            return termination_info

        # Update completion status
        shared_state[user_id].update({
            'status': 'completed',
            'return_code': return_code,
            'end_time': time.time(),
            'runtime': time.time() - shared_state[user_id]['start_time']
        })
        logging.info(f"Process {process.pid} completed successfully")

        if return_code == 0:
            result_file_path = os.path.join(output_dir, "trees.json")
            with open(result_file_path, 'r') as result_file:
                result_content = result_file.read()
            return {"result": result_content, "success": True}
        else:
            logging.error(f"Process failed with stderr: {stderr}")
            return {"success": False, "error": stderr}

    except Exception as e:
        logging.exception(f"Unexpected error in run_aizynthfind for user {user_id}")
        if user_id in shared_state:
            shared_state.pop(user_id)
        return {"success": False, "error": str(e)}

@app.route('/aizynthfind', methods=['POST'])
def aizynthfind():
    try:
        config_content = request.json.get('config')
        smiles = request.json.get('smiles')
        user_id = request.json.get('id')
        
        logging.info(f"New request from user {user_id}")

        # Check for existing process
        if user_id in shared_state and shared_state[user_id].get('status') == 'running':
            logging.info(f"Existing running process found for user {user_id}")
            if not terminate_previous_process(user_id):
                return jsonify({
                    "success": False,
                    "error": "Could not terminate existing process",
                    "status": "termination_failed"
                }), 500

        # Run new process
        result = run_aizynthfind(config_content, smiles, user_id)
        
        if result.get("terminated", False):
            logging.warning(f"Process was terminated for user {user_id}")
            return jsonify(result), 500
        elif result["success"]:
            logging.info(f"Request completed successfully for user {user_id}")
            return jsonify(result), 200
        else:
            logging.error(f"Request failed for user {user_id}: {result.get('error')}")
            return jsonify(result), 500

    except Exception as e:
        logging.exception(f"Unexpected error handling request for user {user_id}")
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/status/<user_id>', methods=['GET'])
def get_status(user_id):
    logging.info(f"Status check for user {user_id}")
    if user_id in shared_state:
        status_info = dict(shared_state[user_id])
        logging.debug(f"Returning status info: {status_info}")
        return jsonify(status_info), 200
    logging.warning(f"No status found for user {user_id}")
    return jsonify({"error": "No process found for this user"}), 404

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)