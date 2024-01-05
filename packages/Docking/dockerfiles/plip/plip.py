from flask import Flask, request, jsonify
import subprocess
import os
import json

app = Flask(__name__)

@app.route('/process_pdb', methods=['POST'])
def process_pdb():
    try:
        raw_data = request.data
        json_data = json.loads(raw_data)
        # Get the PDB file as a string from the request
        pdb_data = json_data.get('pdb_data', '')

        current_directory = os.getcwd()
        # Construct the pdb_file_path
        pdb_file_path = os.path.join(current_directory, 'src', 'input.pdb')
        # Write content to the file
        with open(pdb_file_path, 'w') as pdb_file:
            pdb_file.write(pdb_data)

        # Substitute with the command you need
        plip_command = [
            'python3', 'src/plip/plip/plipcmd.py',
            '-f', pdb_file_path,
            '-o'
        ]
        subprocess.run(plip_command, shell=True)

        # Read and return the contents of the results file
        # Enter your file path
        results_file_path = ''
        with open(results_file_path, 'r') as results_file:
            results_data = results_file.read()

        return jsonify({'results': results_data})

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
