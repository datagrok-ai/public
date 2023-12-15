from flask import Flask, request, jsonify
import subprocess
import os
import tempfile
import json
import hashlib
import re

app = Flask(__name__)

def prepare_autogrid_config(folder_path, receptor_basename, autodock_gpf):
    config = "{}/{}.gpf".format(folder_path, receptor_basename)
    with open(config, "w") as config_file:
        config_file.write(autodock_gpf)
    return config

def calculate_hash(data):
    return hashlib.sha256(data.encode()).hexdigest()

def run_process(command, folder_path, shell=False):
    output_file = 'out.txt'
    error_file = 'err.txt'
    with open(output_file, 'w+') as fout:
        with open(error_file, 'w+') as ferr:
            returncode = subprocess.call(command, stdout=fout, stderr=ferr, shell=shell, cwd=folder_path)
            fout.seek(0)
            output = fout.read()
            ferr.seek(0)
            error = ferr.read()
            os.remove(output_file)
            os.remove(error_file)
    return returncode, output, error

@app.route('/check_opencl', methods=['GET'])
def check_opencl():
    try:
        command = ['clinfo']
        process = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        output, _ = process.communicate()
        if process.returncode == 0:
            return jsonify({'success': True, 'output': output})
        else:
            return jsonify({'success': False, 'error': 'clinfo execution failed'})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/autodock/dock_ligand', methods=['POST'])
def dock():
    raw_data = request.data
    json_data = json.loads(raw_data)

    receptor_value = json_data.get('receptor', '')
    receptor_format = json_data.get('receptor_format', '')
    ligand_value = json_data.get('ligand', '')
    ligand_format = json_data.get('ligand_format', '')
    autodock_gpf = json_data.get('autodock_gpf', '')
    debug_mode = request.args.get('debug', False)

    folder_name = calculate_hash(receptor_value + autodock_gpf)
    folder_path = os.path.join(os.getcwd(), folder_name)

    pattern = r"receptor\s+(\S+\.pdbqt)"
    match = re.search(pattern, autodock_gpf)
    receptor_basename = match.group(1) if match else None

    receptor_name = receptor_basename[:-6]
    ligand_name = 'ligand'

    receptor_path = '{}.{}'.format(receptor_name, receptor_format)
    ligand_path = '{}.{}'.format(ligand_name, ligand_format)
    
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        
        with open('{}/{}'.format(folder_path, receptor_path), 'w') as receptor_file:
            receptor_file.write(receptor_value)

        if 'pdbqt' not in receptor_path:
            subprocess.call(['prepare_receptor4.py', '-r', receptor_path], cwd=folder_path)
        
        autogrid_config = prepare_autogrid_config(folder_path, receptor_name, autodock_gpf)
        subprocess.call(['/usr/local/x86_64Linux2/autogrid4', '-p', autogrid_config, '-l', "%s.autogrid.log" % receptor_name], cwd=folder_path)

    with open('{}/{}'.format(folder_path, ligand_path), 'w') as ligand_file:
        ligand_file.write(ligand_value)

    if 'pdbqt' not in ligand_path:
        subprocess.call(['prepare_ligand4.py', '-F', '-l', ligand_path], cwd=folder_path)

    command = [
        '/opt/autodock-gpu',
        '--ffile', '{}.maps.fld'.format(receptor_name),
        '--lfile', '{}.pdbqt'.format(ligand_name),
        '--nrun', '30',
        '--resnam', '{}-{}'.format(receptor_name, ligand_name)
    ]

    returncode, gpu_output, gpu_error = run_process(command, folder_path, False)
    if returncode != 0:
        error = gpu_output if gpu_output != '' else gpu_error
        response = {
            'error': error
        }
        return jsonify(response)

    convert_command = [
        'cat', '{}-{}.dlg'.format(receptor_name, ligand_name),
        '|', 'grep', '"^DOCKED: "',
        '|', 'cut', '-b', '9-', '>', '{}-{}.pdbqt'.format(receptor_name, ligand_name)
    ]

    returncode, grep_output, grep_error = run_process(convert_command, folder_path, True)
    if returncode != 0:
        error = grep_output if grep_output != '' else grep_error
        response = {
            'error': error
        }
        return jsonify(response)

    with open('{}/{}-{}.pdbqt'.format(folder_path, receptor_name, ligand_name), 'r') as result_file:
        result_content = result_file.read()
    
    response = {
        'poses': result_content
    }

    if debug_mode:
        response['debug_info'] = {
            'gpu_output': gpu_output,
            'gpu_error': gpu_error,
            'grep_output': grep_output,
            'grep_error': grep_error
        }
 
    return jsonify(response)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True)