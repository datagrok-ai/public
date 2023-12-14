from flask import Flask, request, jsonify
import subprocess
import os
import tempfile
import json
import hashlib

app = Flask(__name__)

def prepare_autogrid_config(folder_path, receptor_basename, x, y, z):
    config = "{}/{}.gpf".format(folder_path, receptor_basename)
    with open(config, "w") as config_file:
        config_file.write("receptor {}.pdbqt\n".format(receptor_basename))
        config_file.write("npts {} {} {}\n".format(x, y, z))
        config_file.write("gridfld {}.maps.fld\n".format(receptor_basename))
        config_file.write("spacing 0.375\n")
        config_file.write("receptor_types A C Fe N NA OA SA\n")
        config_file.write("ligand_types A C F NA OA HD\n")
        config_file.write("gridcenter auto\n")
        config_file.write("smooth 0.5\n")
        config_file.write("map {}.A.map\n".format(receptor_basename))
        config_file.write("map {}.C.map\n".format(receptor_basename))
        config_file.write("map {}.F.map\n".format(receptor_basename))
        config_file.write("map {}.NA.map\n".format(receptor_basename))
        config_file.write("map {}.OA.map\n".format(receptor_basename))
        config_file.write("map {}.HD.map\n".format(receptor_basename))
        config_file.write("elecmap {}.e.map\n".format(receptor_basename))
        config_file.write("dsolvmap {}.d.map\n".format(receptor_basename))
        config_file.write("dielectric -0.1465\n")

    return config

def calculate_hash(data):
    return hashlib.sha256(data.encode()).hexdigest()

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

@app.route('/dock', methods=['POST'])
def dock():
    raw_data = request.data
    json_data = json.loads(raw_data)

    receptor_value = json_data.get('receptor', '')
    ligand_value = json_data.get('ligand', '')

    x = request.args.get('x', 100)
    y = request.args.get('y', 100)
    z = request.args.get('z', 100)
    debug_mode = request.args.get('debug', False)
    subprocess_outputs = {}

    folder_name = calculate_hash(receptor_value) + str(x) + str(y) + str(z)
    folder_path = os.path.join(os.getcwd(), folder_name)
    receptor_name = 'receptor'
    ligand_name = 'ligand'

    receptor_path = "%s.pdb" % receptor_name
    ligand_path = "%s.pdb" % ligand_name

    receptor_path_prep = "%s.pdbqt" % receptor_name
    ligand_path_prep = "%s.pdbqt" % ligand_name
    
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        
        with open('{}/{}'.format(folder_path, receptor_path), 'w') as receptor_file:
            receptor_file.write(receptor_value)
        
        subprocess.call(['prepare_receptor4.py', '-r', receptor_path], cwd=folder_path)
        
        autogrid_config = prepare_autogrid_config(folder_path, receptor_name, x, y, z)
        subprocess.call(['/usr/local/x86_64Linux2/autogrid4', '-p', autogrid_config, '-l', "%s.autogrid.log" % receptor_name], cwd=folder_path)

    with open('{}/{}'.format(folder_path, ligand_path), 'w') as ligand_file:
        ligand_file.write(ligand_value)
    
    subprocess.call(['prepare_ligand4.py', '-F', '-l', ligand_path], cwd=folder_path)

    command = [
        '/opt/autodock-gpu',
        '--ffile', '{}.maps.fld'.format(receptor_name),
        '--lfile', '{}.pdbqt'.format(ligand_name),
        '--nrun', '30',
        '--resnam', '{}-{}'.format(receptor_name, ligand_name)
    ]

    process = subprocess.Popen(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, cwd=folder_path)
    gpu_stdout, gpu_stderr = process.communicate()
    subprocess_outputs['gpu_output'] = gpu_stdout
    subprocess_outputs['gpu_error'] = gpu_stderr

    convert_process = subprocess.Popen('cat {}-{}.dlg | grep "^DOCKED: " | cut -b 9- > {}-{}.pdbqt'.format(receptor_name, ligand_name, receptor_name, ligand_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=folder_path)
    grep_stdout, grep_stderr = convert_process.communicate()
    subprocess_outputs['grep_output'] = grep_stdout
    subprocess_outputs['grep_error'] = grep_stderr

    with open('{}/{}-{}.pdbqt'.format(folder_path, receptor_name, ligand_name), 'r') as result_file:
        result_content = result_file.read()
    
    response = {
        'poses': result_content
    }

    if debug_mode is False and result_content == '':
        response = {
            'error': next((value for value in list(subprocess_outputs.values()) if value), 'not found')
        }

    if debug_mode is True:
        response['debug_info'] = subprocess_outputs
 
    return jsonify(response)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True)