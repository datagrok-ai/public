from flask import Flask, request, jsonify
import subprocess
import os
import tempfile

app = Flask(__name__)

def prepare_autogrid_config(receptor_basename, x, y, z):
    config = "%s.gpf" % receptor_basename
    with open(config, "w") as config_file:
        config_file.write("%s %s.pdbqt\n" % (receptor_basename, receptor_basename))
        config_file.write("npts {} {} {}\n".format(x, y, z))
        config_file.write("%s.maps.fld\n" % receptor_basename)
        config_file.write("spacing 0.375\n")
        config_file.write("receptor_types A C Fe N NA OA SA\n")
        config_file.write("ligand_types A C F NA OA HD\n")
        config_file.write("gridcenter auto\n")
        config_file.write("smooth 0.5\n")
        config_file.write("%s.A.map\n" % receptor_basename)
        config_file.write("%s.C.map\n" % receptor_basename)
        config_file.write("%s.F.map\n" % receptor_basename)
        config_file.write("%s.NA.map\n" % receptor_basename)
        config_file.write("%s.OA.map\n" % receptor_basename)
        config_file.write("%s.HD.map\n" % receptor_basename)
        config_file.write("%s.e.map\n" % receptor_basename)
        config_file.write("%s.d.map\n" % receptor_basename)
        config_file.write("dielectric -0.1465\n")

    return config

@app.route('/dock', methods=['POST'])
def dock():
    if 'receptor' not in request.files or 'ligand' not in request.files:
        return jsonify({'error': 'Both receptor and ligand files are required'}), 400

    receptor_file = request.files['receptor']
    ligand_file = request.files['ligand']

    receptor_basename = 'receptor'
    ligand_basename = 'ligand'

    receptor_path = "%s.pdb" % receptor_basename
    ligand_path = "%s.pdb" % ligand_basename

    receptor_file.save(receptor_path)
    ligand_file.save(ligand_path)

    subprocess.call(['prepare_receptor4.py', '-r', receptor_path])
    subprocess.call(['prepare_ligand4.py', '-F', '-l', ligand_path])

    x = int(request.form.get('x', 100))
    y = int(request.form.get('y', 100))
    z = int(request.form.get('z', 100))

    autogrid_config = prepare_autogrid_config(receptor_basename, x, y, z)

    subprocess.call(['/usr/local/x86_64Linux2/autogrid4', '-p', autogrid_config, '-l', "%s.autogrid.log" % receptor_basename])
    
    subprocess.call(['/opt/autodock-gpu', '--ffile {}.maps.fld --lfile {}.pdbqt --nrun 30 --resnam {}-{}.dlg'.format(receptor_basename, ligand_basename, receptor_basename, ligand_basename)])
    return jsonify({'message': 'Preparation completed successfully'}), 200

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, threaded=True)