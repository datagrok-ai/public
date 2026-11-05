# Build-time check of the service: ligand/receptor preparation, autogrid4, and that autodock-gpu
# and clinfo load (the build host has no OpenCL device, so docking itself ends with an error).
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, '/app')
import autodock

here = os.path.dirname(os.path.abspath(__file__))
work = tempfile.mkdtemp()
os.chdir(work)
with open(os.path.join(here, 'ligand.pdb')) as f:
    ligand = f.read()
with open(os.path.join(here, 'rec.gpf')) as f:
    gpf = f.read()

clinfo = subprocess.Popen(['clinfo'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8', 'replace')
print(clinfo)
assert 'Number of platforms' in clinfo

client = autodock.app.test_client()
form = {'receptor': ligand, 'receptor_format': 'pdb', 'ligand': ligand, 'ligand_format': 'pdb',
        'autodock_gpf': gpf, 'pose_count': 1}
res = client.post('/autodock/dock_ligand', data=json.dumps(form)).get_json()
print(res)
assert 'error' in res and 'shared librar' not in res['error'], res

folder = os.path.join(work, autodock.calculate_hash(ligand + gpf))
for name in ['rec.pdbqt', 'ligand.pdbqt', 'rec.maps.fld', 'rec.C.map', 'rec.e.map']:
    assert os.path.getsize(os.path.join(folder, name)) > 0, name

if len(sys.argv) > 1:
    for name in os.listdir(folder):
        shutil.copy(os.path.join(folder, name), sys.argv[1])
shutil.rmtree(work)
print('selftest: OK')
