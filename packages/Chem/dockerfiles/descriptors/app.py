from flask import Blueprint, Flask, request, Response
from chem import *
import logging
import sys
import json

from utils import parallelize

app = Flask('descriptors')
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
app.logger.handlers.clear()
app.logger.addHandler(handler)
app.logger.setLevel(logging.DEBUG)
bp = Blueprint('descriptors', __name__)

headers_text = {'Content-Type': 'text/plain'}
headers_app_json = {'Content-Type': 'application/json'}
headers_app_octet_stream = {'Content-Type': 'application/octet-stream'}
cpu_count = 4


@bp.route('/chem/descriptors/tree', methods=['GET'])
def chem_descriptors_tree():
    return _make_response(json.dumps(get_descriptors_tree()), headers=headers_app_json)


@bp.route('/chem/descriptors', methods=['POST'])
def chem_descriptors():
    app.logger.debug("Started descriptors calculation")
    data = json.loads(request.data)
    molecules = data[request.args.get('molKey', 'molecules')]
    result = parallelize(get_descriptors, [molecules, data[request.args.get('desc', 'descriptors')]], [molecules],
                         cpu_count)
    app.logger.debug("Finished descriptors calculation")
    return _make_response(json.dumps(result), headers=headers_app_json)


@bp.route('/chem/molecules_to_canonical', methods=['POST'])
def chem_molecules_to_canonical():
    molecules = json.loads(request.data)
    result = parallelize(molecules_to_canonical, [molecules], [molecules], available_cores)
    return _make_response(json.dumps(result), headers=headers_app_json)


@bp.route('/chem/molecules_to_inchi', methods=['POST'])
def chem_molecules_to_inchi():
    molecules = json.loads(request.data)
    result = parallelize(molecules_to_inchi, [molecules], [molecules], available_cores)
    return _make_response(json.dumps(result), headers=headers_app_json)


@bp.route('/chem/molecules_to_inchi_key', methods=['POST'])
def chem_molecules_to_inchi_key():
    molecules = json.loads(request.data)
    result = parallelize(molecules_to_inchi_key, [molecules], [molecules], available_cores)
    return _make_response(json.dumps(result), headers=headers_app_json)


@bp.route('/chem/inchi_to_inchi_key', methods=['POST'])
def chem_inchi_to_inchi_key():
    inchi = json.loads(request.data)
    result = parallelize(inchi_to_inchi_key, [inchi], [inchi], available_cores)
    return _make_response(json.dumps(result), headers=headers_app_json)


@bp.route('/chem/inchi_to_smiles', methods=['POST'])
def chem_inchi_to_smiles():
    inchi = json.loads(request.data)
    result = parallelize(inchi_to_smiles, [inchi], [inchi], available_cores)
    return _make_response(json.dumps(result), headers=headers_app_json)


def _make_response(data, headers=None):
    response = Response(data)
    response.headers = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Headers': 'Content-Type',
        'Access-Control-Allow-Methods': 'POST,GET,DELETE,PUT,OPTIONS'
    }
    if headers is not None:
        response.headers.update(headers)
    return response


app.register_blueprint(bp)

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, threaded=True)
