import io
from io import StringIO
import zipfile
from flask import Blueprint, Flask, request, Response
import logging
import sys
import json
import os
import pandas as pd

from modeling import get_all_engines, get_engine_by_type
from chemprop import ChemProp
from settings import Settings

app = Flask('grok_compute')
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
app.logger.handlers.clear()
app.logger.addHandler(handler)
app.logger.setLevel(logging.DEBUG)
bp = Blueprint('grok_compute', __name__)

headers_text = {'Content-Type': 'text/plain'}
headers_app_json = {'Content-Type': 'application/json'}
headers_app_octet_stream = {'Content-Type': 'application/octet-stream'}

@bp.route('/modeling/engines', methods=['GET'])
def modeling_engines():
    engines = get_all_engines()
    return _make_response(json.dumps(engines), headers=headers_app_json)

@bp.route('/modeling/train', methods=['POST'])
def modeling_train():
    id = request.args.get('id', '', type=str)
    type = request.args.get('type', '', type=str)
    table_server_url = request.args.get('table_server_url', '', type=str)
    table_token = request.args.get('table_token', '', type=str)
    predict = request.args.get('predict', '', type=str)
    parameter_values = json.loads(request.data)
    engine = get_engine_by_type(type)
    try:
        model_blob, log, performance = engine.train(id, table_server_url, table_token, predict, parameter_values)
        buffer = io.BytesIO()
        archive = zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED)
        archive.writestr('blob.bin', model_blob)
        archive.writestr('log.txt', log)
        archive.writestr('performance.json', json.dumps(performance))
        archive.close()
        return _make_response(buffer.getvalue(), headers=headers_app_octet_stream), 201
    except Exception as e:
        return _make_response(str(e)), 400

@bp.route('/modeling/train_chemprop', methods=['POST'])
def modeling_train_chemprop():
    id = request.args.get('id', '', type=str)
    type = request.args.get('type', '', type=str)
    table_str = request.args.get('table', '', type=str)
    table = pd.read_csv(StringIO(table_str))
    predict = request.args.get('predict', '', type=str)
    parameter_values = json.loads(request.data)
    engine = get_engine_by_type(type)
    chemprop = ChemProp()
    try:
        model_blob, log = chemprop.train_impl(id, table, predict, parameter_values)
        performance = engine.estimate_performance_impl(id, model_blob, table, predict,
                                                     True) if model_blob == '' else chemprop.estimate_performance_impl(id,
                                                                                                                   model_blob,
                                                                                                                   table,
                                                                                                                   predict,
                                                                                                                   False)
        buffer = io.BytesIO()
        archive = zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED)
        archive.writestr('blob.bin', model_blob)
        archive.writestr('log.txt', log)
        archive.writestr('performance.json', json.dumps(performance))
        archive.close()
        return _make_response(buffer.getvalue(), headers=headers_app_octet_stream), 201
    except Exception as e:
        return _make_response(str(e)), 400


@bp.route('/modeling/predict', methods=['POST'])
def modeling_predict():
    id = request.args.get('id', '', type=str)
    type = request.args.get('type', '', type=str)
    table_server_url = request.args.get('table_server_url', '', type=str)
    table_token = request.args.get('table_token', '', type=str)
    engine = get_engine_by_type(type)
    try:
        prediction = engine.predict(id, request.data, table_server_url, table_token)
        # TODO Return column list
        # TODO Add converter Pandas dataframe -> ColumnList or DataFrame
        return _make_response(json.dumps({'outcome': prediction['pred_0'].tolist()}),
                              headers=headers_app_json), 201
    except Exception as e:
        return _make_response(str(e)), 400

@bp.route('/modeling/predict_chemprop', methods=['POST'])
def modeling_predict_chemprop():
    id = request.args.get('id', '', type=str)
    table_str = request.args.get('table', '', type=str)
    table = pd.read_csv(StringIO(table_str))
    chemprop = ChemProp()
    try: 
        prediction = chemprop.predict_impl(id, request.data, table)
        return _make_response(json.dumps({'outcome': prediction[prediction.columns[0]].tolist()}),
                              headers=headers_app_json), 201
    except Exception as e:
        return _make_response(str(e)), 400

@bp.route('/modeling/models/<id>/<command>', methods=['POST'])
def modeling_model(id, command):
    # commands: zip, status
    return _make_response(None)

@bp.route('/modeling/models', methods=['GET'])
def modeling_models():
    # All models status
    return _make_response(None)

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

app.register_blueprint(bp, url_prefix=Settings.application_root)

if __name__ == '__main__':
    if Settings.keys is not None:
        app.run(host=Settings.host, port=Settings.port, ssl_context=Settings.keys, threaded=True)
    else:
        app.run(host=Settings.host, port=Settings.port, threaded=True)
