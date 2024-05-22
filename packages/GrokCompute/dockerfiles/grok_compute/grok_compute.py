import io
import zipfile
from flask import Blueprint, Flask, request, Response
import logging
import sys
import json

from modeling import get_all_engines, get_engine_by_type


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


app.register_blueprint(bp)


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, threaded=True)
