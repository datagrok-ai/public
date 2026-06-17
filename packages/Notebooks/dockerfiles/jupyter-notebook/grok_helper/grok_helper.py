#!/usr/bin/env python

import json
import sys

from flask import Flask, Response, request, Blueprint

from cache import *
from environments import *
from settings import Settings

Settings.init(os.environ.get('GROK_HELPER_CONFIGURATION', 'configuration.yaml'))

try:
    import nbformat
    from jupyter_contrib_nbextensions import __file__ as contrib_init
    from nbconvert import HTMLExporter, PDFExporter
    from nbconvert.preprocessors import ExecutePreprocessor
except ImportError as e:
    if Settings.notebook_enabled:
        raise e
    else:
        pass

app = Flask(__name__)
bp = Blueprint('helper', __name__)
headers_text = {'Content-Type': 'text/plain'}
headers_app_json = {'Content-Type': 'application/json'}
headers_app_octet_stream = {'Content-Type': 'application/octet-stream'}

cache = Cache(Settings.cache_path)
environments_list = Environments(Settings.environments_path)


@bp.route('/', methods=['OPTIONS'])
def options():
    return _make_response('')


@bp.route('/info', methods=['GET'])
def info():
    return _make_response(json.dumps({
        'name': 'GrokHelper server',
        'version': '1.0.0',
        'python': sys.version_info[0]
    }))


@bp.route('/notebooks/convert', methods=['POST'])
def notebooks_convert():
    format = request.args.get('format', 'html', type=str)
    execute = request.args.get('execute', 'false', type=str) == 'true'
    nb = nbformat.reads(request.data, 4)
    if execute:
        ep = ExecutePreprocessor(timeout=600)
        ep.preprocess(nb, {'metadata': {'path': '/tmp'}})
    extra_template_basedirs = [os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            'nbconvert_templates')]
    if format == 'html':
        exporter = HTMLExporter(
            extra_template_basedirs=extra_template_basedirs,
            template_name='hide-input-html'
        )
    elif format == 'pdf':
        exporter = PDFExporter(
            extra_template_basedirs=extra_template_basedirs,
            template_name='hide-input-pdf'
        )
    else:
        return _make_response('Unknown format: "' + format + '"'), 400
    body, resources = exporter.from_notebook_node(nb)
    return _make_response(body)


@bp.route('/cache/put', methods=['POST'])
def cache_put():
    ce = CacheEntity().from_request(request)
    if ce.id != '':
        return _make_response(cache.put(request.data, ce)), 201
    else:
        return _make_response('Empty "id" parameter'), 400


@bp.route('/cache/exists', methods=['GET'])
def cache_exists():
    return _make_response('true' if cache.exists(CacheEntity().from_request(request)) else 'false')


@bp.route('/cache/get', methods=['GET'])
def cache_get():
    ce = CacheEntity().from_request(request)
    if ce.id != '':
        # TODO Convert to bytearray?
        return _make_response(cache.get(ce), headers=headers_app_octet_stream), 201
    else:
        return _make_response('Empty "id" parameter'), 400


@bp.route('/cache/file_path', methods=['GET'])
def cache_file_path():
    ce = CacheEntity().from_request(request)
    if ce.id != '':
        return _make_response(cache.get_entity_path(ce)), 201
    else:
        return _make_response('Empty "id" parameter'), 400


@bp.route('/cache/dir_path', methods=['GET'])
def cache_dir_path():
    ce = CacheEntity().from_request(request)
    if ce.id != '':
        path = cache.get_entity_path(ce, ext='')
        if not os.path.exists(path):
            os.mkdir(path)
        return _make_response(path), 201
    else:
        return _make_response('Empty "id" parameter'), 400


@bp.route('/cache/remove', methods=['DELETE'])
def cache_remove():
    ce = CacheEntity().from_request(request)
    if cache.exists(ce):
        return _make_response(cache.remove(ce)), 201
    else:
        return _make_response('Entity does not exists'), 400


@bp.route('/cache/clear', methods=['DELETE'])
def cache_clear():
    before_date = request.args.get('before_date', '', type=int)
    cache.clear(before_date)
    _make_response('ok')


@bp.route('/environments', methods=['GET'])
def environments():
    return _make_response(json.dumps(environments_list.list()), headers=headers_app_json)


@bp.route('/environments/<id>/create', methods=['POST'])
def environments_create(id):
    try:
        language = request.args.get('language', 'python', type=str)
        log = environments_list.create(id, request.data.decode('utf-8'), language=language)
        print(log)
        if id not in environments_list.list():
            raise Exception(log)
        return _make_response('created', headers=headers_text), 201
    except Exception as e:
        return _make_response(str(e)), 400


@bp.route('/environments/<id>', methods=['DELETE'])
def environments_remove(id):
    try:
        print(environments_list.remove(id))
        return _make_response('removed', headers=headers_text), 201
    except Exception as e:
        return _make_response(str(e)), 400


@bp.route('/environments/<id>/exists', methods=['POST'])
def environments_exists(id):
    try:
        language = request.args.get('language', 'python', type=str)
        if environments_list.exists(id, request.data.decode('utf-8'), language=language):
            return _make_response('true', headers=headers_text), 200
        return _make_response('false', headers=headers_text), 404
    except Exception as e:
        return _make_response(str(e)), 400


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
