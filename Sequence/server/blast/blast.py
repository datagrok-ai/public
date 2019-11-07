#!/usr/bin/env python

import json
import multiprocessing
from chem import *
from cache import *
from flask import Flask
from flask import Response
from flask import request
from utilites import *
from datagrok.utils import *

import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor
from jupyter_contrib_nbextensions import __file__ as contrib_init


app = Flask(__name__)
cache = None
num_cores = 4
headers_app_json = {'Content-Type': 'application/json'}
headers_app_octet_stream = {'Content-Type': 'application/octet-stream'}


@app.route('/align', methods=['POST'])
def chem_lipinski():
    return _make_response(json.dumps(chem.lipinski(json.loads(request.data))))


if __name__ == '__main__':
    app.run(host=host, port=port, threaded=True)
