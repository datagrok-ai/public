#!/bin/bash

# Upgrade pip
python -m pip install --upgrade pip

cd ../python-api
pip install .
cd ../docusaurus

# Set the path to the Python script and the argument
SCRIPT="../python-api/generate_docs.py"
SRC_DIR="../python-api/datagrok_api"
OUT_DIR="/api/py"

# Run the Python script with the argument
python "$SCRIPT" "$SRC_DIR" "$OUT_DIR"

# Copy the python-api.md file to /api/py as index.md
cp ../help/develop/python-api.md /api/py/index.md

# Copy the api-overview.md file to /api as index.md
cp ../help/develop/api-overview.md /api/index.md
