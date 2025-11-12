@echo off

:: Upgrade pip
python -m pip install --upgrade pip

cd ..\python-api
pip install .
cd ..\docusaurus
:: Set the path to the Python script and the argument
set SCRIPT=..\python-api\generate_docs.py
set SRC_DIR=..\python-api\datagrok_api
set OUT_DIR=api\py

:: Generate docs
python "%SCRIPT%" "%SRC_DIR%" "%OUT_DIR%"

:: Copy the python-api.md file to /api/py as index.md
copy ..\help\develop\python-api.md api\py\index.md

:: Copy the api-overview.md file to /api as index.md
copy ..\help\develop\api-overview.md api\index.md
