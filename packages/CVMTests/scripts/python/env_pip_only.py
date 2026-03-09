#name: PythonEnvPipOnly
#description: Environment with pip-only extra dependencies
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.11, {pip: [chardet, pyyaml]}]
#input: string text
#output: string result

import chardet
import yaml

# Verify chardet works
detected = chardet.detect(text.encode('utf-8'))
encoding = detected['encoding']

# Verify pyyaml works
data = yaml.safe_load('key: ' + text)

result = 'ok:' + data['key']
