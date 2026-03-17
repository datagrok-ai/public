#name: PythonEnvReuse
#description: Tests that existing environment is correctly detected and reused
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.11, {pip: [chardet]}]
#input: string text
#output: string result

import chardet

detected = chardet.detect(text.encode('utf-8'))
result = 'ok:' + str(detected['confidence'])
