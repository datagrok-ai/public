#name: PythonEnvPy312
#description: Environment with Python 3.12 template
#language: python
#environment: channels: [conda-forge], dependencies: [python=3.12, {pip: [httpx]}]
#input: string text
#output: string result

import sys
import httpx

# Verify correct Python version
major, minor = sys.version_info[:2]
assert major == 3 and minor == 12, f"Expected Python 3.12, got {major}.{minor}"

# Verify pip dependency works
result = f"py{major}.{minor}:httpx={httpx.__version__}"
