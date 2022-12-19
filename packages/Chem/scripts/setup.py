#name: SetupScaffold
#description: Installation of scaffoldgraph module
#language: python

import itertools
import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

modules = !pip list
counter = itertools.count(0)
[(next(counter), module) for module in modules if 'ScaffoldGraph' in module]
if str(counter) != 'count(1)':
  install('scaffoldgraph')