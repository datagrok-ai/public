#name: vdjwrapper
#description: vdjtools-1.2.1.jar wrapper
#language: python
#tags: demo, files, panel, cli, metadata-extractor, vdjtools
#require: vdjtools-1.2.1.jar
#input: file file
#output: dataframe metadata

import pandas as pd
import tempfile
from subprocess import Popen, PIPE

working_dir = tempfile.mkdtemp()
params = ['java', '-Xmx20G', '-jar', os.path.join(utils_dir,'vdjtools-1.2.1.jar'),'CalcSpectratype',file, 'out/0']
process = Popen(params, cwd=working_dir, stdout=PIPE, stderr=PIPE)
process.communicate()
metadata = pd.read_csv(os.path.join(working_dir, 'out/0.spectratype.nt.wt'),sep = '\t')