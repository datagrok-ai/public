#name: STR Tool
#description: Calculates string weight
#language: python
#tags: demo, cli
#require: str-tool-macos, str-tool-linux
#input: dataframe data [Input data table]
#input: column input {type:categorical} [Impute data table columns]
#input: bool do_log = false [Use LOG weights]
#output: dataframe output {action:replace(data)} [Output data]
#output: string stdout [STDOUT content]
#output: string stderr [STDERR content]

import os
import tempfile
import pandas as pd
from subprocess import Popen, PIPE


working_dir = tempfile.mkdtemp()

input_path = os.path.join(working_dir, 'input.csv')
output_path = os.path.join(working_dir, 'output.csv')
input.to_csv(input_path)

params = [os.path.join(utils_dir, 'str-tool-linux'), '-i', input_path, '-o', output_path]
if do_log:
    params.append('-l')

process = Popen(params, cwd=working_dir, stdout=PIPE, stderr=PIPE)
stdout, stderr = process.communicate()

if process.returncode == 0:
    output = pd.read_csv(output_path)
