#name: Tika
#description: Tika files metadata extractor
#language: python
#tags: demo, files, panel, cli, metadata-extractor, tika
#require: tika-extractor.jar
#input: file file
#output: map metadata
#condition: file.isFile && file.size < 1e6

import tempfile
from subprocess import Popen, PIPE

working_dir = tempfile.mkdtemp()
params = ['java', '-jar', os.path.join(utils_dir, 'tika-extractor.jar'), file, 'metadata.json']
process = Popen(params, cwd=working_dir, stdout=PIPE, stderr=PIPE)
process.communicate()
metadata = json.loads(open(os.path.join(working_dir, 'metadata.json'), 'r').read())
