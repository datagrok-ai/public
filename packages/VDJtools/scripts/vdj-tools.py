#name: vdjwrapper
#description: vdjtools-1.2.1.jar wrapper
#language: python
#tags: demo, files, panel, cli, metadata-extractor, vdjtools
#require: vdjtools.jar
#input: dataframe clonotypeDf
#input: string shcall
#input: string resultsfile
#output: dataframe metadata

from subprocess import Popen, PIPE
import pandas as pd
import tempfile
import shutil
import os

# create temporary folders to store files
mainFolder = tempfile.mkdtemp()
samplesFolder = tempfile.mkdtemp(suffix='_samples',dir=mainFolder)

# deposit file into the temporary directory
saveTo = os.path.join(samplesFolder, 'clonotypes.txt')
clonotypeDf.to_csv(saveTo, index=None, sep='\t', mode='a')

# create and execute sh call
shcall = shcall.replace('pathToJar', os.path.join(utils_dir, 'vdjtools.jar'))
shcall = shcall.replace('samplePath', saveTo)
shcall = shcall.split()
process = Popen(shcall, cwd=mainFolder, stdout=PIPE, stderr=PIPE)
process.communicate()

# retrive results from temporary storage
metadata = pd.read_table(os.path.join(mainFolder, resultsfile))

# delete temporary folder and all subfolders
shutil.rmtree(mainFolder, ignore_errors=True)