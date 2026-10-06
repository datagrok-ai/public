#name: PythonLibs
#description: displays a list of locally installed python packages
#tags: python, libs
#language: python
#output: dataframe x
#test: PythonLibs().columns.length == 2

# importlib.metadata reads the interpreter that is actually running this script. Shelling
# out to `pip` assumed it was on PATH, which stopped being true once the gateway moved to
# conda, and parsing `pip freeze` breaks on conda entries anyway: they are reported as
# "name @ file:///..." with no "==" to split on.
from importlib.metadata import distributions

x = pd.DataFrame(
    sorted((d.metadata['Name'], d.version) for d in distributions()),
    columns=['lib', 'version'])
