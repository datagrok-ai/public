#name: PythonLibs
#description: displays a list of locally installed python packages
#tags: python, libs
#language: python
#output: dataframe x

x = !pip freeze
x = pd.DataFrame(x)
x.rename(columns={0:'lib'},inplace = True)
x[['lib','version']] = x['lib'].str.split('==',expand=True)
