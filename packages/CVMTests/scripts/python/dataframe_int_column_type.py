#name: Python Int Column
#language: python
#output: dataframe resultInBound
#output: dataframe resultOutBound
resultInBound = pd.DataFrame.from_dict({'col1': [1000, 10000, 100000]})
resultOutBound = pd.DataFrame.from_dict({'col1': [1000, 10000, 2**31]})
