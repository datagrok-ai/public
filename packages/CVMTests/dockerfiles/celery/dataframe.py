import pandas as pd

#name: cvmDataframe
#input: dataframe df
#output: dataframe result
def cvmDataframe(df):
    return df

#name: cvmDataframeNulls
#output: dataframe result
def cvmDataframeNulls():
    return pd.DataFrame({'f': [1.0, float('nan'), 3.0]})

#name: cvmIntInBound
#output: dataframe result
def cvmIntInBound():
    return pd.DataFrame.from_dict({'col1': [1000, 10000, 100000]})

#name: cvmIntOutBound
#output: dataframe result
def cvmIntOutBound():
    return pd.DataFrame.from_dict({'col1': [1000, 10000, 2 ** 31]})

#name: cvmEmptyDataframe
#output: dataframe result
def cvmEmptyDataframe():
    return pd.DataFrame({'col1': pd.Series([], dtype='float64')})
