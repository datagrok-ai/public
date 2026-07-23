#name: ColumnMean
#language: python
#input: dataframe df
#input: string colName
#output: double mean_out
import pandas as pd
mean_out = float(pd.to_numeric(df[colName], errors='coerce').mean())
