#name: mergepredictions
#description: Merge linear regression and clustering data
#language: python
#input: dataframe df_clustering
#input: dataframe df_regression
#output: dataframe df_combined { viewer: scatterPlot(x:"actual", y:"predicted", color:"cluster")}

import pandas as pd

df_combined = df_regression.copy()
df_combined["cluster"] = df_clustering["cluster"]