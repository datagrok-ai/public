#name: KNN Python
#description: Imputes (numerical) missing values using the kNN algorithm
#reference: https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm
#language: python
#tags: demo, hide-suggestions
#sample: demog.csv
#input: dataframe data [Input data table with NA elements]
#input: column_list imputeColumns {type:numerical} [Impute data table columns]
#input: column_list dataColumns [Input data table columns]
#input: int neighbours = 5 [Number of Nearest Neighbours used]
#output: dataframe data_out {action:replace(data)} [Output data table without NA elements]

import numpy as np
import pandas as pd
from fancyimpute import KNN

# Convert categories into numbers
for col in data.columns:
    if (data[col].dtype == 'object'):
        data[col]= data[col].astype('category')
        data[col] = data[col].cat.codes

# Convert to numpy array
columns = data.columns
data = data.as_matrix()

data_out = KNN(k=neighbours).fit_transform(data)

# Convert back to Pandas DataFrame
data_out = pd.DataFrame(data_out, columns=columns)
data_out = data_out[imputeColumns]
