#name: Sample pandas output
#description: Produces dataframe with an integer sequence
#language: python
#input: int size
#output: dataframe table
#meta.queueName: python_docker

import pandas as pd
table = pd.DataFrame([i for i in range(size)], columns=['col'])