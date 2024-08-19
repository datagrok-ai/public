#name: Python Dataframe Graphics Cached
#description: df input/output
#language: python
#input: dataframe df
#output: dataframe resultDf
#output: graphics scatter
#meta.cache: client

import matplotlib.pyplot as plt
import random

resultDf = df
x = [random.randint(1, 9), random.randint(1, 9), random.randint(1, 9), random.randint(1, 9), random.randint(1, 9), 
     random.randint(1, 9), random.randint(1, 9), random.randint(1, 9)]
y = [10, 11, 12, 11, 14, 15, 16, 17]
plt.plot(x, y, 20, marker="x", color="red")
