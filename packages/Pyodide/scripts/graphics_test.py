#name: Pyodide Graphics
#description: graphics output column input
#language: pyodide
#input: dataframe df
#input: column xName {type:numerical; allowNulls:false} [Column x]
#input: column yName {type:numerical; allowNulls:false} [Column y]
#output: graphics scatter [scatter plot]

import matplotlib.pyplot as plt

plt.plot(df[[xName]], df[[yName]])
plt.title('Test')
plt.xlabel(xName)
plt.ylabel(yName)
plt.show()
