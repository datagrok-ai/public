#name: Scatter Plot Python
#language: python
#tags: demo, viewers, hide-suggestions
#input: dataframe t
#input: column xColumnName {type:numerical}
#input: column yColumnName {type:numerical}
#input: column colorColumnName {type:numerical}
#output: graphics

import numpy as np
import matplotlib.pyplot as plt

color = t[colorColumnName].values
cmap = plt.cm.Spectral
norm = plt.Normalize(vmin=min(color), vmax=max(color))

plt.scatter(t[[xColumnName]], y=t[[yColumnName]], color=cmap(norm(color)), alpha=0.5)
plt.xlabel(xColumnName)
plt.ylabel(yColumnName)
plt.legend()
plt.show()
