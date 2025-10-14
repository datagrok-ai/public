#name: Scottish hills ridge plot
#description: Displays ridge plot
#language: pyodide
#tags: demo
#sample: geo/scottish_hills.csv
#input: dataframe df
#output: graphics ridge

import numpy as np
import matplotlib.pyplot as plt

# Sort the data by height
df_sorted = df.sort_values(by='Height', ascending=False)
# Set up the matplotlib figure
plt.figure(figsize=(10, 6))

# Create a color palette
colors = plt.cm.viridis(np.linspace(0, 1, len(df_sorted['Hill Name'].unique())))

# Create the ridge plot
for i, (hill, subset) in enumerate(df_sorted.groupby('Hill Name')):
    heights = subset['Height']
    density, bins = np.histogram(heights, bins=30, density=True)
    bins_center = 0.5 * (bins[:-1] + bins[1:])
    plt.fill_between(bins_center, density + i, i, color=colors[i], alpha=0.6)

plt.title('Ridge Plot of Scottish Hills Heights')
plt.xlabel('Height (m)')
plt.ylabel('Density')
plt.show()
