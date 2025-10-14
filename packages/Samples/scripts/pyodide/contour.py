#name: Scottish hills contour plot
#description: Displays contour plot
#language: pyodide
#tags: demo
#sample: geo/scottish_hills.csv
#input: dataframe df
#output: graphics contour

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Create a grid of coordinates
x = np.linspace(df['Longitude'].min(), df['Longitude'].max(), 100)
y = np.linspace(df['Latitude'].min(), df['Latitude'].max(), 100)
X, Y = np.meshgrid(x, y)

# Interpolate height values on the grid
Z = griddata((df['Longitude'], df['Latitude']), df['Height'], (X, Y), method='cubic')
plt.figure(figsize=(12, 8))

# Create the contour plot with enhanced styling
contour = plt.contourf(X, Y, Z, levels=15, cmap='viridis', alpha=0.75)
contour_lines = plt.contour(X, Y, Z, levels=15, colors='black', linewidths=0.5)
plt.clabel(contour_lines, inline=True, fontsize=8, fmt='%1.0f')

# Add a color bar
cbar = plt.colorbar(contour)
cbar.set_label('Height (m)', rotation=270, labelpad=15)

# Add titles and labels with enhanced styling
plt.title('Contour Plot of Scottish Hills Heights', fontsize=16, fontweight='bold')
plt.xlabel('Longitude', fontsize=12)
plt.ylabel('Latitude', fontsize=12)

# Add grid and customize ticks
plt.grid(True, linestyle='--', alpha=0.5)
plt.xticks(fontsize=10)
