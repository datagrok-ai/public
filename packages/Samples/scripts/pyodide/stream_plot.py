#name: Scottish hills demo
#description: Displays simple scatter plot
#language: pyodide
#tags: demo
#sample: geo/scottish_hills.csv
#input: dataframe data
#output: graphics scatter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Extract necessary columns
latitude = data['Latitude'].values
longitude = data['Longitude'].values
height = data['Height'].values

# Create a grid for the streamplot
X, Y = np.meshgrid(np.linspace(min(longitude), max(longitude), 100),
                   np.linspace(min(latitude), max(latitude), 100))

# Create a vector field for the streamplot
U = np.sin(X)  # Example vector field component
V = np.cos(Y)  # Example vector field component

# Interpolate height data to match the grid shape
from scipy.interpolate import griddata
height_grid = griddata((longitude, latitude), height, (X, Y), method='linear')

# Plot the streamplot
plt.figure(figsize=(10, 7))
stream = plt.streamplot(X, Y, U, V, color=height_grid, linewidth=2, cmap='viridis')
plt.colorbar(stream.lines, label='Height')
plt.title('Streamplot of Scottish Hills')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()
