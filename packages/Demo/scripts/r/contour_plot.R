#name: Contour Plot
#description: Contour plot
#reference: https://www.image.ucar.edu/GSP/Software/Fields/Help/image.plot.html
#language: r
#tags: demo, viewers
#sample: volcano.csv
#input: dataframe data {columns:numerical} [Input data table]
#input: column_list columns {type:numerical} [Input data table columns]
#output: graphics [Contour plot]

library(fields)

# Compose input columns into data frame with required names
data <- data.matrix(data, rownames.force = NA)

# Plots
image.plot(data)
contour(data, add = TRUE)
