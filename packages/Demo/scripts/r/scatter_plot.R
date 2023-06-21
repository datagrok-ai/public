#name: Scatter Plot R
#language: r
#tags: demo, viewers
#input: dataframe t
#input: column xColumnName
#input: column yColumnName
#input: column colorColumnName
#output: graphics

require(ggplot2)

# Compose input columns into data frame with required names
data <- data.frame(x=t[[xColumnName]], y=t[[yColumnName]], color=t[[colorColumnName]])

# Plots
plotScatter <- ggplot(data, aes(x, y, colour=color), xlab=x, ylab=y) +
  labs(x=xColumnName, y=yColumnName) +
  geom_point()
print(plotScatter)
