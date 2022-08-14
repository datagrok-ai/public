#name: Surface Plot
#description: Surface plot
#reference: https://stat.ethz.ch/R-manual/R-devel/library/lattice/html/cloud.html
#language: r
#tags: demo, viewers
#sample: TSLA.csv
#input: dataframe t [Input data table]
#input: column X {type:numerical} [X axis column name]
#input: column Y {type:numerical} [Y axis column name]
#input: column Z {type:numerical} [Z axis column name]
#output: graphics [Surface plot]

require(lattice)
require(stats)

# Compose input columns into data frame with required names
data <- data.frame(x=t[[X]], y=t[[Y]], z=t[[Z]])

# Perform Loess fit
data.loess = loess(z ~ x * y, data=data, degree=2)
data.fit = expand.grid(list(x=seq(min(data$x),max(data$x),length.out=50),
                            y=seq(min(data$y),max(data$y),length.out=50)))
z = predict(data.loess, newdata=data.fit)
data.fit$z.predict = as.numeric(z)

# PLot
plotSurf = wireframe(z.predict ~ x * y, data=data.fit,
                     xlab=X, ylab=Y, zlab=Z,
                     colorkey=TRUE, col="grey", drape=TRUE,
                     col.regions=colorRampPalette(c("red", "yellow"))(100))
print(plotSurf)
