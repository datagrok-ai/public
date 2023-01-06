#name: Fitting DRC
#description: Fit dose-response curve with a 4-parameter model with no constrain
#reference: https://en.wikipedia.org/wiki/Dose–response_relationship
#language: r
#tags: demo, fitting, drc
#sample: acidiq.csv
#input: dataframe table [Input data table]
#input: column x {type:numerical} [X axis column name]
#input: column y {type:numerical} [Y axis column name]
#output: graphics

require(drc)

m.4para = drm(table[[y]] ~ table[[x]], data=table, fct=L.4())
fitPlot = plot(m.4para, xlab=x, ylab=y, lwd=2, type='none', col="skyblue3")
points(table[[x]], table[[y]], col="skyblue3", pch=16)
grid(10, 10, lwd=1)

print(fitPlot)
