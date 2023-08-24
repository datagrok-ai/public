#name: R Graphics
#description: graphics output column input
#language: r
#input: dataframe df
#input: column xName {type:numerical; allowNulls:false} [Column x]
#input: column yName {type:numerical; allowNulls:false} [Column y]
#output: graphics scatter [scatter plot]

v1 <- df[[xName]]
v2 <- df[[yName]]
plot(v1, v2, main = "Main title",
     xlab = "X axis title", ylab = "Y axis title",
     pch = 19)