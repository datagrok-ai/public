#name: Linear Regression
#description: Linear regression
#reference: https://en.wikipedia.org/wiki/Linear_regression
#language: r
#tags: demo
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column x {type:numerical} [X axis column name]
#input: column y {type:numerical} [Y axis column name]
#input: bool interceptZero = FALSE [Force fit intercept zero]
#output: double k [Slope]
#output: double b [Intercept]

# Compose columns into data frame with required names
df = data.frame(
  x = data[[x]],
  y = data[[y]]
)

# Regression
if (interceptZero == TRUE) {
  regress = lm(y ~ 0 + x, data=df)
  k = summary(regress)$coefficients[1]
  b = 0
} else {
  regress = lm(y ~ x, data=df)
  k = summary(regress)$coefficients[2]
  b = summary(regress)$coefficients[1]
}
