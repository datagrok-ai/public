#name: Financial Chart
#description: Charting tool to create standard financial charts given a time series
#reference: http://www.quantmod.com/examples/intro/
#language: r
#tags: demo, viewers
#sample: TSLA.csv
#input: dataframe rates [Financial data table]
#input: column date {type:dateTime; format:MM/dd/yyyy} [Date time column]
#input: column open {type:numerical} [Open rates column]
#input: column high {type:numerical} [High rates column]
#input: column low {type:numerical} [Low rates column]
#input: column close {type:numerical} [Close rates column]
#input: column volume {type:numerical} [Volume column]
#input: column adjusted {type:numerical} [Adjusted rates column]
#output: graphics [Chart series plot with moving average subplot]

require(xts)
require(quantmod)

# Compose input columns into data frame with required names
df.date = as.Date(rates[[date]], "%m/%d/%Y")
df.rates = data.frame(
  Open = rates[[open]],
  High = rates[[high]],
  Low = rates[[low]],
  Close = rates[[close]],
  Volume = rates[[volume]],
  Adjusted = rates[[adjusted]])

# Time series vector
data = xts(df.rates, df.date)

# Plot
chartSeries(data, TA=c(addMACD(), addVo()))
