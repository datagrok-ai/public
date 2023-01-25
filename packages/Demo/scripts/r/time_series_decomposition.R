#name: Time Series Decomposition
#description: Time series decomposition
#reference: https://en.wikipedia.org/wiki/Decomposition_of_time_series
#language: r
#tags: demo, viewers
#sample: births.csv
#input: dataframe data [Input data table]
#input: column dates {type:dateTime; format:MM/dd/yyyy; allowNulls:false} [Date time column]
#input: column observations {type:numerical; allowNulls:false} [Observations column]
#output: graphics [Components plot]

require(stats)
require(xts)

# Get columns data and correct data types
dates = as.Date(data[[dates]], "%m/%d/%Y")
observations = data[[observations]]

# Convert into timeseries vector
idxs <- order(dates)
dates = dates[idxs]
observations = observations[idxs]
frequencies = list(
  "daily" = 365,
  "weekly" = 52,
  "monthly" = 12,
  "quarterly" = 4,
  "yearly" = 1
)
obs <- ts(observations,
          start=c(as.numeric(format(dates[1],'%Y')),
                  as.numeric(format(dates[1],'%m'))),
          frequency=frequencies[[periodicity(dates)$scale]])

# Decompose observation components
components <- decompose(obs)

# Plot components
plot(components)
grid(10, 10)
