#name: ARIMA Forecasting
#description: Time series forecasting using ARIMA
#reference: https://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average
#language: r
#tags: demo, viewers
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column dates {type:dateTime; format:MM/dd/yyyy; allowNulls:false} [Date time column]
#input: column observations {type:numerical; allowNulls:false} [Observations column]
#input: int P = 4 [AR order (P)]
#input: int D = 3 [Degree of differencing (D)]
#input: int Q = 3 [MA order (Q)]
#input: int obsForecast = 25 [Number of observations to forecast]
#output: graphics [Forecast plot]

require(forecast)
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

# Fit an ARIMA model of order P, D, Q
fit <- arima(obs, order=c(P, Q, D))

# Predict next observations
predicted = forecast(fit, obsForecast)
plot(predicted, xlab="Time", ylab="Value")
grid(10, 10)
