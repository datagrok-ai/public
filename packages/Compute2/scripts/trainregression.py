#name: trainregression
#description: Train regression on wine quality data
#language: python
#input: dataframe df_wine
#input: double test_size = 0.2 {caption: Fraction of test data; min:0.1; max:0.5 }
#output: dataframe df_regression { viewer: scatterPlot(x:"actual", y:"predicted")}

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# --- Prepare features and target ---
X = df_wine.drop(columns=["quality", "wine_type"])  # features
y = df_wine["quality"]  # target

# --- Train/test split ---
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# --- Train model ---
model = LinearRegression()
model.fit(X_train, y_train)

# --- Predict and evaluate on the test data ---
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
print(f"R2 score: {r2:.2f}")

# --- Make dataframe actual vs predicted for full-zise data ---
df_regression = pd.DataFrame({
    'actual': y,
    'predicted': model.predict(X),
}, index=df_wine.index)
