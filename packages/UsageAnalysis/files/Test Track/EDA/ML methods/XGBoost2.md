### XGBoost 12

1. Open the DataFrame:
* Open the **cars.csv** file from the Demo files.
* Alternatively, click on the 'star' icon in TestTrack (tooltip: "Open test data") to load the cars.csv dataset.

2. Train the Partial least squares regression:
* Go to **Top Menu > ML > Models > Train Model…**.
* In the **Predict field**, select **“Price”**.
* In the **Features** field, select all columns excluding  “Price” and “Model”.
* In the Model Engine dropdown, select “Eda: XGBoost”. 
3. In the model result:
- use sliders to vary **Rate, Lambda, Alpha**
- use clickers to vary **Iterations, Max depth**


**Expected Results:**

XGBoost Result data should be clickable and interactive. XGBoost model should be trained successfully without errors.



---
{
"order": 5,
"datasets": ["System:DemoFiles/cars.csv"]
}