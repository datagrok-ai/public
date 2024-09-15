### Softmax

1. Open the DataFrame:
* Open the **iris.csv** file from the Demo files.
* Alternatively, click on the 'star' icon in TestTrack (tooltip: "Open test data") to load the iris.csv dataset.
2. Train the Partial least squares regression:
* Go to **Top Menu > ML > Models > Train Model…**.
* In the **Predict field**, select **“Species”**.
* In the Features field, select all columns excluding “Species” and “col 1”.
* Apply checkbox "One-hot encode string features"
* In the Model Engine dropdown, select “Eda: Softmax”. 
3. In the model result vary **Hyperparameters**

**Expected Results:**

Softmax Result data should be clickable and interactive. Softmax model should be trained successfully without errors.



---
{
"order": 3,
"datasets": ["System:DemoFiles/iris.csv"]
}