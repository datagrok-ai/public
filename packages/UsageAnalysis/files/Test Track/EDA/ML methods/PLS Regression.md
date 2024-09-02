### Partial least squares regression

1. Open the DataFrame:
* Open the **cars.csv** file from the Demo files.
* Alternatively, click on the 'star' icon in TestTrack (tooltip: "Open test data") to load the cars.csv dataset.
2. Train the Partial least squares regression:
* Go to Top Menu > ML > Models > Train Model….
* In the Predict field, select “Price”.
* In the Features field, select all columns except “Price” and “Model”.
* Apply checkbox "One-hot encode string features"
* In the Model Engine dropdown, select “Eda: PLS Regression”. 
3. Check viewers: Loadings **scatterplot** and Regression coefficients **bar chart**, Explained variances **barchart**, scores **scatterplot**.

**Expected Results:**

Partial least squares regression Result data should be clickable and interactive. Partial least squares regression model should be trained successfully without errors, and all selected features should be included in the model.



---
{
"order": 2,
"datasets": ["System:DemoFiles/cars.csv"]
}