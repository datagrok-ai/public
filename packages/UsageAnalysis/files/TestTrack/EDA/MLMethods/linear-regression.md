---
feature: eda
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [eda.cp.train-predictive-model, eda.int.train-then-apply-linear-regression]
realizes: [eda.linear-regression, eda.softmax, eda.xgboost, eda.pls-regression]
realized_as:
  - linear-regression-spec.ts
related_bugs: []
---

### Linear Regression

1. Open the DataFrame:
* Open the **cars.csv** file from the Demo files.
* Alternatively, click on the 'star' icon in TestTrack (tooltip: "Open test data") to load the cars.csv dataset.
2. Train the Linear Regression Model:
* Go to Top Menu > ML > Models > Train Model….
* In the Predict field, select “Price”.
* In the Features field, select all columns except “Price” and “Model”.
* In the Model Engine dropdown, select “Eda: Linear Regression”. 

**Expected Results:**

Linear Regression modelResult data should be clickable and interactive.
The Linear Regression model should be trained successfully without errors, and all selected features should be included in the model.



---
{
"order": 1,
"datasets": ["System:DemoFiles/cars.csv"]
}