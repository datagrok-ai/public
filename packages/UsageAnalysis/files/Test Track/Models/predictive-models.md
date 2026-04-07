#### 1. Train

1. Open **Browse > Files > Demo > Sensors > accelerometer.csv**
2. Go to **ML > Models > Train Model...**. A dialog opens 
1. Set:
   * **Features** to `accel_y`, `accel_z`, `time_offset`
   * **Model Engine** to `EDA: PLS`
   * **Components** to 3
   * **The result of modeling should appear**
3. Click "SAVE" to save the model as "Accelerometer_model_PLS"
5. Switch the **Model Engine** to EDA: Linnear Regression
6. Save the model as "Accelerometer_model_LR"

#### 2. Apply

1. Open **Browse > Files > Demo > Sensors > accelerometer.csv**
11. Go to  **ML > Models > Apply Model...** 
1. Select "Accelerometer_model_PLS" model - **Inputs in "Apply predictive model" form should be set correctly (`accel_y`, `accel_z`, `time_offset`)**
1. Click OK - **the predictive model result should appear in the last column**
12. Go to  **Top menu > ML > Models > Apply Model...** 
1.  Run prediction "Accelerometer model LR" - **check the inputs and the result**

#### 3. Apply on a new dataset

1. Go to **Tools > Dev > Open test dataset**
1. Set 1000 rows, 10 cols and "random walk as demo table". Click OK
14. Go to **ML > Models > Apply model**
1. Select "Accelerometer_model_LR" - **check the inputs and the result**

#### 4. Delete 

1. Go to **Browse > Platform > Predictive models**
5. Locate the "Accelerometer_model_LR" and "Accelerometer_model_PLS" models
3. Check the **Context Panel** tabs 
2. Delete the "Accelerometer_model_LR" and "Accelerometer_model_PLS" models

---
{
  "order": 6
}