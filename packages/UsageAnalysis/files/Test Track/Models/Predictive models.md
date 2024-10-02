**Train**

1. Open demo file: Browse > Files > Demo > Sensors > accelerometer.csv
2. Top menu > ML > Models > Train Model... Predictive model form opens. Fill it according to sample:

* Table field: accelerometer
* for Predict set: accel_x
* for Features set: accel_y, accel_z, time_offset
* Model Engine: EDA: PLS
* Components: 3

3. The result of modeling should appear. Press "Train" button.
4. Save Predictive model with name "Accelerometer model PLS".
5. Go to Toolbox menu > Models. Press "Train" button.
6. Fill the form according to sample:

* Table field: accelerometer
* for Predict set: accel_x
* for Features set: accel_y, accel_z, time_offset
* Model Engine: EDA: Linnear Regression

7. The result of modeling should appear. Press "Train" button. 
8. Save Predictive model with name "Accelerometer model LR".
9. File > Close All.

**Apply**

10. Open demo file once again: Browse > Files > Demo > Sensors > accelerometer.csv
11. Top menu > ML > Models > Apply Model... Apply predictive model form opens. Choose "Accelerometer model PLS" model. Inputs in "Apply predictive model" form should be set correctly (accel_x - accel_x). Press OK button. The predictive model result should appear.
12. Toolbox > Models. Run prediction "Accelerometer model LR". Inputs in "Apply predictive model" form should be set correctly (accel_x - accel_x). The predictive model result should appear.
13. Go to Top menu > Tools > Dev > Open test dataset. Set 1000 rows, 10 cols and "random walk as demo table". Preaa OK button. 
14. Go to ML > Models > Apply model. Choose "Accelerometer model LR" for new table and run it. The predictive model result should appear.

**Delete** 

15. Delete "Accelerometer model LR" and "Accelerometer model PLS" models in Toolbox. 

---
{
  "order": 6
}