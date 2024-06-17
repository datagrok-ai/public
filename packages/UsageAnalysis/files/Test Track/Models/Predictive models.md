You can train models using different methods and servers "H2o" or "OpenCPU". Previously created models can be applied with data similar in structure to training

1. Open "demog.csv" file.

2. Top menu > ML > Use "Missing Values Imtutation" for select and delete all nulls.

3. Top menu > Select > Use "Random" for select 50% random values. 
* Press "Extract Selected Ros". Extract 50% random values to new table. 
* Rename this table to "Train Data".
* 50% of the random values are extracted into a new table with name "Train Data".

4. Delete extracted rows from table "demog" (Edit > Remove > Selected rows) and rename this table to "Test Data".
Input table "demog" is divided into two tables "Train Data" and "Test Data".

5. Open "Train" from ML > Models menu.

6. Enter the following parameters:

* Name: Deep_Learning_H2o_Test
* Description: Test model with using method "Deep Learning" on H2o server
* Method: Deep Learning
* Table: "Train Data"
* Features: HEIGHT, WEIGHT (columns)
* Outcome: AGE (column)

7. The new model is trained. During training, the status bar is shown. When completed, a notification of the successful completion of the new model creation is shown.

8. Use all available methods for training new models on OpenCPU server. When creating new models, use the following rules:

* Name: {method_name}_OpenCPU_Test
* Description: Test model with using method {method_name} on OpenCPU server
* Remaining parameters is default

9. Open "Test Data" table. Apply all created models. 
* This can be done from different places, "Models" tab on toolbox, "Models" tab in table context panel or "Apply Model" dialog from "Train" from ML > Models menu.
* All created models should apply successfully and correctly.

---
{
  "order": 6
}