<!-- TITLE: Tests: Predictive Models -->
<!-- SUBTITLE: -->

# Tests: Predictive Models

You can train [models](../plugins/predictive-modeling.md) using different methods and servers 
"[H2o](http://h2o.ai)" or "[OpenCPU](https://topepo.github.io/caret/index.html)". 
Previously created models can be applied with data similar in structure to training

## Testing scenario

1. Open "demog.csv" file

1. Use "[Select Missing Values](../dialogs/select-missing-values.md)" for select and delete all nulls

1. Use "[Select Random](../dialogs/select-random-rows.md)" for select 50% random values. 
Extract 50% random values to new table. Rename this table to "Train Data"
   * 50% of the random values ​​are extracted into a new table with name "Train Data"

1. Delete extracted rows from table "demog" and rename this table to "Test Data"
   * Input table "demog" is divided into two tables "Train Data" and "Test Data"

1. Open "Train" from **Tools | Predictive Model** menu
   
1. Enter the following parameters:

       Name: Deep_Lerning_H2o_Test
       Description: Test model with using method "Deep Learning" on H2o server
       Method: Deep Learning
       Table: "Train Data"
       Feauchers: HEIGHT, WEIGHT (columns)
       Outcome: AGE (column)

1. Leave the remaining parameters as default.

1. Click "Train" button
   * The new model is trained. During training, the status bar is shown. When completed, a notification of the 
     successful completion of the new model creation is shown

1. Use all available methods for training new models on H2o server. 
When creating new models, use the following rules:

       Name: {method_name}_H2o_Test
       Description: Test model with using method {method_name} on H2o server
       
       Remaining parameters is default

1. Use all available methods for training new models on OpenCPU server. When creating new models, use the following rules:
    
       Name: {method_name}_OpenCPU_Test
       Description: Test model with using method {method_name} on OpenCPU server

       Remaining parameters is default
 
1. Open "Test Data" table. Apply all created models. 
This can be done from different places, "Models" tab on toolbox, "Models" tab in table property panel or "Apply Model" dialog from **Tools | Predictive Modeling**.
   * All created models are applied successfully and correctly

See also: 
  * [Predictive modeling](../plugins/predictive-modeling.md)
  * [Predictive modeling browser test](../tests/predictive-models-browser-test.md)
  * [Predictive modeling tutorial](../tutorials/predictive-modeling.md)
  * [OpenCPU](https://www.opencpu.org/)
  * [H2O](http://h2o.ai/)
s