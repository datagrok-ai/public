1. Open demog.csv
2. Go to **ML > Models > Train Model...**. A dialog opens 
1. Set:
   * **Predict** to SEX
   * **Features** to WEIGHT and HEIGHT
1. Select the **Impute missing** checkbox - **a dialog opens**
1. Click RUN - **check the result**
1. Unselect the **Impute missing** checkbox 
4. Select the **Ignore missing** checkbox
1. Select the **Predict Probability** checkbox - **check the result**
1. Save the model as "TestDemog"

Repeat training the model predicting **WEIGHT** by **HEIGHT** (numeric data by numeric)

---
{
  "order": 1
}