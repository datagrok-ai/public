1. Open demog.csv
2. Delete missing values (using **Select** | **Missing values…** on the menu ribbon)
3. On Toolbox, expand **Models** tab and click the **Train** button
4. Select **H2O** for **Model Engine** field
5. Enter `test_h2o_model` to the **Name** field
6. Select the **SEX** column for the **Predict** field
7. Select the **HEIGHT** and **WEIGHT** columns for the **Features** field
8. Click the **TRAIN** button

Repeat training the model predicting **WEIGHT** by **HEIGHT** (numeric data by numeric)
---
{
  "order": 1
}