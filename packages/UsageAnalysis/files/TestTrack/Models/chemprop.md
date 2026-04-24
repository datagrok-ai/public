### Chemprop model

#### 1. Train

1. Open smiles.csv
2. Go to **ML > Models > Train model**
3. Set **Predict** to `Ring Count`
4. Set **Features** to `canonical_smiles` - **a model starts training, check the progress bar for correct info**
5. After the training is complete, check the interactive dashboard
5. Change values for **Activation**, **Split_type**, **Epochs** and click TRAIN
5. After the training is complete, save model with name "test_chemprop"
6. Change **Metric** to `auc` and click TRAIN - **A balloon should appear explaining why the model can’t be trained with the given value**

#### 2. Apply

1. Close All
9. Open smiles_only.csv. 
10. Go to **ML > Models > Apply model** and select **test_chemprop**
11. After the end of applying process, make sure that the column with the prediction is added to table and is has name Ring Count (2)

#### 3. Container

1. Go to **Browse** > **Platform** > **Dockers** and locate the container named **chem-chemprop**.
8. Right-click  the container and select **Stop**. Once the container has stopped, repeat the steps and select **Run** to restart it.

#### 4. Browse

1. Go to **Browse > Platform > Predictive models**
2. Locate the test_chemprop model - **check the Context Panel tabs**
3. Share the model
4. Delete the model
---
{
  "order": 4
}