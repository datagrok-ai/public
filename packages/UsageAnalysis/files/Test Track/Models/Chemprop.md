1. Open smiles.csv.
2. On Toolbox, expand **Models** tab and click the **Train** button
3. Set Predict to Ring Count.
4. Set Features to canonical_smiles.
5. Click Train.
6. Save model with name "test_chemprop".
6. After training is complete, right-click on Views menu, then click **Close all**.
7. Navigate to **Browse** > **Platform** > **Dockers** and locate the container named **chem-chemprop**.
8. Right-click on the container or use `Ctrl + Click` to open the context menu. Select **Stop** from the menu to stop the container. Once the container has stopped, repeat the steps and choose **Run** to restart it.
9. Open smiles_only.csv. 
10. On Toolbox, expand **Models** tab and click the **Play** icon near **test_chemprop**
11. After the end of applying process, make sure that the column with the prediction is added to table and is has name Ring Count (2).
---
{
  "order": 4
}