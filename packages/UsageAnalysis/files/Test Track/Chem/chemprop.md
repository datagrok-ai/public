1. Open the file (Chem -> mol1k.sdf).
2. ML > Models > Train Model...
3. Select "molecule" column as features and "pIC50_HIV_Integrase" for prediction. 
4. Run the training, when the model is built run prediction for the same dataset.
5. Make sure that column with predictions is nearly equal to "pIC50_HIV_Integrase" (use the scatterplot).


---
{
  "order": 13,
  "datasets": [
    "System:AppData/Chem/mol1K.sdf"
  ]
}
