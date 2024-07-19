1. Open the file (Chem -> mol1k.sdf). 
2. Run Chemprop - select "molecule" column as features and "pIC50_HIV_Integrase" for prediction. 
3. Run the training, when the model is built run prediction for the same dataset.
4. Make sure that column with predictions is nearly equal to "pIC50_HIV_Integrase" (use the scatterplot).


---
{
  "order": 13,
  "datasets": [
    "System:AppData/Chem/mol1K.sdf"
  ]
}