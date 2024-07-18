1. Open:
   * Linked datasets (check for smiles, molV2000, molV3000 formats)
   * smiles_2_columns.csv (TODO: add to linked datasets)
1. Run Chem > Calculate > Descriptors 
2. Select arbitrary values in the dialog.
3. Click OK.
4. Make sure that the column with calculated values is added to the table.
Do the asme for each section in the Calculate menu ()
---
{
  "order": 1,
   "datasets": [
    "System:DemoFiles/chem/SPGI.csv",
    "System:DemoFiles/chem/smiles.csv",
    "System:AppData/Chem/mol1K.sdf",
    "System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf"
  ]
}
