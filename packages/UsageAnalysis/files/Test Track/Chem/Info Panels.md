1. Open linked datasets
2. Use smiles.cvs
3. Click the canonical_smiles column's header.
4. Go to the Context Pane and check all info panels (tabs).
  1. for testing panel Chemistry -> Rendering use chembl_scaffolds.cvs
  2. Choose 'Scaffold' as a Scaffold column and check 'Highlight scaffold'.
  3. Molecules in "Smiles" column should be aligned by scaffold and scaffold should be highlighted 

***

1. Open linked datasets (need to check smiles (smiles.csv), molV2000 (ml1K.sdf), molV3000 (ApprovedDrugs2015.sdf), smarts (SMARTS_example_temp) formats)
2. Make sure that structures of molecules are rendered.
3. Click a cell with a structure.
4. Check that all necessary panels are displayed on the **Context Panel.**
5. Expand each tab on the **Context Panel**.
6. Make sure the content for each info panel is displayed correctly.
---
{
  "order": 1,
  "datasets": [
    "System:DemoFiles/chem/smiles.csv",
    "System:AppData/Chem/chembl-scaffolds.csv",
    "System:AppData/Chem/mol1K.sdf",
    "System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf",
    "System:AppData/UsageAnalysis/test_datasets/SMARTS_example_temp.csv"
  ]
}
