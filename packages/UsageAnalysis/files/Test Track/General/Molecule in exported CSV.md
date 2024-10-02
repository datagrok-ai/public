### Molecule column in exported CSV
1. Open the Dataset: Open the SPGI dataset.
2. Select Columns: select a few columns, ensuring that the 'Structure' column is included.
3. Export as CSV:
* Press the arrow icon next to the 'Save' button and choose "As CSV (options)...". 
* The Save as CSV form opens.
* Select the options: "Molecules as SMILES" and "Selected columns only".
* Export the dataset as a CSV file. The *.csv file should be downloaded.
4. Verify CSV File:
* Go to File > 'Close all' to close the current dataset.
* Open the downloaded *.csv file.
* Check that the molecule column is present and correctly visualized as SMILES.
* **Expected Results**: All selected columns are exported, and the 'Structure' column is correctly exported as SMILES.
---
{
  "order": 5,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}