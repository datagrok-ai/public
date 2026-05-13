# Prompt examples

## Function chaining

| prompt                                                                            | interpretation |   |
|-----------------------------------------------------------------------------------|----------------|---|
| Adme properties for aspirin                                                       |                |   |
| Gasteiger for CHEMBL1234                                                          |                |   |
| Gasteiger for Apremilast                                                          |                |   |
| convert C[C@H](N)C(=O)N(C)CC(=O)N[C@@H](Cc1ccc(N)cc1)C(=O)O to helm               |                |   |
| convert GROKPEP-000002 sequence to molecule                                       |                |   |
| Convert to SMILES: GATTACA                                                        |                |   |
| Visualization of the aspirin synthesis path                                       |                |   |
| generate different variations of CHEMBL1235 and calculate molecular mass and logS |                |   |
| generate molecular structures of GROKPEP-000001 variants                          |                |   |
| ongoing studies in rats                                                           |                |   |
| optimize PSA of GROKPEP-000001                                                    |                |   |
|                                                                                   |                |   |
                                                              |                                                                                   |                |

## Database queries


| prompt                                        | interpretation |   |
|-----------------------------------------------|----------------|---|
| ADCs with ic50 above 400                      |                |   |
| ADCs with caspase activity > 8                |                |   |
| Which compounds are active against CHEMBL1827 |                |   |
| Shigella bioactivity                          |                |   |
| PK for LEVOFLOXACIN                           |                |   |
| CHO assay data                                |                |   |
| Northwind: sales by country                   |                |   |
| CHEMBL: targets involving thrombi             |                |   |


### DB identifiers

| prompt                  | interpretation |   |
|-------------------------|----------------|---|
| CHEMBL1234, CHEMBL13566 |                |   |
| GROKPEP001              |                |   |

## Context functions


| Context: prompt                                               | interpretation |   |
|---------------------------------------------------------------|----------------|---|
| MolCol: Show me the chemical space                            |                |   |
| MolCol: Calculate Lipinski descriptors                        |                |   |
| Scatter plot: zoom in, regression line, color by age          |                |   |
| Table View: add a histogram, close scatter plot               |                |   |
| Shell: show users view, switch to demog, show recent projects |                |   |


## Code generation

| prompt                                                                                                                                            | interpretation |                                                                                                                    |
|---------------------------------------------------------------------------------------------------------------------------------------------------|----------------|--------------------------------------------------------------------------------------------------------------------|
| JS: Add a view with the sketcher on the left, and the (predicted ADME properties                                                                  | calculated MW  | Gasteiger Partial Charges) on the right. As user sketches the molecule, the right part should update interactively |                |   |
| JS: given a dataframe “sales” parameter with “date” and “amount” columns, create a “cumsum” column with the cumulative sum of sales for that date |                |                                                                                                                    |
| For demog table, compute the mean and sd of HEIGHT group by RACE, SEX                                                                             |                |                                                                                                                    |
| Python: write a function that accepts a molecule and renders Gasteiger partial charges as a result                                                |                |                                                                                                                    |


## Visualizations

| prompt                                        | interpretation |   |
|-----------------------------------------------|----------------|---|
| Add a scatter plot, show height vs weight, show regression line, zoom in, color by age                      |                |   |
| Show me the distribution of activities per compound class              |                |   |
| Which compounds are active against CHEMBL1827 |                |   |
| Shigella bioactivity                          |                |   |
| PK for LEVOFLOXACIN                           |                |   |
| CHO assay data                                |                |   |