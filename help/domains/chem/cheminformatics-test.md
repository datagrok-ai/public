<!-- TITLE: Tests: Cheminformatics -->
<!-- SUBTITLE: -->

# Tests: To inchi (inchi key)

Performs conversion from smiles to InChi and InChi key formats

## Testing scenario: inchi

1. Open *"smiles"* table

1. Convert *smiles* column to InChi from its context menu - **Chem | To InChi**
   (or **Chem | To InChi** action in *"Actions"* tab on [Property Panel](../../overview/navigation.md#properties)) for *
   smiles* column)

* Process bar shows conversion process
* New column *InChI* added to table
* *InChI* column contains values in InChi respectively values in smiles for each row

1. Convert *smiles* column to InChi Key from its context menu -
   **Chem | To InChi Key** (or **Chem | To InChi Key** action in *"Actions"* tab
   on [Property Panel](../../overview/navigation.md#properties)) for *smiles* column)

* Process bar shows conversion process
* New column *InChI Key* added to table
* *InChI Key* column contains values in InChi respectively values in smiles for each row

# Tests: map identifiers

Retrieves chemical identifiers for the specified source, using UniChem database.

## Testing scenario: map

1. Open *"smiles"* table

1. Open "Map Identifiers" dialog from context menu of *smiles* column - **Chem | Map Identifiers**
   (or **Chem | Map Identifiers** action in *"Actions"* tab
   on [Property Panel](../../overview/navigation.md#properties))
   for *smiles* column)

* "Get Identifiers" dialog is open
* There are three fields in dialog: "Table", "Column", "Source"
* Until value for "Source" field not selected, it's highlighted with red and OK is not active for clicking

1. Select "actor" value for "Source" field and execute dialog

* Process bar shows process
* New column *actor* added to table
* *actor* column contains matching chem IDs for molecules in table
* If row is empty, then corresponding molecule is missing in source database

1. Repeat step 2 for all available sources

1. Open **Tools | Console**

* Console displays all ChemMapIdentifiers() functions performed in previous steps

1. Open "Map Identifiers" dialog history

* History contains all previous runs with their parameters

See also:

* [Cheminformatics](cheminformatics.md)
* [Chem info panels test](chem-info-panels-test.md)
