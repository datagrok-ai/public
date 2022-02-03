<!-- TITLE: Tests: Substructure search -->
<!-- SUBTITLE: -->

# Tests: substructure search

Substructure Search Performs substructure pattern search in list of SMILES.

## Testing scenarios

1. Open *smiles.csv* file
2. Click on *"smiles"* column

3. Open *"Substructure Search"* dialog in *"Active"* tab on [Property Panel](../../overview/navigation.md#properties)
   (or use context menu for *"smiles"* column, **Chem | Substructure Search**)

    * *"Substructure Search"* dialog is open

4. Click on "Pattern" field

    * "Molecule Structure Editor" is open.

5. In "Molecule Structure Editor" create molecule that matches SMILE: C1(=CC=CC=C1)

    * In "Pattern" field is displayed created molecule structure in editor

6. Execute *"Substructure Search"* dialog

    * Status bar with progress is displayed
    * New column has been added to table
    * Name of new column corresponds to molecule structure created in "Molecule Structure Editor" (
      Matches: C1(CCCCC1))
    * New column type is boolean
    * In new column all values are ```false``` (because SMARTs are searched by default in table and *"smiles"* table has
      SMILES)

7. Open *"Substructure Search"* dialog

8. Select last entry from history

    * Dialog fields are filled with last used values

9. Mark "Is Smarts" field as ```false```

10. Execute *"Substructure Search"* dialog

    * Status bar with progress is displayed
    * New column has been added to table
    * Name of new column corresponds to molecule structure created in "Molecule Structure Editor" (
      Matches: C1(=CC=CC=C1))
    * New column type is boolean
    * In new column values with ```true``` value meet the rows, in structure of molecules of which created structure was
      found

11. Test non-functional modules (UI, help, navigation, console)
