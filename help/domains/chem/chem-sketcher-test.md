<!-- TITLE: Tests: Chem Sketcher -->
<!-- SUBTITLE: -->

# Tests: chem sketcher

Sketcher is tool for drawing molecule structure or describing it in form of smiles, InCHI, InCHI key, etc.

## Testing scenario

1. Open *smiles.csv* file

1. Open sketcher from **Chem | Sketcher**

* "Sketcher" dialog is open
* There is field for entering smiles, InCHIs, etc.
* In dialog there is drawing editor for molecule

1. Click on "Filter" button on "Sketcher" dialog

* No rows fell under filter

1. Cancel filter by *Esc* key

1. Click on "Select" button on "Sketcher" dialog

* No rows is selected

1. Draw molecule in the editor that matches ```CC1=CC(=O)C=CC1=O``` in smiles format

1. Click on "Select" button on "Sketcher" dialog

* One row is selected in table that corresponds to molecule drawn in Sketcher

1. Click on "Filter" button on "Sketcher" dialog

* One row is filtered that corresponds to molecule drawn in Sketcher

1. Cancel filter and selection by *Esc* key

1. Enter ```[N+](=O)([O-])C1=CNC(=N)S1``` in sketcher text field

* Molecule was drawn in editor corresponding to entered smiles

1. Add molecule "Favorites" using "Hamburger" menu near sketcher text field

* Added molecule appeared in favorites list under "Hamburger" menu

1. Enter ```InChI=1S/C4H8N2O2/c1-3(5-7)4(2)6-8/h7-8H,1-2H3``` in sketcher text field

* Molecule was drawn in editor corresponding to entered InCHI

1. Repeat step 9 with Inchi Key, ChEMBL ID and etc.

1. Select arbitrary molecule in table and add its to "Favorites" (**Property Panel | Actions | Add to favorites**)

1. Open sketcher from **Chem | Sketcher** and expand "Hamburger" menu near sketcher text field

* Added to favorites molecule from table appeared in list

1. Select molecule witch added from table from the list under "Hamburger" menu

* Text field of Sketcher filled with molecule description in smiles format
* Selected molecule structure was drawn in editor

See also:

* [Cheminformatics Test](cheminformatics-test.md)
