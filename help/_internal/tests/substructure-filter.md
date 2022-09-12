<!-- TITLE: Tests: Substructure Filter -->
<!-- SUBTITLE: -->

# Tests: Substructure Filter

For filtering by columns that contain molecules Datagrock supports a special type of filter - Substructure filter.

## Adding Substructure Filter to the Filter Viewer

1. Open table with molecules (e.g. "smiles.csv")

1. Add Filters Viewer. Click on specical icon on Toolbox near "Search" section

* Filters is open. Filter for column with moluceles ("smiles" column) is displayed first in the list because it is special filter

1. Remove special substructure filter for "smiles" column by clicking "✕" icon for the corresponding column filter

* Special substructure filter for "smiles" column has been removed

1. Click on "hamburger" icon to open Filters menu

* Menu is open

1. Expand the first level menu "Add External"

* Second level submenu contains "Substructure Filter"

1. Click on "Substructure Filter" item from submenu

* Dialog for selecting columns for Substructure filter has opened

1. In column selection dialog, check the "smiles" column and click OK

* Special Subctructure filter has been added for the "smiles" column in Filter Viewer

1. Click on Sketcher field inside filter for "smiles" column

* Chemickal Sketcher is open

1. Draw molecule in the sketcher that corresponds to the structure "smiles: C1CCCCC1"

* Filtering happens as you draw structure (94 rows should be filtered for "smiles" table)

1. Click on OK button in Sketcher dialog

* Filtering is still applied. Sketcher field for the "smiles" column filter displays the drawn structure

1. Click on "Turn on\off" chekbox for "smiles" column filter (To turn the Substructure filter off)

* Substructure Filter appears to be inactive and is no longer applied to the table

1. Turn the filter back on by clicking on "Turn on\off" chekbox

* Filtering in table is applied. Filter is shown as active

1. Remove special substructure filter for "smiles" column by clicking "✕" icon for the corresponding column filter

* Special substructure filter for "smiles" column has been removed from Filters and filtering reset for table

## Dragging molecular column into Filters

Pre-requisites: Reproduce this case after last step of previous test. It is necessary to have open Filters and that there are no special substructure filters in there.

1. Drag-and-drop the column with molecules ("smiles" column) to Filters Viewer

* Special Subctructure filter has been added for the "smiles" column in Filter Viewer
* Filter is displayed as active

1. Repeat steps 8-13 from previous test (Adding Substructure Filter to the Filter Viewer). Expected results should be same and repeat exactly same for added Substructure filter using column drag-and-drop

## Using Subctructure Filtering in column "hamburger" menu

Pre-requisites: Open table with molecules (e.g. "smiles.csv")

1. Click on "hamburger" menu icon for column with molecules ("smiles" column) 

* Properties panel for the "smiles" column will open
* "Filter" tab will be expanded
* The "Filter" tab displays Sketcher for drawing a structure for filtering 

1. Draw molecule in the sketcher that corresponds to the structure "smiles: c1ccccc1"

* Filtering happens as you draw structure (924 rows should be filtered for "smiles" table)

1. Close "hamburger" menu for column with molecules ("smiles" column)

* Filtering is still applied for table

1. Open the Property Panel for molecular column ("smiles" column). Click on its title and it will become the current object

1. Expand "Filter" tab on PP for molecular column

* Filter set in the "hamburger" column menu is also displayed here

1. Press the "ESC" key on your keyboard

* Filters have been reset and are no longer applied to table


See also:

* [Cheminformatics](../../domains/chem/cheminformatics.md)
* [Filters](../../visualize/viewers/filters.md
