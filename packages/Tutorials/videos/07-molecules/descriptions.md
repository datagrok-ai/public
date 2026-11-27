# Molecules

Playlist: **Molecules**

For a chemist new to Datagrok. Everything shown runs in the browser (Chem's RDKit and OpenChemLib); nothing that
needs Python, a container or an outside service (Descriptors, Curate, Mutate, IUPAC names, 3D, Identifiers).

Plan, one idea per clip, in this order:

1. **View and edit molecules in a table**: smiles.csv opened from Browse, SMILES drawn as structures, a molecule's
   Properties in the context panel; double-click a structure, the sketcher, paste SMILES, OK. smiles. (Merged on
   2026-11-22 from "Molecules in a table" and "Draw or edit a molecule": each was under 40 s.)
3. **Find and highlight a substructure**: Chem > Search > Substructure Search keeps the molecules with a fragment,
   reset brings them back, then the column's Chemistry > Highlight colors another fragment in every molecule.
   smiles. (Merged on 2026-11-22 from "Find molecules that contain a fragment", 0:31, and "Highlight a fragment in
   every molecule", 0:37.)
4. **Find similar and diverse molecules, and map the chemical space**: the two search viewers, following the current
   molecule, then Chem > Analyze > Chemical Space, a click on its map moving the similarity search, colored by a
   property. smiles. (Merged on 2026-11-22 from "Find similar and diverse molecules", 0:33, and "Map the chemical
   space", 0:31.)
5. **Calculate properties**: Chem > Calculate > Chemical Properties, new columns, then a filter on MW from the
   Search box and a sort by MW (extended on 2026-11-22 from 0:29). smiles.
7. **Activity cliffs**: similar molecules with very different activity. activity_cliffs (Chem demo).
9. Not delivered yet: **R-group analysis** (filmed, but the trellis plot it adds draws the same pie in every cell,
   which would confuse; feature kept as 08-r-group-analysis.feature).
10. **The Chem menu**: a tour of its groups (Search, Calculate, Analyze, Transform, Admetica), one line each,
   nothing opened. smiles. Delivered as 09; feature 10-chem-menu-tour.feature.

## 01 - View and edit molecules in a table (1:00)

Title: View and edit molecules in a table in Datagrok

Description: Open a file of SMILES (here Files > Demo > chem > smiles.csv in Browse) and Datagrok draws the
molecules as structures in the grid. Click one, open the context panel with the toggle at the bottom right, and
Chemistry > Properties lists its weight, LogP, polar surface area, formula and more; click another and the panel
follows. Double-click a structure to open the sketcher, draw or paste a SMILES into the field at the top, and OK
puts the new molecule into the cell: here row 1 becomes aspirin, and the panel shows its formula, C9H8O4.

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/01-molecules-in-a-table.feature

## 02 - Find and highlight a substructure (1:11)

Title: Substructure search and highlighting in Datagrok

Description: Two ways to work with a fragment in a column of molecules. Chem > Search > Substructure Search opens
the sketcher: draw the fragment or type its SMILES (here a carboxylic acid, C(=O)O), press OK, and only the
molecules that contain it stay, 314 of 1000, with the match drawn in each; the fragment stays on the filter panel,
and the panel's reset icon brings every molecule back. To mark a fragment without hiding any rows, click the
column's header, open the context panel and use Chemistry > Highlight > Sketch: here pyridine rings are colored
wherever they occur, and all 1000 rows stay. Uses the smiles demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/02-substructure-search.feature

## 03 - Find similar and diverse molecules, and map the chemical space (1:05)

Title: Similarity search, diversity search and chemical space in Datagrok

Description: Chem > Search > Similarity Search adds a viewer with the current molecule first and the most similar
ones after it, each with its score (Tanimoto on Morgan fingerprints); make another row current and the search
follows, or click a card to make that molecule current. Chem > Search > Diversity Search adds a set of molecules
picked to be as different from each other as possible. Chem > Analyze > Chemical Space then plots every molecule so
that similar structures sit close together (UMAP, computed in your browser): point at a dot to see its structure,
click it to make it current (the similarity search follows), and color the map by any column. Uses the smiles demo
dataset.

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/03-similarity-and-diversity.feature

## 04 - Calculate properties for every molecule (0:57)

Title: Calculate molecular properties in Datagrok, then filter and sort

Description: Chem > Calculate > Chemical Properties adds molecular properties as new columns of the table. Tick
the Chemical Properties calculator and the properties you need (molecular weight, LogP, polar surface area, H-bond
donors and acceptors, rotatable bonds, and more) and press OK; they are computed in your browser for every
molecule. The new columns work like any other: in the Search box (Ctrl+F) type a condition such as MW > 500 and
press Enter to keep only those molecules (211 of 1000 here), and double-click a column header to sort by it. Uses
the smiles demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/04-calculate-properties.feature

## 05 - Find activity cliffs (0:42)

Title: Find activity cliffs in Datagrok

Description: Chem > Analyze > Activity Cliffs finds pairs of similar molecules whose activity differs a lot. Check
the molecule and activity columns it picked and press OK: the molecules are mapped by similarity, and a line joins each cliff
pair. Show only cliffs hides the rest, pointing at a line shows both molecules with their activities, and the
cliffs button at the top lists every pair in a table, with how far their activities differ (drag the border
above it to give it room). Uses Chem's activity_cliffs demo file (Files > Demo > chem).

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/05-activity-cliffs.feature

## 06 - What the Chem menu does (0:52)

Title: A tour of the Chem menu in Datagrok

Description: A walk through the Chem menu, which appears in a table view once a column holds molecules. Search
finds molecules by a fragment, by similarity or for diversity; Calculate adds properties, scores and alerts as new
columns; Analyze has the structure-activity tools such as Chemical Space, R-Groups Analysis, Activity Cliffs and
Scaffold Tree; Transform converts notations, cleans up structures and runs reactions; Admetica predicts ADMET
properties on the server. Other videos in this playlist show several of these items step by step. Uses the smiles
demo dataset.

Feature: public/packages/Tutorials/bdd/features/howto/07-molecules/06-chem-menu-tour.feature
