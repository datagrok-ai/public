1. Open smiles. csv (Browse > Files > Demo Files/chem)
2. Double-click a molecule.
3. In the hamburger menu of the sketcher, click Favorites > Add to Favorites.
8. In the sketcher, enter **C1CCCCC1** to molecular input field.
9. In the hamburger menu, check the content in **Recent** and **Favorites**.
10. In the hamburger menu, click **Copy as SMILES**.
11. Go to the molecular input field and press CTRL+V.
12. In the hamburger menu, click **Copy as MOLBLOC**.
13. Go to the molecular input field and press CTRL+V.
14. Repeat steps 3-13 for all available sketcher types for selection from hamburger menu

Check:
* [#1608](https://github.com/datagrok-ai/public/issues/1608): Chem: Render a structure in a tooltip for the sketch box in the filters panel if the depicted substructure is large. Use InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H
*  [#2448](https://github.com/datagrok-ai/public/issues/2448): Some structures are displayed incorrectly when highlighted if they are in SMILES format:
   * Open SMILES_highlighted.csv dataset
   * Make these structures highlighted: either use scaffold tree to highlight or use structure filter
   
   Expected result: highlighted structure is displayed in the same way as non-highlighted (with stereochemistry kept)
---
{
  "order": 5
}