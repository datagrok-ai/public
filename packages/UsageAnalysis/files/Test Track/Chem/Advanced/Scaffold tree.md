Open linked dataset

Scaffold tree highlighting issue

1. Add a scaffold with coloring
2. Add a SF (Structure Filter) with a substructure that intersects with the scaffold
3. Disable the scaffold tree filter - its highlighting should not be displayed

   Note: Check this behavior with cloned view

***
1. Add a scaffold tree filter (by hand, load it, from the layout).
2. Set color to scaffolds.
3. Disable scaffold tree filter.
3. Close/open the entire panel.
3. Check the highlight in the grid.

***
1. Open the filter panel.
2. View > layout > clone.
2. On the original view, add a scaffold tree filter (enable color, do not check the checkbox!).
3. Disable the scaffold tree filter.
3. Close the filter panel.
3. Open the filter panel - the scaffold should not disappear from the panel.

***
1. Add a scaffold tree filter and SF.
2. Disable both.
3. Close the filter panel.
4. Open the filter panel.
5. Enable the scaffold tree
4. Enable SF - grid should not be empty.

***

#### Scaffold tree shouldn't reorient molecules in the grid

1. Open chembl-scaffolds.csv dataset from Demo > Chem > Chembl, click the **Smiles** column's header.
2. Go to the context panel > Chemistry > Rendering.
3. Select **Scaffold** in Scaffolds column dropdown - each molecule in dataset is reoriented by corresponding scaffold.
4. Open the Filter panel and add a scaffold tree filter.
3. Add benzene structure (c1ccccc1) to scaffold tree filter and color it.

- Expected result: benzene is highlighted with corresponding color but molecules and alignment is kept (in the grid).

#### Scaffold tree works after fixing an invalid structure

1. Open new table. Add scaffold tree.
2. Add new root structure
3. Start drawing some valid structure > table is filtered, drawn structure is highlighted
4. Edit structure to make it invalid > invalid parts are highlighted in editor
5. Fix the structure so that it is valid again
- Expected results: user can fix invalid structure and continue working with scaffold tree.

- Close all.

#### Applying Scaffold Tree Colors to Scatterplot
1. Load the smiles_50.csv dataset.Generate the Scaffold Tree viewer.
2. Apply Colors in Scaffold Tree: 
  - Select the first 3 molecules.
  - Assign a different color to each.
  - Ensure the checkboxes next to the selected scaffolds are checked.
3. Add Scatterplot viewer: 
  - As Color set 'Scaffold Tree_1' column
  - Expected result: Scatterplot should inherit colors form Scaffold Tree.
4. Save you data as new project and new layout.
  - Reopen your saved project and layout. 
  - Expected result: Coloring in Scatterplot remains the same; it should not reset or change.
5. Change name for 'Scaffold Tree_1' column. (e.g., to 'MyScaffoldTree')
  - Column name changes. 
  - The Scaffold Tree viewer name updates accordingly.Scaffold Tree`s name should be updated. 
  - The Color setting in the Scatterplot viewer also updates to reflect the new column name.

---
{
"order": 2,
"datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
