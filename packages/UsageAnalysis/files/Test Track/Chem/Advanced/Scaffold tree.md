Open linked dataset

#### Scaffold tree highlighting issue

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
  - As Color set 'Scaffold Tree_1_canonical_smiles_colors' column
  - Expected result: Scatterplot should inherit colors form Scaffold Tree.
4. Save you data as new project and new layout.
  - Reopen your saved project and layout. 
  - Expected result: Coloring in Scatterplot remains the same; it should not reset or change.
5. Change Scaffold`s name in Description > Title for 'Scaffold Tree_1'. (e.g., to 'MyScaffoldTree')
  - Column name prefix changes, the ending is always **_canonical_smiles_colors** . 
  - The Scaffold Tree viewer name updates accordingly. Scaffold Tree`s name should be updated. 
  - The Color setting in the Scatterplot viewer also updates to reflect the new column name.
6. Open filter panel, Add filter > Scaffold tree filter. Select column dialog opens. **Only** canonical_smiles column is available to choose.


#### Modifying scaffold structure

1. Open the file **SPGI.csv**.
2. Navigate to **Chem > Analyze > Scaffold Tree**.  
   **Expected:** The **Scaffold Tree viewer** should be added.
3. In the Scaffold Tree toolbar, click the **plus (+)** icon to add new root structure.
4. In the sketcher, draw the structure: **C1CCCC1**. Click **Add**.  
   **Expected:** The scaffold is added as a root node in the tree.
5. Filter the dataframe using the newly added scaffold.  
   **Expected:** The dataframe displays **185 rows**, corresponding to the **185 hits** for this scaffold.
6. Hover over the newly added scaffold node and click **Edit scaffold** from the floating toolbar.  
7. Replace the structure with **C1CCCCC1** and click **Save**.  
   **Expected:**  
   - The scaffold node is updated.  
   - The dataframe is now filtered to **332 rows**, reflecting the **332 hits** for the new structure.

#### Handling empty scaffold values

1. Hover over any existing scaffold node and click **Edit scaffold** from the floating toolbar.  
2. In the sketcher, **delete the structure**, then click **Save**.  
   **Expected:**  
   - The scaffold is updated to an empty structure.  
   - The dataframe is filtered to **0 rows**, indicating **0 hits**.
3. In the Scaffold Tree toolbar, click the **plus (+)** icon.
4. In the sketcher, do not draw anything and click **Add**.  
   **Expected:**  
   - A new scaffold node is added to the tree.  
   - The scaffold has **0 hits**.  
   - The dataframe remains filtered to **0 rows**.

#### Other tests

1. Check that Magic wand icon is disabled for Scaffold tree on non-molecular data (for example, demog dataset).
2. Check that additional color columns will be added only after first Scaffold tree coloring and be removed when deleting a tree or closing the viewer.
3. Check the limit of 500 structure categories for Scaffold tree generation (you can use any molecule data, just remove extra rows, but do not use test datasets from top panel > Tools > Dev > Open test dataset):
  - if we have less than 500 categories - we generate a tree;
  - is we have 500 and more categories - magic wand icon is disabled.


---
{
"order": 2,
"datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
