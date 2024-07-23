Open linked datasets

Scaffold tree highlighting issue

1. Add a scaffold with coloring
2. Add a SF (Structure Filter) with a substructure that intersects with the scaffold
3. Disable the scaffold tree filter - its highlighting should not be displayed

   Note: Check this behavior with cloned view

***

1. Add a scaffold tree filter and SF.
2. Disable both.
3. Close the filter panel.
4. Open the filter panel.
5. Enable the scaffold tree
4. Enable SF - grid should not be empty.

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

[#2526](https://github.com/datagrok-ai/public/issues/2526) Scaffold tree shouldn't reorient molecules in the grid:

1. Open chembl-scaffolds.csv dataset, click the **Smiles** column's header.
2. Go to the context panel > Chemistry > Rendering.
3. Select **Scaffold** in Scaffolds column dropdown - each molecule in dataset is reoriented by corresponding scaffold.
4. Open the Filter panel and add a scaffold tree filter.
3. Add benzene structure to scaffold tree filter and color it.

Expected result: benzene is highlighted with corresponding color but molecules and alignment is kept (in the grid).

***

[#2626](https://github.com/datagrok-ai/public/issues/2626):

1. Open linked datasets.
2. Add a scaffold tree filter.
3. Draw some structure containing H atoms ([H]N).
4. Add a child node.

***

[#2139](https://github.com/datagrok-ai/public/issues/2139):

1. Open linked datasets.
2. Add a scaffold tree filter.
3. Add a new root structure.
4. Draw an invalid structure.
5. Fix the structure - check: user can fix invalid structure and continue working with scaffold tree.

***

[#2504](https://github.com/datagrok-ai/public/issues/2504) - Question mark on the Filter Panel - check:

1. Add scaffold tree filter.
3. Add scaffold tree viewer.
3. On the Filter Panel hover over the question mark.

---
{
"datasets": [
"System:DemoFiles/chem/SPGI.csv"
]
}
