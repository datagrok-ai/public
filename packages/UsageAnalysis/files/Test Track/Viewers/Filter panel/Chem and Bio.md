#### Refresh with chem filter

1. Open **SPGI** from the storage (not locally)
2. Click the `Core` column header
3. Go to **Context Pane > Chemistry > Rendering**
4. Set **Filter Type** to `Categorical`
5. Open the **Filter Panel** and verify:
    - The **Structure** filter uses the Sketcher
    - The **Core** filter contains categories (molecules)
6. In the **Structure** filter, draw the `c1cc2ccccc2cc1` molecule. 
   - Use different search options (e.g., *Included in*, *Exact*, etc.) and verify the filtering
7. Filter by the **Core** filter
8. Go to **Toolbox > File** and click **Refresh** — verify that filters are not cleared and data remains filtered
9. Go to **View > Layout > Clone View**
10. On the cloned view, in the **Filter Panel**, close the **Structure** filter tab — data should remain filtered on both views
11. Go to **Toolbox > File** and click **Refresh** — both views should remain filtered


#### ~ Columns behavior 

1. Run **Chem > Analyze > R-group analysis**
2. Open the **Filter Panel**  - verify that hidden special columns (name starts with "~") are not visible in the filter panel
3. Open the **Order Or Hide Columns** dialog - verify that the special columns are also not visible there

#### Bio

1. Open peptides.csv
1. Open the **Filter Panel**
1. For the **AlignedSequence** filter, enter `T-T-Y-K-N-Y-V` - verify the result
1. Set another separator value, e.g. `/` or `.`

---
{
"order": 2,
"datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/bio/peptides.csv"]
}
