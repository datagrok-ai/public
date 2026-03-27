1. Open demog.csv
2. Open the **Filter Panel**
3. **Hamburger menu** > **Add filter** > **Hierarchical Filter**
4. On the filter's header, click the hierarchy icon (tree-like icon with nodes). A column picker dialog opens showing all columns with checkboxes. Click "None" to uncheck all, then check `SEX`, `RACE`, `SEVERITY` — click OK. The tree displays hierarchy: SEX → RACE → SEVERITY.
5. Expand nodes: click the `>` arrow next to `F` → expand `Caucasian` → SEVERITY values (None, High, Low, Medium, Critical) become visible with counts.
6. Select a parent node (`Caucasian`): click its checkbox — all child SEVERITY nodes become checked, only Caucasian rows pass the filter, other RACE values show count 0. Verify row count matches Caucasian F rows.
7. Deselect child nodes: uncheck `Low` and `Medium` under SEVERITY — the parent `Caucasian` checkbox changes to a partial/indeterminate state (filled square instead of checkmark). Verify filtered row count decreases.
8. Click the hierarchy icon — column picker dialog opens. Click SEVERITY row, then click and drag it up above RACE. Release. Click OK. Verify the filter header shows `SEX / SEVERITY / RACE` (new order).
9. Collaborative filtering: in the hierarchical filter, click the text `None` (not the checkbox) under SEVERITY to select only that value — row count = 1815. Verify that categorical filters `SEX` and `SEVERITY` show non-zero counts only for the selected categories (`F` = 1815, `None` = 1815). Then in the `DIS_POP` categorical filter check only `AS` and `Indigestion` — verify row count = 339.
10. Save the layout via Toolbox → Layouts → SAVE.
11. Close the Filter Panel — all filters removed, all rows visible.
12. Apply the saved layout by clicking the layout thumbnail — verify the filter panel with hierarchical filter is restored and filtered row count matches the state before close.

---
{
"order": 7,
"datasets": ["System:DemoFiles/demog.csv"]
}
