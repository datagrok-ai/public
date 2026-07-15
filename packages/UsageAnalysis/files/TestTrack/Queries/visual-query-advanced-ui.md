---
feature: queries
target_layer: playwright
coverage_type: edge
priority: p2
realizes: []
realized_as: []
related_bugs: []
---

# Visual Query Advanced — manual UI checks

This is the **manual companion** to `visual-query-advanced.md`. The autotest
`visual-query-advanced.test.ts` runs the parameterised-query runtime path using a
pre-existing fixture (`postgres customers in @country`) and verifies the param
refresh + add result viewer + save project flow. This file covers everything that
DOM automation can't reach.

The autotest covers (via the fixture):
* Right-click → Run, parameter dialog with default value, OK, result view appears.
* Toolbar parameter input → type new value → REFRESH; row count changes.
* Add a Bar Chart viewer to the result view via Toolbox > Viewers icon.
* Click Save → Save-project dialog → set name → OK → dismiss auto-Share dialog.
* Server-side assertion: project persisted; in-memory both Grid + Bar chart present.

The manual steps below cover four carve-outs the autotest cannot exercise:

1. **Building a brand-new visual query through the visual builder** — Group by /
   Where (with Expose-as-parameter checkbox) / Aggregate / Pivot / Order by columns.
   The column selectors are canvas-rendered and the drag-drop uses Dart pointer-handlers
   that don't respond to JS-synthesised events. Same blocker as
   `new-visual-query-ui.md` step 5–11.
2. **Sharing the query** — the Share dialog and its share-target picker.
3. **Adding a layout (viewers + color coding + format + row size) inside the
   Layout tab** — same Layout-tab blocker as `query-layout-ui.md`.
4. **Re-opening a saved project from the Browse tree to verify the layout was
   restored** — SPA route doesn't resolve to the project's TableView in headless
   Playwright; works in a real human browser session.
5. **Delete rows / columns + Refresh with Enrich on** — verifying that Enrich
   restores deleted rows and columns without changing the layout.

## Pre-conditions

* You're logged in to a Datagrok server with `Postgres > NorthwindTest` (dev) /
  `Postgres > Northwind` (public) connection visible in **Browse > Databases**.
* You have a Northwind table with at least a few thousand rows so layout choices
  (color coding, row size, viewers) are visually meaningful.

## Steps

### Building the visual query (covers carve-out 1)

1. Right-click a Northwind table (e.g. `customers`) and select **New Visual Query...**
2. Configure the builder so the result has a parameterised `Where` clause:
   * **Where**: `customers.companyname contains C` and tick the
     **Expose as function parameter** checkbox.
   * **Group by**: `companyname`.
   * **Aggregate**: `sum(orders.freight)`.
   * **Pivot**: `orders.shipcity` (optional).
   * **Order by**: `companyname`.
3. Open the **Debug** tab and click the **bug icon** to debug the query — there
   should be no errors.
4. Click **Run query** in the Toolbox > Actions, OR press **Play** in the ribbon
   — the result view opens.
5. Set a custom name in the **Name** field (e.g. `test_visual_query`).
6. Click **Save** — the query is persisted.

### Sharing the query (covers carve-out 2)

7. Right-click the saved query in **Browse > Databases > Postgres > NorthwindTest**
   and select **Share...** (or use the SAVE dropdown's Share option).
8. In the Share dialog, pick a target group (e.g. `All Users` if you have
   permission), and confirm.
9. Verify the query now has the share marker in the Browse tree / Context Panel.

### Adding a post-process snippet (covers the part that
`query-postprocessing-ui.md` already documents)

10. Right-click the saved query → **Edit...**
11. Switch to the **Post-Process** tab and add a snippet on line 7, e.g.:
    ```javascript
    grok.shell.info(result.rowCount);
    ```
12. **Save**.

### Building a layout (covers carve-out 3)

13. Switch to the **Layout** tab and click **Run query** so the preview
    TableView populates.
14. Add some viewers (drag from Toolbox > Viewers): a Scatterplot, a Histogram,
    plus any other.
15. Apply **color coding** to a numeric column from the Grid context menu.
16. Change a column's **format** (e.g. number format) from the column header menu.
17. Change the grid **row size** (Grid > Settings > Row size).
18. **Save** the query.

### Verify post-process and saved layout fire on preview

19. **Close All**.
20. Click the saved query in the Browse tree to open its preview — verify:
    * The saved layout is restored (all viewers from step 14, color coding,
      format, row size).
    * A green info balloon with the row count appears (post-process ran).
21. Right-click the saved query → **Edit...**, then for each tab in the editor
    where Run is meaningful (Post-Process, Layout), click **Run query** and verify
    the green balloon fires.

### Edit + change layout cycle

22. Still in **Edit...**, switch to the **Query** tab and change a setting
    (e.g. add or remove an Order by column).
23. Switch to the **Layout** tab and add another viewer.
24. **Save**.
25. **Close All**.
26. Click the query to preview — verify the new layout (extra viewer) restored
    AND the post-process balloon still fires.

### Project save/reopen (covers carve-out 4 + 5)

27. Run the query (right-click → Run, OK with default param).
28. On the **Toolbar**, change the parameter value and click **REFRESH**.
29. Add some more viewers via Toolbox > Viewers.
30. Delete a few rows from the grid (right-click rows → Delete) and a column or
    two (right-click column header → Remove).
31. Toggle the **Enrich** option on (Toolbox > Source / Refresh dropdown — the
    exact location depends on the build) and click **Refresh** — verify:
    * The layout does NOT change (viewers stay).
    * Deleted rows and columns ARE restored.
32. Click **Save** → Save project, set a unique name, OK. Dismiss the Share
    dialog if it appears.
33. **Close All**.
34. Open the saved project from the Browse tree (Projects section) — verify the
    layout (all viewers, color coding, format, row size) is restored intact.

## What to look for

* The visual builder commits the column you intended in each section, no off-by-one
  selections.
* Save → Close → reopen preview shows the saved layout exactly.
* Post-process info balloon fires consistently from any Run path.
* Refresh with Enrich on restores deleted rows/columns without rearranging viewers.
* Saved project's Open opens the TableView with the full layout intact.
* No red error balloons or console errors during any step.
