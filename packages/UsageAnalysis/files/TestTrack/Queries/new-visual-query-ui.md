---
feature: queries
realizes: [views.queries]
priority: p2
target_layer: manual-only
coverage_type: edge
manual_only_reason: |
  The Visual Query builder's column selectors and pointer-drag interactions
  must be driven by hand.
related_bugs: []
---

# New Visual Query — manual UI checks

This is the **manual companion** to `new-visual-query.md`. It covers the parts that the
autotest cannot exercise — the Visual Query builder's column selectors and
pointer-drag interactions must be driven by hand.

The manual steps below build a *new* parameterised visual query end-to-end and verify
the Debug tab, multi-section builder behaviour, and the post-save runtime.

## Pre-conditions

* You're logged in to a Datagrok server with `Postgres > NorthwindTest` (dev) /
  `Postgres > Northwind` (public) connection visible in **Browse > Databases**.

## Steps

1. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest**
   > **Schemas** > **public**. A list of tables opens.
2. Right-click the **customers** table and select **New Visual Query...** from the
   context menu — the Visual Query editor opens.
3. Set **Group by** to `companyname`.
4. Start setting **Order by** — verify that only `companyname` column can be selected
   (validation: Order by columns are restricted to what's already in Group by).
5. Go to **Data** and select `orders.shipcity` and `orders.shipcountry`.
6. Set **Where** to `customers.companyname contains C` and select the
   **Expose as function parameter** checkbox.
7. Set **Aggregate** to `sum(orders.freight)`.
8. Add **Group by**: `orders.shipcountry`.
9. Set **Pivot** to `orders.shipcity`.
10. Add **Order by**: `orders.shipcountry`.
11. Re-order the **Order by** fields (drag to swap) — verify the order in the list
    updates and is reflected in the resulting query.
12. Open the **Debug** tab and press the **bug icon** to **debug query** — no errors
    should appear.
13. **Toolbox > Actions > Run query...** — query result should open in a new view.
14. Save the query. Close All.
15. Run the saved query — the parameter dialog from step 6 (`Expose as function
    parameter`) should appear; check that the parameter is wired correctly.
16. On the **Toolbar**, change the parameter value and click **REFRESH** — the result
    should re-run with the new parameter.
17. Close all.

## What to look for

* Each builder section commits the column you intended (no off-by-one selections).
* Order by reorder produces the matching SQL ordering.
* Debug tab returns without errors and shows the generated query.
* After Save + reopen, the saved Visual Query opens in the Visual Query editor (not
  as a plain SQL query) with all builder fields preserved.
* The exposed parameter from the Where clause is what drives the runtime parameter
  input, and REFRESH re-runs against the new value.

## Cleanup

* Delete the visual query saved in step 14: right-click it under
  **Browse > Databases > Postgres > NorthwindTest** → **Delete...** → confirm.
* Close All.

---
{
  "order": 11,
  "datasets": []
}
