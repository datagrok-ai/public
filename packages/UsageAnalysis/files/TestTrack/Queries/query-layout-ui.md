# Query Layout — manual UI checks

This is the **manual companion** to `query-layout.md`. It covers the parts that the
autotest `query-layout.test.ts` cannot exercise because
adding viewers to the Layout-tab preview requires either a drag-drop interaction
(canvas + Dart pointer-handlers — does not respond to synthesised events) or a JS
handle on the nested `TableView` instance inside the Layout tab, which is not exposed
through `grok.shell.tv` / `grok.shell.v`.

The autotest covers: creating a query, opening the Layout tab, running Play to
populate the preview, Save, and round-trip reopen so the Layout tab still attaches
its preview TableView. It does **not** add viewers, save layouts with viewers, save
projects, or test the post-reopen rehydration of layouts.

The manual steps below cover the full lifecycle of a saved query layout, project
save/reopen with that layout, and the Toolbox > File > Refresh interaction.

## Pre-conditions

* You're logged in to a Datagrok server with `Postgres > NorthwindTest` (dev) /
  `Postgres > Northwind` (public) connection.
* A query exists on the connection that returns multiple columns suitable for a few
  viewers (e.g. `PostgresAll` from the test track, or any `select * from products`-
  style query).

## Steps

1. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest**.
2. Right-click the **PostgresAll** query and select **Edit...**
3. Go to the **Layout** tab.
4. Click **Run query** — the preview TableView populates.
5. Add some viewers. Mix the layout:
   * Add at least two viewers side-by-side (no overlap).
   * Add at least one viewer **docked over another** — verify docking works as expected.
6. **Save** — the layout is persisted with the query.
7. **Close All**.
8. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest** and click the
   **PostgresAll** query — the **preview** should open with the new layout (the
   viewers and their docking arrangement from step 5).
9. Run the query — the **result** view should also open with the new layout.
10. Add some more viewers to the result view.
11. **Save the project** (top-right SAVE button → Save project).
12. **Close All**.
13. Open the saved project — check that the layout is restored: all viewers from
    steps 5 + 10, in the same arrangement.
14. Go to **Toolbox > File** and click **Refresh** — the layout should NOT change
    (Refresh re-runs the query but preserves the layout).

## What to look for

* Layout viewers persist correctly across Save → Close → reopen of the query
  (preview view shows the same docking configuration).
* Layout viewers persist across Save project → Close → reopen of the project.
* Refresh updates the data without resetting the layout (no viewers disappear or
  rearrange).
* No console errors / red balloons during any of the steps.
