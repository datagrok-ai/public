---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.rename-dependent-entity-reopen, projects.int.derive-then-save-inside-project]
realizes: [views.projects, data.menu.join-tables, views.space]
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - context-menu-rename-query
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
realized_as:
  - complex-derived-tables-spec.ts
  - complex-rename-spec.ts
related_bugs: [GROK-19212, GROK-21026]
---

# Complex — renaming and moving what a project depends on

A saved project keeps working when the things it is built from change:

1. The source tables of a pivot and a join are renamed (GROK-19212).
2. A query whose result is joined with another table is renamed
   (GROK-21026).
3. The project, its query and its script are moved into a Space.

Each part makes its own project. Saving tables from many sources,
adding tables, copies with and without Data sync, and sharing are
covered by `complex-integration.md`, `complex-augment.md`,
`complex-save-copy.md`, `complex-move.md` and the
`projects-lifecycle-*.md` files.

## Setup

1. Log in as the test user.
2. Names in this test: projects `complexRename`, `complexJoin` and
   `complexMove`; queries `complexQuery` and `complexMoveQuery`;
   script `complexScript`; Space `complexSpace`. Names have letters
   and digits only.
3. The Datagrok database is at **Browse > Databases > Postgres >
   Datagrok > Schemas > public**; on dev it is **Datagrok > Catalogs >
   datagrok > public**.

## Part 1 — renaming the sources of derived tables (GROK-19212)

1. **Open two files.**
   - Go to **Browse > Files > Demo > northwind**.
   - Double-click `orders.csv`.
   - **Verify:** the `orders` view opens with 830 rows.
   - Double-click `customers.csv`.
   - **Verify:** the `customers` view opens with 91 rows.

2. **Add a pivot table.**
   - Click the `orders` view tab.
   - In **Toolbox > Viewers**, click **Pivot table**.
   - Set **Group by** to `ShipCountry`.
   - Set **Aggregate** to `count(OrderID)`.
   - Click **ADD** in the top-right corner of the viewer.
   - **Verify:** the view *orders aggregation* opens.
   - Write down its row count from the status bar.

3. **Join the two tables.**
   - Open **Data > Join Tables...**.
   - In **Tables**, set the first table to `orders`.
   - Set the second table to `customers`.
   - Set both **Key Columns** to `CustomerID`.
   - Set **Join Type** to `inner`.
   - Click **OK**.
   - **Verify:** the join result opens as a new view.
   - Write down its row count.

4. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `complexRename` as the name.
   - Leave **Data sync** ON for every table.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.

5. **Rename the source tables.**
   - Right-click the `orders` view tab and choose **Table >
     Rename...**.
   - Change the name to `ordersR`.
   - Click **OK**.
   - **Verify:** the view tab reads `ordersR`.
   - Right-click the `customers` view tab and choose **Table >
     Rename...**.
   - Change the name to `customersR`.
   - Click **OK**.

6. **Save the renamed tables.**
   - Click **SAVE** on the ribbon.
   - Leave **Save original project** selected.
   - **Verify:** every table shows **Data sync** ON.
   - Click **OK**.

7. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `complexRename` into the search box.
   - Click the refresh icon.
   - Double-click the `complexRename` tile.
   - Known bug GROK-19212 (reopened): the pivot and the join do not
     come back; the task bar keeps showing *2 tasks* and no dialog
     appears.
   - **Verify:** four views open: `ordersR` with 830 rows, `customersR`
     with 91 rows, the pivot and the join with the row counts written
     down in steps 2 and 3.
   - **Verify:** no error balloon and no **Data loading error** dialog
     appear.
   - Click **SAVE** on the ribbon.
   - **Verify:** every table shows **Data sync** ON and a **CREATION
     SCRIPT**.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

## Part 2 — renaming a query used by a join (GROK-21026)

1. **Create the query.**
   - Go to the Datagrok database (see Setup).
   - Right-click `entity_types` and choose **New SQL Query...**.
   - Replace the **Name** with `complexQuery`.
   - Click **Save** on the ribbon.
   - Close the query editor.

2. **Open the query result and the table.**
   - Under **Datagrok**, double-click `complexQuery`.
   - **Verify:** the view `complexQuery` opens with rows.
   - Write down its row count.
   - In the Datagrok database, right-click `entity_types` and choose
     **Get All**.
   - **Verify:** the view `entity_types` opens with the same row count.

3. **Join them.**
   - Open **Data > Join Tables...**.
   - In **Tables**, set the first table to `entity_types`.
   - Set the second table to `complexQuery`.
   - Set both **Key Columns** to `id`.
   - Click **OK**.
   - **Verify:** the join result opens as a new view.

4. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `complexJoin` as the name.
   - Leave **Data sync** ON for every table.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

5. **Rename the query.**
   - Under **Datagrok**, right-click `complexQuery` and choose
     **Rename...**.
   - Change the name to `complexQueryRenamed`.
   - Click **OK**.

6. **Reopen.**
   - Go to **Browse > Dashboards**.
   - Type `complexJoin` into the search box.
   - Click the refresh icon.
   - Double-click the `complexJoin` tile.
   - Known bug GROK-21026: the join does not come back. The **Data
     loading error** dialog says *Could not resolve table
     "complexQueryRenamed"*; on dev no dialog appears and *Opening
     project* stays in the task bar.
   - **Verify:** three views open: the query result, `entity_types`
     and the join.
   - **Verify:** no **Data loading error** dialog appears.
   - Right-click the left sidebar and select **Close All**.

## Part 3 — moving the project and its sources into a Space

1. **Create the Space.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - Replace *New Space* with `complexSpace`.
   - Click **OK**.

2. **Create the query and the script.**
   - In the Datagrok database, right-click `entity_types` and choose
     **New SQL Query...**.
   - Replace the **Name** with `complexMoveQuery`.
   - Click **Save** on the ribbon.
   - Close the query editor.
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click **NEW** and choose **JavaScript Script...**.
   - Replace the template with the text below.
   - Click **SAVE** in the editor.

   ```
   //name: complexScript
   //language: javascript
   //output: dataframe df
   df = await grok.data.getDemoTable('demog.csv');
   ```

3. **Run both.**
   - Under **Datagrok**, double-click `complexMoveQuery`.
   - **Verify:** the view `complexMoveQuery` opens with rows.
   - In **Scripts**, right-click `complexScript` and choose **Run...**.
   - **Verify:** the `demog` view opens with 5,850 rows.

4. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter `complexMove` as the name.
   - Leave **Data sync** ON for both tables.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

5. **Move the project into the Space.**
   - Go to **Browse > Dashboards**.
   - Type `complexMove` into the search box.
   - Right-click the `complexMove` tile and choose **Move to
     Space...**.
   - Set **Space** to `complexSpace`.
   - Click **OK**.

6. **Move the query and the script into the Space.**
   - Drag `complexMoveQuery` from under **Datagrok** onto
     `complexSpace` in the Browse tree.
   - In the **Move entity** dialog, select **Move**.
   - Click **YES**.
   - Drag `complexScript` from **Scripts** onto `complexSpace`.
   - Select **Move**.
   - Click **YES**.
   - Click **Spaces > complexSpace**.
   - **Verify:** `complexMove`, `complexMoveQuery` and `complexScript`
     are listed.

7. **Open the moved project.**
   - In `complexSpace`, double-click `complexMove`.
   - **Verify:** the query result and `demog` open, `demog` with 5,850
     rows.
   - **Verify:** no error dialog appears.
   - Click **SAVE** on the ribbon.
   - **Verify:** both tables show **Data sync** ON and a **CREATION
     SCRIPT**.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

## Cleanup

- In **Browse > Dashboards**, right-click the `complexRename` tile and
  choose **Delete Project**. Click **DELETE** and wait until the dialog
  closes. Repeat for `complexJoin`.
- In **Spaces > complexSpace**, right-click `complexMove` and choose
  **Delete Project**. Click **DELETE** and wait until the dialog closes.
- In **Spaces > complexSpace**, right-click `complexMoveQuery` and
  choose **Delete**. Click **DELETE**.
- Right-click `complexScript` and choose **Delete**. Click **YES**.
- Under **Datagrok**, right-click `complexQueryRenamed` and choose
  **Delete**. Click **DELETE**.
- In **Browse > Spaces**, right-click `complexSpace` and choose
  **Delete Space**.
- **Verify:** the dialog says *Delete space "complexSpace"? This will
  delete space and its related data…*.
- Click **DELETE**.

## Expected results

- Renaming the source of a pivot or a join from its view tab keeps the
  derived tables in the project.
- Renaming a query keeps a join over its result working.
- A project whose query and script are moved into a Space, together
  with the project itself, opens with all its data.

## Automation notes

- In Part 2 the query name has no spaces on purpose: GROK-21026 shows
  only when the result table's name equals the query's internal name.
