---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.cp.upload-save-reopen-golden]
realizes: [views.projects, viewers.pivot-viewer, data.menu.aggregate-rows, data.menu.join-tables, views.space, views.databases]
realized_as:
  - complex-integration-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — one project from many sources

Opens tables from nine sources in one workspace, saves them together
as one project with Data sync, and checks that everything comes back
after a reopen. The sources are a file, a file from a Space, a database
table, a saved query and a script, plus a pivot, an aggregate, a join
and a clone built from them.

## Setup

1. Log in as the test user.
2. Names in this test: project `integration`, Space `integration`,
   script `integrationScript`.
3. **Create the Space.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - Enter `integration` as the name.
   - Click **OK**.
   - Go to **Browse > Files > Demo > northwind**.
   - Drag `customers.csv` onto the `integration` Space in the Browse
     tree.
   - In the **Move entity** dialog, select **Copy**.
   - Click **YES**.
4. **Create the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click **NEW** and choose **JavaScript Script...**.
   - Replace the template with the text below.
   - Click **SAVE** in the editor.

   ```
   //name: integrationScript
   //language: javascript
   //output: dataframe df
   df = await grok.data.getDemoTable('demog.csv');
   ```

## Scenario

1. **Open a file.**
   - Go to **Browse > Files > Demo > northwind**.
   - Double-click `orders.csv`.
   - **Verify:** the `orders` view opens with 830 rows.

2. **Open a file from the Space.**
   - Go to **Browse > Spaces > integration**.
   - Double-click `customers.csv`.
   - **Verify:** the `customers` view opens with 91 rows.

3. **Open a database table.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `products` and choose **Get All**.
   - **Verify:** the `products` view opens.

4. **Run a saved query.**
   - Go to **Browse > Databases > Postgres > NorthwindTest**.
   - Double-click the query **PostgresAll**.
   - **Verify:** the result view opens with 830 rows.

5. **Run the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Right-click `integrationScript` and choose **Run...**.
   - **Verify:** the `demog` view opens with 5,850 rows.

6. **Add a pivot table.**
   - Click the `orders` view tab.
   - In **Toolbox > Viewers**, click **Pivot table**.
   - In the viewer, set **Group by** to `ShipCountry`.
   - Set **Aggregate** to `count(OrderID)`.
   - Click **ADD** in the top-right corner of the viewer.
   - **Verify:** the view *orders aggregation* opens.

7. **Aggregate rows.**
   - Click the `demog` view tab.
   - Open **Data > Aggregate Rows...**.
   - **Verify:** the **Pivot table** panel opens.
   - Set **Group by** to `RACE`.
   - Set **Aggregate** to `avg(AGE)`.
   - Click **ADD**.
   - **Verify:** a new view with the aggregated table opens.

8. **Join two tables.**
   - Open **Data > Join Tables...**.
   - In **Tables**, set the first table to `orders`.
   - Set the second table to `customers`.
   - In **Key Columns**, set the key of `orders` to `CustomerID`.
   - Set the key of `customers` to `CustomerID`.
   - Set **Join Type** to `inner`.
   - Click **OK**.
   - **Verify:** the join result opens as a new view.

9. **Clone a table.**
   - Click the `products` view tab.
   - Click `products` at the left of the status bar.
   - In the **Context Panel**, expand **Actions**.
   - Click **Clone**.
   - **Verify:** the view `products (2)` opens.

10. **Check the workspace.**
    - On the left sidebar, click the **Dashboards** icon.
    - **Verify:** **New Dashboard** lists nine tables.
    - **Verify:** no other project node is listed.

11. **Save.**
    - Click **SAVE** on the ribbon.
    - **Verify:** each of the nine tables has its own **Data sync**
      toggle.
    - Enter `integration` as the name.
    - Leave **Data sync** ON for every table.
    - Click **OK**.
    - **Verify:** the balloon *Project "integration" uploaded.* appears.
    - In the **Share** dialog, click **CANCEL**.

12. **Reopen.**
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `integration` into the search box.
    - Click the refresh icon.
    - Click the `integration` tile.
    - In the **Context Panel**, expand **Content**.
    - **Verify:** **Content** lists nine tables.
    - Double-click the `integration` tile.
    - **Verify:** nine views open, each with rows.
    - **Verify:** no error balloon appears.

13. **Cleanup.**
    - Right-click the left sidebar and select **Close All**.
    - Right-click the `integration` tile and choose **Delete Project**.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - In **Browse > Spaces**, right-click `integration` and choose
      **Delete Space**.
    - Click **DELETE**.
    - In **Browse > Platform > Functions > Scripts**, right-click
      `integrationScript` and choose **Delete**.
    - Click **YES**.

## Expected results

- Tables from all sources can live in one project.
- The project saves with Data sync and reopens all tables, including
  the derived ones and the clone.
- Sources do not interfere with each other.
