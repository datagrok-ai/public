---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.upload-save-reopen-golden, projects.int.derive-then-save-inside-project]
realizes: [views.projects, data.menu.link-tables, data.menu.aggregate-rows, file.menu.save.tables-as-project, viewers.pivot-viewer]
realized_as:
  - uploading-spec.ts
pyramid_layer: source-matrix
ui_coverage_responsibility:
  - save-project-dialog
  - link-tables-dialog
  - join-tables-dialog
  - pivot-table-add-to-workspace
  - aggregate-rows-add-to-workspace
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - selection-to-filter-linking-direction-asymmetry
  - spgi-v2-infinity-csv-source-ambiguity
  - project-naming-collision-risk-in-ci
  - add-to-workspace-semantics
  - original-step-renumbering-across-sync-on-sync-off-blocks
  - order-1-collision-with-upload-project-md
scope_reductions: []
related_bugs:
  - GROK-19103
  - GROK-18345
---

# Uploading — save projects built from different sources

Saves projects built from different data sources, once with **Data
sync ON** and once with **Data sync OFF**, and checks that each project
reopens with its tables, links and derived tables intact.

Each scenario is run once for every row of its table. A word in bold
in the steps, such as **Project**, means the value from that column.

## Setup

1. Log in as the test user.
2. **Create the Space `uploading`.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - Enter `uploading` as the name.
   - Click **OK**.
   - Go to **Browse > Files > Demo > northwind**.
   - Drag `customers.csv` onto the `uploading` Space in the Browse
     tree.
   - In the **Move entity** dialog, select **Copy**.
   - Click **YES**.
   - Drag `orders.csv` onto the `uploading` Space.
   - In the **Move entity** dialog, select **Copy**.
   - Click **YES**.

The tables used here share the customer key: `CustomerID` in the
files and `customerid` in the query **PostgresAll** (830 rows).
`customers.csv` has 91 rows, and its first two rows are `ALFKI` and
`ANATR`.

## Scenario 1: Two linked tables

| Project | Data sync | Creation script | Folder 1 | Item 1 | Folder 2 | Item 2 | Filtered |
|---|---|---|---|---|---|---|---|
| `TestCase1Sync` | ON | shown | Files > Demo > northwind | `customers.csv` | Files > Demo > northwind | `orders.csv` | 10 |
| `TestCase1NoSync` | OFF | not shown | Files > Demo > northwind | `customers.csv` | Files > Demo > northwind | `orders.csv` | 10 |
| `TestCase2Sync` | ON | shown | Databases > Postgres > NorthwindTest | `PostgresAll` | Databases > Postgres > NorthwindTest | `PostgresAll` | 11 |
| `TestCase2NoSync` | OFF | not shown | Databases > Postgres > NorthwindTest | `PostgresAll` | Databases > Postgres > NorthwindTest | `PostgresAll` | 11 |
| `TestCase3Sync` | ON | shown | Files > Demo > northwind | `customers.csv` | Databases > Postgres > NorthwindTest | `PostgresAll` | 10 |
| `TestCase3NoSync` | OFF | not shown | Files > Demo > northwind | `customers.csv` | Databases > Postgres > NorthwindTest | `PostgresAll` | 10 |
| `TestCase4Sync` | ON | shown | Spaces > uploading | `customers.csv` | Spaces > uploading | `orders.csv` | 10 |
| `TestCase4NoSync` | OFF | not shown | Spaces > uploading | `customers.csv` | Spaces > uploading | `orders.csv` | 10 |
| `TestCase5Sync` | ON | shown | Spaces > uploading | `customers.csv` | Files > Demo > northwind | `orders.csv` | 10 |
| `TestCase5NoSync` | OFF | not shown | Spaces > uploading | `customers.csv` | Files > Demo > northwind | `orders.csv` | 10 |
| `TestCase6Sync` | ON | shown | Spaces > uploading | `customers.csv` | Databases > Postgres > NorthwindTest | `PostgresAll` | 10 |
| `TestCase6NoSync` | OFF | not shown | Spaces > uploading | `customers.csv` | Databases > Postgres > NorthwindTest | `PostgresAll` | 10 |

1. **Open the first table.**
   - In **Browse**, go to **Folder 1**.
   - Double-click **Item 1**.
   - **Verify:** a table view opens with rows.

2. **Open the second table.**
   - At the bottom of the left panel, click the **Browse** tab.
   - In **Browse**, go to **Folder 2**.
   - Double-click **Item 2**.
   - **Verify:** a second table view opens with rows.

3. **Link the tables.**
   - Open **Data > Link Tables...**.
   - In **Tables**, set the first table to the one opened from
     **Item 1**.
   - Set the second table to the one opened from **Item 2**.
   - In **Key Columns**, set the key of the first table to its
     customer key.
   - Set the key of the second table to its customer key.
   - Set **Link Type** to **selection to filter**.
   - Click **LINK**.
   - Click **CLOSE**.

4. **Check the link.**
   - Click the first row of the first table.
   - Shift-click its second row.
   - Click the view tab of the second table.
   - **Verify:** the status bar shows **Filtered:** **Filtered**.

5. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter **Project** as the name.
   - Set **Data sync** for both tables to **Data sync**.
   - Click **OK**.
   - **Verify:** the balloon *Project "**Project**" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

6. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type **Project** into the search box.
   - Click the refresh icon.
   - Double-click the **Project** tile.
   - **Verify:** both tables open with rows.
   - **Verify:** no error balloon appears.
   - Click **SAVE** on the ribbon.
   - **Verify:** the **CREATION SCRIPT** block under each table is
     **Creation script**.
   - Click **CANCEL**.

7. **Check the link again.**
   - Click the first row of the first table.
   - Shift-click its second row.
   - Click the view tab of the second table.
   - **Verify:** the status bar shows **Filtered:** **Filtered**.

8. **Close.**
   - Right-click the left sidebar and select **Close All**.

## Scenario 2: File with a pivot table

| Project | Data sync | Creation script |
|---|---|---|
| `TestCase7Sync` | ON | shown |
| `TestCase7NoSync` | OFF | not shown |

1. **Open the table.**
   - Go to **Browse > Files > Demo > northwind**.
   - Double-click `orders.csv`.
   - **Verify:** the `orders` view opens with 830 rows.

2. **Add the pivot table.**
   - In **Toolbox > Viewers**, click **Pivot table**.
   - Set **Group by** to `ShipCountry`.
   - Set **Pivot** to `ShipVia`.
   - Set **Aggregate** to `count(OrderID)`.
   - Click **ADD** in the top-right corner of the viewer.
   - **Verify:** the view *orders aggregation* opens.

3. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter **Project** as the name.
   - Set **Data sync** for both tables to **Data sync**.
   - Click **OK**.
   - **Verify:** the balloon *Project "**Project**" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

4. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type **Project** into the search box.
   - Click the refresh icon.
   - Double-click the **Project** tile.
   - **Verify:** the views `orders` and *orders aggregation* open with
     rows.
   - **Verify:** no error balloon appears.
   - Click **SAVE** on the ribbon.
   - **Verify:** the **CREATION SCRIPT** block under each table is
     **Creation script**.
   - Click **CANCEL**.

5. **Close.**
   - Right-click the left sidebar and select **Close All**.

## Scenario 3: Database table with aggregated rows (GROK-19103)

| Project | Data sync | Creation script |
|---|---|---|
| `TestCase8Sync` | ON | shown |
| `TestCase8NoSync` | OFF | not shown |

1. **Open the table.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `orders` and choose **Get All**.
   - **Verify:** the `orders` view opens with 830 rows.

2. **Aggregate rows.**
   - Open **Data > Aggregate Rows...**.
   - **Verify:** the **Pivot table** panel opens.
   - Set **Group by** to `customerid`.
   - Set **Aggregate** to `count(orderid)`.
   - Click **ADD**.
   - **Verify:** a new view with one row per customer opens.

3. **Save.**
   - Click **SAVE** on the ribbon.
   - Enter **Project** as the name.
   - Set **Data sync** for both tables to **Data sync**.
   - Click **OK**.
   - **Verify:** the balloon *Project "**Project**" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

4. **Reopen.**
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type **Project** into the search box.
   - Click the refresh icon.
   - Double-click the **Project** tile.
   - **Verify:** both views open with rows.
   - Click **SAVE** on the ribbon.
   - **Verify:** the **CREATION SCRIPT** block under each table is
     **Creation script**.
   - Click **CANCEL**.
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** only the **Project** node is listed; the aggregated
     table has no separate project.

5. **Close.**
   - Right-click the left sidebar and select **Close All**.

## Cleanup

1. **Delete the projects.**
   - Go to **Browse > Dashboards**.
   - Type `TestCase` into the search box.
   - Right-click a `TestCase…` tile and choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Repeat for every `TestCase…` tile.

2. **Delete the Space.**
   - In **Browse > Spaces**, right-click `uploading` and choose **Delete
     Space**.
   - Click **DELETE**.

## Expected results

- Every case saves with both Data sync ON and OFF.
- After reopening, all tables of the case are there with data.
- Links between tables still filter after reopening.
- Pivot and aggregate results are saved inside the project.
