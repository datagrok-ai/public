---
feature: projects
ui_coverage_split_to:
  - complex-ui.md
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.rename-dependent-entity-reopen, projects.cp.share-with-unshared-deps, projects.cp.share-spaces-datasync, projects.cp.view-and-use-failure-state, projects.int.derive-then-save-inside-project]
realizes: [views.projects, data.menu.join-tables, views.space]
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - context-menu-rename-project
  - context-menu-rename-query
  - context-menu-rename-script
  - pcmdShareProject
  - share-dialog-permissions-editor
  - logout-login-as-second-user
  - data-sync-refresh-verification
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.projects.dragDropOntoDashboard
  - helpers.playwright.session.logoutAndLoginAs
  - helpers.playwright.projects.move
  - helpers.playwright.projects.rename
unresolved_ambiguities:
  - another-user-second-user-single-identity-vs-two
  - order-of-access-level-grants-in-step-12
  - table-from-a-script-semantics-in-step-1
  - pivot-table-add-and-aggregate-rows-add-ui-controls
  - move-to-any-file-share-then-to-any-space-in-step-10
  - rename-ui-access-points-in-step-9
  - verify-data-sync-updates-correctly-in-step-11
  - grok-19728-in-chain-rev-3-emission-rationale
scope_reductions: []
ui_companion: complex-ui.md
realized_as:
  - complex-derived-tables-spec.ts
  - complex-rename-spec.ts
  - complex-share-second-user-spec.ts
related_bugs: [GROK-19212, GROK-19103, GROK-19403, GROK-18345, GROK-19728, github-3550]
---

# Complex — the full project lifecycle

One long run through everything a project goes through:

1. Tables from many sources are saved together.
2. More tables are added.
3. Copies are saved with and without Data sync.
4. Tables and the entities behind them are renamed and moved.
5. The project is shared and opened by a second user.

The run follows the reproduction paths of six known bugs:

- **GROK-19212.** A renamed table loses its Data sync on reopen.
- **GROK-19103.** A join result is saved as a separate, broken project.
- **GROK-18345.** A shared project with a Space table does not open
  for the recipient.
- **GROK-19403.** A shared project does not open because its script is
  not shared.
- **GROK-19728.** A view-only user can edit a failing creation script.
- **github-3550.** Renaming a query breaks the project.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. All names in this test start with `complex`.
3. **Space.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - Enter `complex` as the name.
   - Click **OK**.
   - Go to **Browse > Files > Demo**.
   - Drag `demog.csv` onto the `complex` Space in the Browse tree.
   - In the **Move entity** dialog, select **Copy**.
   - Click **YES**.
4. **Query.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `orders` and choose **New SQL Query...**.
   - Replace the **Name** with `complexQuery`.
   - Click **Save** on the ribbon.
   - Close the query editor.
5. **Script.**
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

## Scenario

1. **Open a file.**
   - Go to **Browse > Files > Demo > northwind**.
   - Double-click `orders.csv`.

2. **Open a file from the Space.**
   - Go to **Browse > Spaces > complex**.
   - Double-click `demog.csv`.

3. **Open a database table.**
   - Go to **Browse > Databases > Postgres > NorthwindTest > Schemas >
     public**.
   - Right-click `customers` and choose **Get All**.

4. **Run the query.**
   - Under **NorthwindTest**, double-click `complexQuery`.
   - **Verify:** the result view opens with 830 rows.

5. **Run the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Right-click `complexScript` and choose **Run...**.
   - **Verify:** a second `demog` view opens with 5,850 rows.

6. **Add a pivot table.**
   - Click the `orders` view tab.
   - In **Toolbox > Viewers**, click **Pivot table**.
   - Set **Group by** to `ShipCountry`.
   - Set **Aggregate** to `count(OrderID)`.
   - Click **ADD** in the top-right corner of the viewer.

7. **Aggregate rows.**
   - Click the `customers` view tab.
   - Open **Data > Aggregate Rows...**.
   - Set **Group by** to `country`.
   - Set **Aggregate** to `count(customerid)`.
   - Click **ADD**.

8. **Join two tables (GROK-19103).**
   - Open **Data > Join Tables...**.
   - In **Tables**, set the first table to `customers`.
   - Set the second table to the `complexQuery` result.
   - Set both **Key Columns** to `customerid`.
   - Click **OK**.
   - **Verify:** the join result opens as a new view.
   - On the left sidebar, click the **Dashboards** icon.
   - **Verify:** no project node other than **New Dashboard** is listed.

9. **Clone a table.**
   - Click the view tab of the script result.
   - Click the table name at the left of the status bar.
   - In the **Context Panel**, expand **Actions**.
   - Click **Clone**.

10. **Save with Data sync.**
    - Click **SAVE** on the ribbon.
    - **Verify:** every table, including the derived ones, has its own
      **Data sync** toggle.
    - Enter `complex` as the name.
    - Leave **Data sync** ON for every table.
    - Click **OK**.
    - In the **Share** dialog, click **CANCEL**.

11. **Open four more tables.**
    - In **Browse > Spaces > complex**, double-click `demog.csv`.
    - In **Browse > Files > Demo > northwind**, double-click
      `customers.csv`.
    - Under **NorthwindTest**, double-click `complexQuery`.
    - In **NorthwindTest > Schemas > public**, right-click `products`
      and choose **Get All**.

12. **Add them to the project by drag and drop.**
    - On the left sidebar, click the **Dashboards** icon.
    - **Verify:** the four new tables are listed under **New
      Dashboard**.
    - Collapse the `complex` node.
    - Drag the first new table onto the `complex` node.
    - In the **Move entity** dialog, click **YES**.
    - Repeat the drag and **YES** for the other three tables.
    - Click **SAVE** next to the `complex` node.
    - In the **Save project** dialog, click **OK**.
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `complex` into the search box.
    - Click the refresh icon.
    - Double-click the `complex` tile.
    - **Verify:** the four added tables are open.

13. **Save a copy without Data sync.**
    - Click **SAVE** on the ribbon.
    - Select **Save a copy**.
    - Enter `complexNoSync` as the name.
    - Switch **Data sync** OFF for every table.
    - Click **OK**.
    - **Verify:** the open project node is `complexNoSync`.

14. **Re-save the copy with Data sync.**
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `complexNoSync` into the search box.
    - Click the refresh icon.
    - Double-click the `complexNoSync` tile.
    - Click **SAVE** on the ribbon.
    - Leave **Save original project** selected.
    - Switch **Data sync** ON for every table.
    - Click **OK**.

15. **Rename all tables (GROK-19212).**
    - Click the first table's view tab.
    - Click the table name at the left of the status bar.
    - In the **Context Panel**, expand **Actions**.
    - Click **Rename...**.
    - Add `R` to the table name.
    - Click **OK**.
    - Repeat for every table of the project.

16. **Save a copy with Data sync.**
    - Click **SAVE** on the ribbon.
    - Select **Save a copy**.
    - Enter `complexRenamed` as the name.
    - Leave **Data sync** ON for every table.
    - Click **OK**.
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `complexRenamed` into the search box.
    - Click the refresh icon.
    - Double-click the `complexRenamed` tile.
    - **Verify:** every table has its new name ending in `R`.
    - **Verify:** every table has rows.

17. **Rename the project.**
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, right-click the `complexRenamed` tile
      and choose **Rename...**.
    - Change the name to `complexRenamedV2`.
    - Click **OK**.

18. **Rename the query (github-3550).**
    - Under **NorthwindTest**, right-click `complexQuery` and choose
      **Rename...**.
    - Change the name to `complexQueryRenamed`.
    - Click **OK**.

19. **Rename the script.**
    - In **Scripts**, right-click `complexScript` and choose
      **Edit...**.
    - Change the first line to `//name: complexScriptRenamed`.
    - Click **SAVE** in the editor.

20. **Move the project into the Space.**
    - In **Browse > Dashboards**, right-click the `complexRenamedV2`
      tile and choose **Move to Space...**.
    - Set **Space** to `complex`.
    - Click **OK**.

21. **Move the query and the script into the Space.**
    - Drag `complexQueryRenamed` onto the `complex` Space in the Browse
      tree.
    - In the **Move entity** dialog, select **Move**.
    - Click **YES**.
    - Drag `complexScriptRenamed` onto the `complex` Space.
    - Select **Move**.
    - Click **YES**.

22. **Check the moved project.**
    - Go to **Browse > Spaces > complex**.
    - Double-click `complexRenamedV2`.
    - **Verify:** every table has rows.
    - **Verify:** no error dialog appears.
    - **Verify:** the pivot, aggregate and join tables are open.
    - Click **SAVE** on the ribbon.
    - **Verify:** every table shows **Data sync** ON and a **CREATION
      SCRIPT**.
    - Click **CANCEL**.
    - Right-click the left sidebar and select **Close All**.

23. **Check the first project.**
    - Go to **Browse > Dashboards**.
    - Type `complex` into the search box.
    - Click the refresh icon.
    - Double-click the `complex` tile.
    - **Verify:** every table has rows.
    - **Verify:** no error dialog appears.
    - Click **SAVE** on the ribbon.
    - **Verify:** every table shows **Data sync** ON and a **CREATION
      SCRIPT**.
    - Click **CANCEL**.
    - Right-click the left sidebar and select **Close All**.

24. **Share for viewing.**
    - In **Browse > Spaces > complex**, right-click `complexRenamedV2`
      and choose **Share...**.
    - Type the second user into **User, group, or email**.
    - Pick the second user from the suggestion list.
    - Leave **View and use** selected.
    - Click **OK**.

25. **Share with full access.**
    - In **Browse > Dashboards**, right-click the `complex` tile and
      choose **Share...**.
    - Type the second user into **User, group, or email**.
    - Pick the second user from the suggestion list.
    - Click **View and use** next to the second user.
    - In the privilege tree, tick **Full access**.
    - Click outside the tree.
    - Click **OK**.

26. **The second user opens the view-only project (GROK-18345,
    GROK-19403).**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `complexRenamedV2` into the search box.
    - Click the refresh icon.
    - Double-click the `complexRenamedV2` tile.
    - **Verify:** every table has rows, including the Space table and
      the script table.
    - Click **SAVE** on the ribbon.
    - **Verify:** **Save original project** is disabled.
    - **Verify:** **Save a copy** is selected.
    - Click **CANCEL**.
    - Right-click the left sidebar and select **Close All**.

27. **The second user saves the full-access project.**
    - Go to **Browse > Dashboards**.
    - Type `complex` into the search box.
    - Click the refresh icon.
    - Double-click the `complex` tile.
    - **Verify:** every table has rows.
    - Click **SAVE** on the ribbon.
    - **Verify:** **Save original project** is enabled and selected.
    - Click **OK**.
    - **Verify:** the balloon *Project "complex" uploaded.* appears.
    - Right-click the left sidebar and select **Close All**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

28. **Break the creation script (GROK-19728).**
    - In **Browse > Platform > Functions > Scripts**, right-click
      `complexScriptRenamed` and choose **Edit...**.
    - Add the line `throw new Error('intentional break');` before the
      `df = …` line.
    - Click **SAVE** in the editor.
    - Go to **Browse > Dashboards**.
    - Type `complexRenamedV2` into the search box.
    - Click the refresh icon.
    - Double-click the `complexRenamedV2` tile.
    - **Verify:** the **Data loading error** dialog offers **OPEN
      ANYWAY**, **EDIT SCRIPT...** and **CLOSE PROJECT**.
    - Click **CLOSE PROJECT**.

29. **The second user sees the broken script.**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Reload the browser tab.
    - Go to **Browse > Dashboards**.
    - Type `complexRenamedV2` into the search box.
    - Click the refresh icon.
    - Double-click the `complexRenamedV2` tile.
    - **Verify:** the **Data loading error** dialog says *Ask the project
      owner to fix the script*.
    - **Verify:** the dialog offers only **OPEN ANYWAY** and **CLOSE
      PROJECT**.
    - Click **CLOSE PROJECT**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

### Cleanup

- In **Browse > Dashboards**, right-click each `complex…` tile and
  choose **Delete Project**.
- Click **DELETE** and wait until the dialog closes before the next
  project.
- Under **NorthwindTest**, right-click `complexQueryRenamed` and
  choose **Delete**.
- Click **DELETE**.
- In **Scripts**, right-click `complexScriptRenamed` and choose
  **Delete**.
- Click **YES**.
- In **Browse > Spaces**, right-click `complex` and choose **Delete
  Space**.
- Click **DELETE**.

## Expected results

- Every entity works after being renamed and moved.
- Data sync reloads tables as expected.
- Users with different access levels see and can change the project
  according to their rights.
