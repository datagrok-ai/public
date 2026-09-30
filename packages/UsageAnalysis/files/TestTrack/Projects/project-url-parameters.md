---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [views.projects, views.functions]
pyramid_layer: integration
produced_from: manual-draft
related_bugs: [GROK-20930, GROK-20929, GROK-20006, GROK-20511]
---

# Project URL parameters — the dashboard link, Toolbox > Source and saving new values

A dashboard built on a query with a parameter can expose that parameter
in its link: `/p/<namespace>.<name>?<alias>=<value>` opens the dashboard
with that value. In the open result, **Toolbox > Source** shows the
query parameters, a **REFRESH** button, a copy icon for the link and a
sliders icon that chooses which parameters the link carries. Before the
result is saved, the copy icon gives a link that runs the query; once it
is saved as a dashboard, the same icon gives the dashboard's link.

This test saves such a dashboard, gives the parameter a short alias,
opens the dashboard by links with and without the parameter, saves new
values into the project and into a copy, and takes the parameter out of
the link.

## Setup

1. Log in as the test user.
2. Names in this test: query `urlParamQuery`, projects `urlParamProj`
   and `urlParamCopy`.
3. **Create the query.**
   - Go to **Browse > Databases > Postgres > Datagrok > Schemas >
     public** (on dev: **Datagrok > Catalogs > datagrok > public**).
   - Right-click `entity_types` and choose **New SQL Query...**.
   - Replace the **Name** with `urlParamQuery`.
   - Replace the query text with:

     ```
     --input: string typeName = "Project"
     select name from entity_types where name = @typeName
     ```

   - Click **Save** on the ribbon.
   - Close the query editor.

## Scenario

1. **Run the query.**
   - Under **Datagrok**, double-click `urlParamQuery`.
   - **Verify:** the `urlParamQuery` view opens with one row, `Project`.

2. **The query link before saving.**
   - Open **Toolbox > Source** and hover the copy icon next to
     **REFRESH**.
   - **Verify:** the tooltip link has `/func/`, the query name
     `urlParamQuery`, `typeName=Project` and `run=true`.

3. **Open the Save dialog.**
   - Click **SAVE** on the ribbon.
   - Enter `urlParamProj` as the name.
   - **Verify:** **Data sync** is ON for `urlParamQuery`.
   - Click **URL Parameters** under the table.
   - **Verify:** `typeName` is ticked, and its alias box reads
     `typeName`.
   - **Verify:** the **Share link:** line ends with
     `.UrlParamProj?typeName=Project`. Known bug GROK-20930: the line
     does not follow the name typed into the dialog.

4. **Give the parameter an alias.**
   - Replace the alias `typeName` with `type`.
   - **Verify:** the **Share link:** line now ends with `?type=Project`.
   - Click **OK**.
   - **Verify:** the balloon *Project "urlParamProj" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.

5. **The link after saving.**
   - In **Toolbox > Source**, hover the same copy icon next to
     **REFRESH**.
   - **Verify:** the link is now `/p/`, the project's namespace and
     `.UrlParamProj?type=Project` — the dashboard, not the query, and
     without `run=true`.
   - Click the copy icon, open a new tab, paste the link and press
     Enter.
   - **Verify:** the dashboard opens with one row, `Project`.
   - Close the tab.

6. **The sliders icon right after saving (GROK-20929).**
   - Known bug GROK-20929: the icon appears only after the project is
     reopened.
   - **Verify:** in **Toolbox > Source**, next to **REFRESH** there is a
     sliders icon; hovering it shows *Choose which parameters the
     dashboard link carries*.
   - Right-click the left sidebar and select **Close All**.

7. **Find the project link.**
   - Go to **Browse > Dashboards**.
   - Type `urlParamProj` into the search box.
   - Click the `urlParamProj` tile once.
   - In the **Context Panel**, under **Details**, click **Links...**.
   - **Verify:** **URL** is the server address, `/p/`, the project's
     namespace (the part of **Grok name** before the colon), and
     `.UrlParamProj`.
   - Click the copy icon next to **URL**.

8. **Open the link with another value.**
   - Open a new browser tab.
   - Paste the URL, add `?type=Package` to its end and press Enter.
   - **Verify:** the `urlParamQuery` view opens with one row, `Package`.
   - Open **Toolbox > Source**.
   - **Verify:** the `typeName` field reads `Package`.
   - **Verify:** no error balloon appears.
   - Close the tab.

9. **Open the link without the parameter.**
   - Open a new browser tab.
   - Paste the URL without any parameter and press Enter.
   - **Verify:** the view opens with one row, `Project` (the saved
     value).
   - Close the tab.

10. **A parameter name that is not in the link is ignored.**
    - Open a new browser tab.
    - Paste the URL, add `?typeName=Package` (the parameter's own name,
      not its alias) and press Enter.
    - **Verify:** the view opens with one row, `Project`.
    - Close the tab.

11. **Reopen the project.**
    - In **Browse > Dashboards**, double-click the `urlParamProj` tile.
    - Open **Toolbox > Source**.
    - **Verify:** the sliders icon is there.

12. **Change the value and copy the link.**
    - In **Source**, set `typeName` to `Package`. Do not click
      **REFRESH**.
    - Hover the copy icon next to **REFRESH**.
    - **Verify:** the tooltip reads *Copy link:* and a link ending with
      `.UrlParamProj?type=Package`.
    - Click **REFRESH**.
    - **Verify:** the grid has one row, `Package`.
    - Click the copy icon.
    - **Verify:** the icon briefly turns into a green check.
    - Open a new browser tab, paste the link, press Enter.
    - **Verify:** the view opens with one row, `Package`.
    - Close the tab.

13. **Save the new value into the project (GROK-20006).**
    - Click **SAVE** on the ribbon.
    - Leave **Save original project** selected.
    - Expand **CREATION SCRIPT**.
    - **Verify:** the call shows `"Package"`.
    - Click **OK**.
    - Right-click the left sidebar and select **Close All**.
    - In **Browse > Dashboards**, double-click the `urlParamProj` tile.
    - **Verify:** the grid has one row, `Package`.

14. **Save a copy with a third value.**
    - In **Source**, set `typeName` to `Project` and click **REFRESH**.
    - Click **SAVE**, select **Save a copy**, enter `urlParamCopy`, click
      **OK**.
    - Right-click the left sidebar and select **Close All**.
    - Open `urlParamCopy` from **Browse > Dashboards**.
    - **Verify:** one row, `Project`.
    - Right-click the left sidebar and select **Close All**. Open
      `urlParamProj`.
    - **Verify:** one row, `Package` (the original is unchanged).

15. **Take the parameter out of the link.**
    - In `urlParamProj`, open **Toolbox > Source** and click the sliders
      icon.
    - **Verify:** the menu lists `typeName` with a check mark.
    - Click `typeName`.
    - **Verify:** the hint *Save dashboard to apply changes* appears
      under the form.
    - Hover the copy icon.
    - **Verify:** the link in the tooltip has no `?type=`.
    - Click **SAVE**, leave **Save original project**, click **OK**.
    - **Verify:** the hint disappears.

16. **The link no longer takes the parameter.**
    - Open a new browser tab with the URL copied in step 7 plus
      `?type=Project`.
    - **Verify:** the view opens with one row, `Package` (the saved
      value; the parameter is ignored).
    - Close the tab.

## Cleanup

- Right-click the left sidebar and select **Close All**.
- In **Browse > Dashboards**, right-click the `urlParamProj` tile and
  choose **Delete Project**. Click **DELETE**. Wait until the dialog
  closes. Repeat for `urlParamCopy`.
- Under **Datagrok**, right-click `urlParamQuery` and choose **Delete**.
  Click **DELETE**.

## Expected results

- A dashboard on a query with scalar parameters exposes them in its
  link on the first save.
- The link uses the alias, not the parameter name.
- A value in the link replaces the saved one; without it the saved
  value is used.
- Before saving, the Source link runs the query; after saving, it opens
  the dashboard with its parameters.
- The Source copy icon gives the dashboard link with the current
  values, even before **REFRESH**.
- New parameter values are saved into the project and into a copy
  independently.
- A parameter taken out with the sliders icon stops working in the link
  after the dashboard is saved.

## Automation notes

- That `entity_types` has a row named `Package` on every stand was not
  checked in source; if it does not, use another value of its `name`
  column wherever this test uses `Package`.
- The project URL for step 16 is the one from **Context Panel >
  Details > Links... > URL** of the `urlParamProj` tile (step 7).
- In the Save dialog, the new alias is written to the project's table
  settings, and the open table's own copy of the settings is updated
  only when the dialog is built (`project_entity_move.dart`). That the
  Source link in step 5 already carries `type` and not `typeName` right
  after the first save was not confirmed in source; record what appears
  on the first run.
