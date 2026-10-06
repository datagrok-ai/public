---
feature: projects
ui_coverage_split_to:
  - projects-copy-clone-ui.md
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.save-copy-with-link-mode]
realizes: [views.projects, sharing.share-dialog]
realized_as:
  - projects-copy-clone-spec.ts
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-copy-clone.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - all-created-projects-in-step-1
  - save-copy-ui-control-location
  - save-with-personal-view-customizations-mode-trigger
  - recipient-open-verification-depth-step-5
  - grok-19750-invariant-which-viewers-count-as-intact
  - grok-19403-cross-cutting-cover-via-step-5
scope_reductions: []
ui_companion: projects-copy-clone-ui.md
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - save-copy-with-link-dialog
  - save-copy-with-clone-dialog
  - save-personal-view-customizations-dialog
  - pcmdShareProject
  - share-dialog-recipients
  - context-panel-sharing-tab
ui_coverage_delegated_to: null
related_bugs: [GROK-19750, GROK-19103, GROK-19403]
---

# Projects — Save modes: original, copy with link, copy with clone, personal view

A project is edited and saved in each of the save modes. Each result
must reopen correctly. The key check is GROK-19750: saving a copy
**with link** must not change the original project.

The three projects this test leaves behind are used by
`project-url.md`. Run that test next, or delete the projects yourself.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. Projects in this test: `copyClone` and its copies `copyCloneLink`
   and `copyCloneClone`.
3. **Create the original.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - In **Toolbox > Viewers**, click **Bar chart**.
   - Click **SAVE** on the ribbon.
   - Enter `copyClone` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "copyClone" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

## Scenario

1. **Preview.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the `copyClone` tile.
   - **Verify:** the tile shows a thumbnail.
   - **Verify:** the **Context Panel** shows the name, the author and
     **Content** with `demog`.

2. **Share.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Right-click the `copyClone` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.
   - Click the `copyClone` tile.
   - In the **Context Panel**, expand **Sharing**.
   - **Verify:** the second user is listed.

3. **Save the original.**
   - Double-click the `copyClone` tile.
   - In **Toolbox > Viewers**, click **Scatter plot**.
   - Click **SAVE** on the ribbon.
   - Leave **Save original project** selected.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.

4. **Save a copy with link.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyClone` tile.
   - In **Toolbox > Viewers**, click **Line chart**.
   - Click **SAVE** on the ribbon.
   - Select **Save a copy**.
   - **Verify:** the name changes to *Copy of copyClone*.
   - Enter `copyCloneLink` as the name.
   - Next to `demog`, switch **Clone** to **Link**.
   - **Verify:** the hint *Local data changes will not be saved*
     appears.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.

5. **Check the copy with link.**
   - Go to **Browse > Dashboards**.
   - Type `copyCloneLink` into the search box.
   - Click the refresh icon.
   - Double-click the `copyCloneLink` tile.
   - **Verify:** it shows the grid, the line chart, the bar chart and
     the scatter plot.
   - **Verify:** the grid has rows.
   - Right-click the left sidebar and select **Close All**.

6. **Check the original (GROK-19750).**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyClone` tile.
   - **Verify:** it shows the grid, the bar chart and the scatter plot.
   - **Verify:** there is no line chart.
   - Right-click the left sidebar and select **Close All**.

7. **Save a copy with clone.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyClone` tile.
   - In **Toolbox > Viewers**, click **Histogram**.
   - Click **SAVE** on the ribbon.
   - Select **Save a copy**.
   - Enter `copyCloneClone` as the name.
   - Leave **Clone** selected next to `demog`.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `copyCloneClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyCloneClone` tile.
   - **Verify:** it shows the histogram, the bar chart and the scatter
     plot.
   - **Verify:** `demog` has 5,850 rows.
   - Right-click the left sidebar and select **Close All**.

8. **Save personal view customizations.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyClone` tile.
   - Right-click the `AGE` column header and choose **Sort >
     Ascending**.
   - Right-click any column header and choose **Order or Hide
     Columns...**.
   - Uncheck `DIS_POP`.
   - Click **OK**.
   - Click **SAVE** on the ribbon.
   - Select **Save personal view customizations**.
   - **Verify:** the dialog has no name field.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - Double-click the `copyClone` tile.
   - **Verify:** the grid is sorted by `AGE` ascending.
   - **Verify:** the `DIS_POP` column is hidden.
   - Right-click the left sidebar and select **Close All**.

9. **Share the copies.**
   - Go to **Browse > Dashboards**.
   - Type `copyCloneLink` into the search box.
   - Right-click the `copyCloneLink` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.
   - Go to **Browse > Dashboards**.
   - Type `copyCloneClone` into the search box.
   - Right-click the `copyCloneClone` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

10. **The second user opens the projects.**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `copyClone` into the search box.
    - Click the refresh icon.
    - Double-click the `copyClone` tile.
    - **Verify:** it shows the grid, the bar chart and the scatter plot.
    - **Verify:** `demog` has 5,850 rows and no error dialog appears.
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `copyCloneLink` into the search box.
    - Click the refresh icon.
    - Double-click the `copyCloneLink` tile.
    - **Verify:** it shows the grid, the line chart, the bar chart and
      the scatter plot.
    - **Verify:** `demog` has 5,850 rows and no error dialog appears.
    - Right-click the left sidebar and select **Close All**.
    - Go to **Browse > Dashboards**.
    - Type `copyCloneClone` into the search box.
    - Click the refresh icon.
    - Double-click the `copyCloneClone` tile.
    - **Verify:** it shows the grid, the histogram, the bar chart and
      the scatter plot.
    - **Verify:** `demog` has 5,850 rows and no error dialog appears.
    - Right-click the left sidebar and select **Close All**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

## Expected results

- Each save mode produces a project that reopens correctly.
- Saving a copy with link does not change the original (GROK-19750).
- A copy with clone has its own copy of the data.
- Personal view customizations are restored for their author.
- Shared copies open for the recipient.
