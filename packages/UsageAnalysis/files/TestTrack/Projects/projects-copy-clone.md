---
feature: projects
ui_coverage_split_to:
  - projects-copy-clone-ui.md
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.save-copy-with-link-mode, projects.cp.url-parameterized-share]
realizes: [views.projects, sharing.share-dialog]
realized_as:
  - projects-copy-clone-spec.ts
  - project-url-spec.ts
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
  - context-panel-links-url-copy
  - new-tab-open-url
ui_coverage_delegated_to: null
related_bugs: [GROK-19750, GROK-19103, GROK-19403]
---

# Projects — Save modes: original, copy with link, copy with clone, personal view; deleting the original

A project is edited and saved in each of the save modes. Each result
must reopen correctly, from **Dashboards** and from its link. The key
check is GROK-19750: saving a copy **with link** must not change the
original project. Personal view customizations stay with their author:
the second user opens the same project without them, and **Reset** in
**Context Panel > Custom views** drops them. At the end the original is
deleted: **Delete Project** removes it for everyone together with the
tables it owns, while the file it reads stays. The copy with clone has
its own tables and keeps working; the copy with link only refers to the
original's tables, so they are gone from it too.

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
   - In the **Context Panel**, expand **Custom views**.
   - **Verify:** it says *No personal view customizations*.

2. **Share.**
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
   - Click **CLOSE**.
   - Click the filter icon on the ribbon (tooltip *Toggle filters*).
   - In the `SEX` filter, clear the check box next to `F`.
   - **Verify:** the status bar shows **Filtered: 2,607**.
   - Click **SAVE** on the ribbon.
   - Select **Save personal view customizations**.
   - **Verify:** the **Name** field is greyed out and reads `copyClone`;
     the description field is greyed out too.
   - **Verify:** the **Presentation mode** switch is not shown.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.

9. **Reopen with the customizations.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the refresh icon.
   - **Verify:** no new tile such as *Copy of copyClone* appeared: this
     save mode makes no copy.
   - Click the `copyClone` tile.
   - In the **Context Panel**, expand **Custom views**.
   - **Verify:** it says *Project has personal view customizations* and
     has a **Reset** button.
   - Double-click the `copyClone` tile.
   - **Verify:** the warning balloon *Project has personal view
     customizations* appears.
   - **Verify:** the grid is sorted by `AGE` ascending.
   - **Verify:** the `DIS_POP` column is hidden.
   - **Verify:** the `SEX` filter keeps only `M`, and the status bar
     shows **Filtered: 2,607**.
   - Click **SAVE** on the ribbon.
   - **Verify:** **Save personal view customizations** is selected.
   - Click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

10. **Open each project from its link.**
    - Go to **Browse > Dashboards**.
    - Type `copyClone` into the search box.
    - Click the `copyClone` tile once.
    - In the **Context Panel**, under **Details**, click **Links...**.
    - **Verify:** the dialog **Links to copyClone** shows **ID**, **Grok
      name**, **Markup** and **URL**.
    - Click the copy icon next to **URL**.
    - Open a new browser tab.
    - Paste the URL into the address bar.
    - **Verify:** the URL is the server address, then `/p/`, the
      project's namespace (the part of **Grok name** before the colon)
      and `.CopyClone`.
    - Press Enter.
    - **Verify:** the new tab shows the grid, the bar chart and the
      scatter plot, and no line chart.
    - **Verify:** your personal view customizations are applied: the
      grid is sorted by `AGE` ascending and `DIS_POP` is hidden.
    - **Verify:** no error balloon appears.
    - Close the tab.
    - Repeat for `copyCloneLink`: the URL ends with `.CopyCloneLink`, and
      the tab shows the grid, the line chart, the bar chart and the
      scatter plot.
    - Repeat for `copyCloneClone`: the URL ends with `.CopyCloneClone`,
      and the tab shows the grid, the histogram, the bar chart and the
      scatter plot.

11. **Share the copies.**
    - Go to **Browse > Dashboards**.
    - Type `copyCloneLink` into the search box.
    - Right-click the `copyCloneLink` tile and choose **Share...**.
    - Type the second user into **User, group, or email**.
    - Pick the second user from the suggestion list.
    - Leave **View and use** selected.
    - Click **OK**.
    - Type `copyCloneClone` into the search box.
    - Right-click the `copyCloneClone` tile and choose **Share...**.
    - Type the second user into **User, group, or email**.
    - Pick the second user from the suggestion list.
    - Leave **View and use** selected.
    - Click **OK**.

12. **The second user opens the projects.**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Go to **Browse > Dashboards**.
    - Type `copyClone` into the search box.
    - Click the refresh icon.
    - Click the `copyClone` tile.
    - In the **Context Panel**, expand **Custom views**.
    - **Verify:** it says *No personal view customizations*.
    - Double-click the `copyClone` tile.
    - **Verify:** it shows the grid, the bar chart and the scatter plot.
    - **Verify:** `demog` has 5,850 rows and no error dialog appears.
    - **Verify:** the owner's personal view customizations are not
      applied: the grid is not sorted by `AGE`, the `DIS_POP` column is
      visible, and no rows are filtered out.
    - **Verify:** the balloon *Project has personal view customizations*
      does not appear.
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

13. **Reset the personal view customizations.**
    - Go to **Browse > Dashboards**.
    - Type `copyClone` into the search box.
    - Click the refresh icon.
    - Click the `copyClone` tile.
    - In the **Context Panel**, expand **Custom views** and click
      **Reset**.
    - **Verify:** the text and the **Reset** button disappear.
    - Click the `copyCloneLink` tile, then the `copyClone` tile again.
    - **Verify:** **Custom views** says *No personal view
      customizations*.

14. **The project opens as saved.**
    - Double-click the `copyClone` tile.
    - **Verify:** the balloon *Project has personal view customizations*
      does not appear.
    - **Verify:** the grid is not sorted by `AGE`, and `DIS_POP` is
      visible.
    - Click **SAVE** on the ribbon.
    - **Verify:** **Save original project** is selected.
    - Click **CANCEL**.
    - Right-click the left sidebar and select **Close All**.

15. **The copy with link refers to the original's table.**
    - In **Browse > Dashboards**, type `copyClone` into the search box
      and click the refresh icon.
    - Click the `copyCloneLink` tile.
    - In the **Context Panel**, expand **Content**.
    - **Verify:** `demog` is listed with a link icon; hovering it shows
      *This entity is not included to this project, but linked.*

16. **Delete the original.**
    - Right-click the `copyClone` tile and choose **Delete Project**.
    - **Verify:** the dialog *Are you sure? Delete project "copyClone"?*
      opens.
    - Click **DELETE**.
    - Wait until the dialog closes.
    - Click the refresh icon.
    - **Verify:** the `copyClone` tile is gone; `copyCloneLink` and
      `copyCloneClone` are still shown.

17. **The file stays.**
    - Go to **Browse > Files > Demo**.
    - **Verify:** `demog.csv` is still listed.

18. **The copy with clone still works.**
    - In **Browse > Dashboards**, type `copyClone` into the search box
      and double-click the `copyCloneClone` tile.
    - **Verify:** `demog` opens with 5,850 rows.
    - **Verify:** no error dialog appears.
    - Right-click the left sidebar and select **Close All**.

19. **The copy with link lost the original's table.**
    - Click the `copyCloneLink` tile.
    - In the **Context Panel**, expand **Content**.
    - **Verify:** `demog` is no longer listed.
    - Double-click the `copyCloneLink` tile.
    - **Verify:** the project does not open as an empty dashboard
      without a word: a dialog or a balloon says that some of its data
      could not be loaded.
    - Close any dialog, then right-click the left sidebar and select
      **Close All**.

## Cleanup

- In **Browse > Dashboards**, right-click the `copyCloneLink` tile and
  choose **Delete Project**.
- Click **DELETE**.
- Wait until the dialog closes.
- Repeat for `copyCloneClone`, and for `copyClone` if the run stopped
  before step 16.

## Expected results

- Each save mode produces a project that reopens correctly, from
  **Dashboards** and from the link in **Links...**.
- Saving a copy with link does not change the original (GROK-19750).
- A copy with clone has its own copy of the data.
- **Save personal view customizations** makes no copy and keeps the
  project name. The sort, the hidden column and the filter come back
  for their author, also when the project is opened from its link.
- The second user opens the same project without the owner's personal
  view customizations.
- Shared copies open for the recipient.
- **Custom views** shows whether you have personal view customizations
  for the project.
- **Reset** removes them; the project then opens as its owner saved it,
  and the Save dialog no longer starts in the personal mode.
- Deleting a project removes it for everyone; the file it read stays.
- A copy with clone is independent and opens after the original is
  deleted.
- A copy with link loses the tables it referred to, and says so when
  opened.

## Automation notes

- The ribbon filter icon, its tooltip *Toggle filters*, the **Custom
  views** texts, the **Reset** button and the balloon *Project has
  personal view customizations* are read from `table_view.dart`,
  `project_meta.dart` and `shell_project.dart`. The `SEX` filter's
  check-box gesture was not checked in source.
- Step 13 needs "another tile" in the gallery: the search `copyClone`
  also lists `copyCloneLink`; any other tile shown works too.
- That the original's tables are deleted with it follows from the
  server code (`projects_repository.dart`: tables and views in the
  project's namespace are deleted) and from the help page
  `help/datagrok/concepts/project/dashboard.md`. What the copy with link
  shows on opening (step 19) was not confirmed from source; record what
  appears on the first run.
