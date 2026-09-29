---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p2
realizes: [views.projects]
pyramid_layer: integration
produced_from: manual-draft
related_bugs: []
---

# Presentation mode on a saved project

**Presentation mode** is a switch in the **Save project** dialog. A
project saved with it opens with only its views: the menus, the
toolbox, the sidebar and the status bar are hidden, and a link at the
top right returns to design mode. `?mode=presentation` in any link
opens the platform in presentation mode too. This test turns the switch
on for a project that is already saved, and for a new project saved
from the **Dashboards** panel.

## Setup

1. Log in as the test user.
2. Names in this test: projects `presentProj`, `presentPlain` and
   `presentNew`.
3. **Save two projects.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - In **Toolbox > Viewers**, click **Scatter plot**.
   - Click **SAVE** on the ribbon, enter `presentProj`, click **OK**. In
     the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.
   - Double-click `demog.csv` again, click **SAVE**, enter
     `presentPlain`, click **OK**, and **CANCEL** in **Share**.
   - Right-click the left sidebar and select **Close All**.

## Scenario

1. **The project opens in design mode.**
   - Go to **Browse > Dashboards**.
   - Type `present` into the search box.
   - Double-click the `presentProj` tile.
   - **Verify:** the top menu, the left sidebar, the toolbox and the
     status bar are shown.

2. **Turn presentation mode on.**
   - Click **SAVE** on the ribbon.
   - Leave **Save original project** selected.
   - **Verify:** the **Presentation mode** switch is shown under the
     project picture.
   - Hover **Presentation mode**.
   - **Verify:** the tooltip says only the visualization will be
     visible and that you can switch back by the "Design mode" link on
     the top right.
   - Switch **Presentation mode** ON.
   - Click **OK**.
   - Right-click the left sidebar and select **Close All**.

3. **Reopen in presentation mode.**
   - Go to **Browse > Dashboards**.
   - Type `present` into the search box.
   - Double-click the `presentProj` tile.
   - **Verify:** the grid and the scatter plot are shown; the top menu,
     the left sidebar, the toolbox and the status bar are hidden.
   - **Verify:** the balloon *Press F7 to go back to the design mode*
     appears.
   - **Verify:** the link *back to design mode* is shown at the top
     right.

4. **Back to design mode and again.**
   - Click **back to design mode**.
   - **Verify:** the menus, the sidebar, the toolbox and the status bar
     are back, and the balloon *Press F7 to go back to the presentation
     mode* appears.
   - Press F7.
   - **Verify:** presentation mode is on again.
   - Press F7.
   - **Verify:** design mode is back.
   - Right-click the left sidebar and select **Close All**.

5. **The project link opens in presentation mode.**
   - Click the `presentProj` tile once.
   - In the **Context Panel**, under **Details**, click **Links...**.
   - Click the copy icon next to **URL**.
   - Open a new browser tab, paste the URL and press Enter.
   - **Verify:** the project opens in presentation mode.
   - Close the tab.

6. **`?mode=presentation` on a project saved without it.**
   - Click the `presentPlain` tile once and copy its **URL** from
     **Links...**.
   - Open a new browser tab, paste the URL and press Enter.
   - **Verify:** the project opens in design mode.
   - Add `?mode=presentation` to the end of the address and press
     Enter.
   - **Verify:** the project opens in presentation mode.
   - Close the tab.

7. **A new project saved from the Dashboards panel.**
   - Go to **Browse > Files > Demo** and double-click `demog.csv`.
   - On the left sidebar, click the **Dashboards** icon.
   - Click **SAVE** next to **New Dashboard**.
   - Enter `presentNew` as the name.
   - Click the **Description** field and type `Presentation test`.
   - Switch **Presentation mode** ON.
   - Click **OK**.
   - In the **Share** dialog, click **OK**.
   - Right-click the left sidebar and select **Close All**.
   - Go to **Browse > Dashboards**, type `present` into the search box
     and click the refresh icon.
   - **Verify:** the `presentNew` tile shows its name and *Presentation
     test*.
   - Double-click the `presentNew` tile.
   - **Verify:** the project opens in presentation mode.
   - Click **back to design mode**.
   - Right-click the left sidebar and select **Close All**.

## Cleanup

- In **Browse > Dashboards**, delete `presentProj`, `presentPlain` and
  `presentNew` (**Delete Project** > **DELETE**, wait for each).

## Expected results

- **Presentation mode** can be turned on when saving a project that is
  already on the server; the project then opens in presentation mode,
  from **Dashboards** and from its link.
- *back to design mode* and F7 switch between the modes.
- `?mode=presentation` opens any project in presentation mode.
- A new project saved from the **Dashboards** panel with **Presentation
  mode** ON keeps its name and description and opens in presentation
  mode.

## Automation notes

- The switch in the Save dialog starts OFF every time the dialog opens,
  also for a project saved with it ON (`project_meta.dart`); the test
  does not check its starting state after step 2 and does not turn the
  mode off again.
- The help page calls the link "Design mode"; the page shows it as
  *back to design mode* (`toggle_design_mode.dart`).
