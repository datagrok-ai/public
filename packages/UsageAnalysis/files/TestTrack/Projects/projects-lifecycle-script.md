---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [projects.cp.share-with-unshared-deps, projects.cp.view-and-use-failure-state, projects.op.rename_external_dep, projects.op.share_with_recipient_open, projects.op.rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-script-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-script.md
migration_date: 2026-05-04
related_bugs:
  - GROK-19403
  - GROK-19728
---

# Projects — lifecycle of a script-based project

A project built from the output of your own script is saved and
shared **without** sharing the script. Then the script is renamed and
broken. Two bugs are checked:

- **GROK-19403.** A recipient of a project whose script was not
  shared must never get an empty table without a message. Today the
  project share also gives the recipient access to the script, so
  the data loads.
- **GROK-19728.** A recipient with **View and use** access must not be
  able to edit the creation script, even when it fails.

## Setup

1. Two accounts: the **owner** (test user) and a **second user**.
2. Names in this test: script `lifecycleScript`, project
   `lifecycleScriptProj`.

## Scenario

1. **Create the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click **NEW** and choose **JavaScript Script...**.
   - **Verify:** the editor opens with a *Hello World* template.
   - Replace the whole text with the text below.
   - Click **SAVE** in the editor.
   - **Verify:** the balloon *Script saved.* appears.

   ```
   //name: lifecycleScript
   //language: javascript
   //output: dataframe df
   df = await grok.data.getDemoTable('demog.csv');
   ```

2. **Run it into the workspace.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Click the refresh icon.
   - Right-click `lifecycleScript` and choose **Run...**.
   - **Verify:** the `demog` view opens with 5,850 rows.

3. **Save the project.**
   - Click **SAVE** on the ribbon.
   - **Verify:** the dialog lists `demog` with a **CREATION SCRIPT** and
     the note *Some tables require this script for data sync.*
   - Enter `lifecycleScriptProj` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

4. **Share the project only.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleScript` into the search box.
   - Right-click the `lifecycleScriptProj` tile and choose **Share...**.
   - Type the second user into **User, group, or email**.
   - Pick the second user from the suggestion list.
   - Leave **View and use** selected.
   - Click **OK**.

5. **The recipient opens it (GROK-19403).**
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the second user's credentials.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleScriptProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleScriptProj` tile.
   - **Verify:** the table opens with 5,850 rows.
   - **Verify:** no error dialog appears.
   - Right-click the left sidebar and select **Close All**.
   - Click your avatar at the bottom of the left sidebar.
   - Click **Logout** in the profile view.
   - Sign in with the owner's credentials.

6. **Rename the script.**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Right-click `lifecycleScript` and choose **Edit...**.
   - Change the first line to `//name: lifecycleScriptRenamed`.
   - Click **SAVE** in the editor.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleScriptProj` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleScriptProj` tile.
   - **Verify:** `demog` opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.

7. **Rename the project.**
   - Go to **Browse > Dashboards**.
   - Right-click the `lifecycleScriptProj` tile and choose
     **Rename...**.
   - Change the name to `lifecycleScriptProjRenamed`.
   - Click **OK**.
   - Go to **Browse > Dashboards**.
   - Type `lifecycleScriptProjRenamed` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleScriptProjRenamed` tile.
   - **Verify:** `demog` opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.

8. **Break the script (GROK-19728).**
   - Go to **Browse > Platform > Functions > Scripts**.
   - Right-click `lifecycleScriptRenamed` and choose **Edit...**.
   - Add the line `throw new Error('intentional break');` before the
     `df = …` line.
   - Click **SAVE** in the editor.

9. **The owner opens the broken project.**
   - Go to **Browse > Dashboards**.
   - Type `lifecycleScriptProjRenamed` into the search box.
   - Click the refresh icon.
   - Double-click the `lifecycleScriptProjRenamed` tile.
   - **Verify:** the **Data loading error** dialog says the project
     *could not load some of its data* and shows `Error: intentional
     break`.
   - **Verify:** the dialog offers **OPEN ANYWAY**, **EDIT SCRIPT...**
     and **CLOSE PROJECT**.
   - Click **CLOSE PROJECT**.

10. **The recipient opens the broken project.**
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the second user's credentials.
    - Reload the browser tab.
    - Go to **Browse > Dashboards**.
    - Type `lifecycleScriptProjRenamed` into the search box.
    - Click the refresh icon.
    - Double-click the `lifecycleScriptProjRenamed` tile.
    - **Verify:** the **Data loading error** dialog says *Ask the
      project owner to fix the script*.
    - **Verify:** the dialog offers only **OPEN ANYWAY** and **CLOSE
      PROJECT**.
    - Click **CLOSE PROJECT**.
    - Click your avatar at the bottom of the left sidebar.
    - Click **Logout** in the profile view.
    - Sign in with the owner's credentials.

11. **Cleanup.**
    - Go to **Browse > Dashboards**.
    - Right-click the `lifecycleScriptProjRenamed` tile and choose
      **Delete Project**.
    - Click **DELETE**.
    - Wait until the dialog closes (about 20 seconds).
    - In **Scripts**, right-click `lifecycleScriptRenamed` and choose
      **Delete**.
    - **Verify:** the dialog asks *Delete script
      "lifecycleScriptRenamed"?*.
    - Click **YES**.

## Expected results

- A script-based project reopens and re-runs the script.
- The recipient gets the data although the script is not shared.
- A recipient with **View and use** access cannot edit the creation
  script, even when it fails.
- Renaming the script and then the project leaves the project
  working.
