---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: [views.projects, views.space]
realized_as:
  - complex-move-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Move a project into a Space

Checks that a saved project can be moved into a Space and, from
there, into another Space, and that it still opens with its data after
each move.

## Setup

1. Log in as the test user.
2. Names in this test: project `moveTest`, Spaces `moveA` and
   `moveB`.

## Scenario

1. **Create two Spaces.**
   - In **Browse**, right-click **Spaces** and choose **Create
     Space...**.
   - In the **Create Space** dialog, replace *New Space* with `moveA`.
   - Click **OK**.
   - Right-click **Spaces** and choose **Create Space...**.
   - Replace *New Space* with `moveB`.
   - Click **OK**.
   - Expand **Spaces**.
   - **Verify:** `moveA` and `moveB` are listed.

2. **Save the project.**
   - Go to **Browse > Files > Demo**.
   - Double-click `demog.csv`.
   - Click **SAVE** on the ribbon.
   - Enter `moveTest` as the name.
   - Leave **Data sync** ON.
   - Click **OK**.
   - **Verify:** the balloon *Project "moveTest" uploaded.* appears.
   - In the **Share** dialog, click **CANCEL**.
   - Right-click the left sidebar and select **Close All**.

3. **Move it to the first Space.**
   - Go to **Browse > Dashboards**.
   - Type `moveTest` into the search box.
   - Right-click the `moveTest` tile and choose **Move to Space...**.
   - In the **Move to space** dialog, set **Space** to `moveA`.
   - Click **OK**.
   - Click **Spaces > moveA**.
   - **Verify:** `moveTest` is listed in `moveA`.

4. **Open it from the Space.**
   - Double-click `moveTest` in `moveA`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - **Verify:** no error balloon appears.
   - Right-click the left sidebar and select **Close All**.

5. **Move it to the second Space.**
   - In `moveA`, right-click `moveTest` and choose **Move to Space...**.
   - Set **Space** to `moveB`.
   - Click **OK**.
   - Click **Spaces > moveB**.
   - **Verify:** `moveTest` is listed in `moveB`.
   - Click **Spaces > moveA**.
   - **Verify:** `moveTest` is not listed in `moveA`.

6. **Open it again.**
   - Click **Spaces > moveB**.
   - Double-click `moveTest`.
   - **Verify:** the `demog` view opens with 5,850 rows.
   - Right-click the left sidebar and select **Close All**.

7. **Cleanup.**
   - In `moveB`, right-click `moveTest` and choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Right-click `moveA` and choose **Delete Space**.
   - **Verify:** the dialog says *Delete space "moveA"? This will
     delete space and its related data…*.
   - Click **DELETE**.
   - Right-click `moveB` and choose **Delete Space**.
   - Click **DELETE**.

## Expected results

- A project can be moved into a Space and between Spaces.
- The project opens with its data after every move.
