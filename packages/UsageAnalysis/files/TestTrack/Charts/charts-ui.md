---
feature: charts
realizes_atlas: [charts.cp.click-segment-to-select-or-filter]
realizes: [charts.viewer.tree]
priority: p0
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  The Tree viewer does not report where its branches are, so branch clicks
  with Shift must be done by hand.
related_bugs: []
---

# Charts: Tree branch clicks (manual)

This file lists only checks needing a human hand on the Tree viewer's branches.

## Setup

1. Open `System:DemoFiles/demog.csv`.
2. On the Menu Ribbon, click **Add viewer** and select **Tree**.
3. Click the **Gear** icon of the Tree. In the Context Panel, click the **...** button of
   **Hierarchy**. In the **Select columns...** dialog click **None**, check **CONTROL**, **SEX**
   and **RACE**, and drag the rows so the order is CONTROL, SEX, RACE. Click **OK**.

## Scenarios

### 1. Shift+Click selects several branches

1. In the Tree, hold **Shift** and click these three branches:
   - `All → false → F → Asian`
   - `All → false → F → Black`
   - `All → false → M → Asian`
2. **Verify:** the status bar shows 174 selected rows, and the grid highlights them.
3. Hold **Shift** and click the branch `All → true → F → Black`.
4. **Verify:** the status bar shows 176 selected rows.

### 2. Shift+Click extends the selection under a filter

Continue from scenario 1 after step 2 (174 rows selected).

1. In the **Filter Panel**, keep only **true** in the **CONTROL** filter.
2. **Verify:** no selected row passes the filter (0 filtered and selected).
3. Hold **Shift** and click the Tree branch `All → true → F → Black`.
4. **Verify:** 2 selected rows pass the filter.
5. In the **Filter Panel**, remove the **CONTROL** filter card.
6. **Verify:** the status bar shows 176 selected rows.

## Cleanup

Close all.

## Expected results

- Shift+Click on a branch adds that branch's rows to the selection.
- A branch picked under a filter adds only the rows under it; the earlier selection stays.

## Automation notes

- Tree branch clicks stay manual until the Tree viewer reports its branches as areas.

---
{
  "order": 35,
  "datasets": ["System:DemoFiles/demog.csv"]
}
