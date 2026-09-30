---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
realizes: [charts.viewer.sankey, charts.viewer.chord]
related_bugs:
  - GROK-17772
  - GROK-18035
  - GROK-18036
  - GROK-18042
  - GROK-18048
---

# Sankey and Chord viewers

Sankey draws flows from **Source** values to **Target** values, weighted by **Value**. Chord
draws links between **From** and **To** categories. Both follow the table filter.

## Setup

A clean session. Each scenario opens its own data.

## Scenarios

### 1. Sankey offers only usable columns (GROK-18048, GROK-18042, GROK-18036)

1. Open `System:DemoFiles/demog.csv`.
2. On the Menu Ribbon, click **Add viewer** and select **Sankey**.
3. **Verify:** the Sankey is painted with **Source** SEX, **Target** RACE, **Value** AGE; no
   errors.
4. Click the **Gear** icon. Open the **Source** column list.
5. **Verify:** CONTROL, AGE and STARTED are not offered; SEX, RACE and DIS_POP are; there is no
   empty choice.
6. Open the **Target** column list.
7. **Verify:** AGE and CONTROL are not offered; there is no empty choice.
8. Set **Source** to **RACE**, **Target** to **DIS_POP**, **Value** to **WEIGHT**.
9. **Verify:** the Sankey repainted; no errors.

### 2. Sankey follows the filter, including a filter with no rows (GROK-18035)

1. Continue from scenario 1.
2. Open the Filter Panel and keep only **F** in the **SEX** filter.
3. **Verify:** 3243 rows pass the filter; the Sankey repainted; no errors.
4. Hover several Sankey flows.
5. **Verify:** no errors.
6. Keep only **true** in the **CONTROL** filter and only **Asian** in the **RACE** filter.
7. **Verify:** 0 rows pass the filter; the Sankey shows no flows; no errors.
8. Reset the filters.
9. **Verify:** 5850 rows pass the filter; the Sankey is painted.

### 3. Chord redraws on filter without a click (GROK-17772)

1. Close all. Open `System:DemoFiles/demog.csv` and add a **Chord**.
2. **Verify:** the Chord is painted with **From** SEX and **To** RACE; no errors.
3. Click the **Gear** icon. Set **From** to **RACE** and **To** to **DIS_POP**.
4. **Verify:** the Chord repainted; no errors.
5. Open the Filter Panel and keep only **Asian** in the **RACE** filter.
6. **Verify:** without clicking the viewer, the Chord has less ink than before; no errors.
7. Reset the filters.
8. **Verify:** the Chord has more ink than before; no errors.

## Cleanup

Close all.

## Expected results

- Sankey **Source** and **Target** list only string columns and cannot be left empty.
- Sankey redraws on every filter change, shows nothing when no row passes, and hovering after a
  filter raises no error.
- Chord redraws as soon as the filter changes, without a click on the viewer.

## Automation notes

- Flows and chords are SVG; the viewers have no widget status for nodes, links or chords, so the
  readings are ink, the filter count and errors. "Shows no flows" can be read as the absence of
  flow paths in the Sankey's SVG.
- Scenario 1 steps 5 and 7 need a step that reads which columns a column list offers.

---
{
  "order": 36,
  "datasets": ["System:DemoFiles/demog.csv"]
}
