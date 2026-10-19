---
feature: charts
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
realizes_atlas: [charts.int.timelines-legend-click-to-filter-stability]
realizes: [charts.viewer.timelines]
realized_as:
  - timelines-spec.ts
related_bugs:
  - GROK-18608
  - GROK-19033
  - GROK-19034
  - GROK-19535
  - GROK-20800
---

# Timelines viewer

The Timelines viewer draws one lane per subject (**Split By**) with an interval per event from
**Start** to **End**, colored by **Color**. Its legend narrows what the viewer draws; the table
filter does not change.

## Setup

1. Open `System:AppData/Charts/ae.csv` (143 adverse events, 71 subjects; AESOC has 15 values:
   SKIN AND SUBCUTANEOUS TISSUE DISORDERS 35, NERVOUS SYSTEM DISORDERS 33, GASTROINTESTINAL
   DISORDERS 23, ...).
2. On the Menu Ribbon, click **Add viewer** and select **Timelines**.

## Scenarios

### 1. The viewer draws, and the legend lists the colors (GROK-20800, GROK-19034)

1. **Verify:** the Timelines viewer is painted; no errors.
2. Click the **Gear** icon. In the Context Panel, open the **Color** column list.
3. **Verify:** only string columns are offered (no AESTDY, AEENDY, AESEQ).
4. Set **Color** to **AESOC**.
5. **Verify:** the legend is visible and lists 15 items.
6. Set **Legend Visibility** to **Always**.
7. **Verify:** the legend is visible and lists 15 items.
8. Set **Legend Visibility** to **Never**.
9. **Verify:** the legend is hidden.
10. Set **Legend Visibility** to **Auto**.
11. **Verify:** the legend is visible and lists 15 items.

### 2. Legend clicks narrow the viewer without blanking it (GROK-19033, GROK-19535, GROK-18608)

1. Continue with **Color** = **AESOC**.
2. Click `SKIN AND SUBCUTANEOUS TISSUE DISORDERS` in the legend.
3. **Verify:** the viewer has less ink than before and is still painted; no errors; 143 rows pass
   the table filter.
4. Ctrl+Click `NERVOUS SYSTEM DISORDERS` in the legend.
5. **Verify:** the viewer has more ink than before; no errors.
6. Ctrl+Click `NERVOUS SYSTEM DISORDERS` again, then click
   `SKIN AND SUBCUTANEOUS TISSUE DISORDERS` again, so no item is chosen.
7. **Verify:** the viewer is painted with all events again; no errors.
8. Set **Split By** to **AESEV**.
9. **Verify:** the viewer repainted; the legend still lists 15 items; no errors.
10. Click `CARDIAC DISORDERS` in the legend.
11. **Verify:** the viewer is painted; no errors; 143 rows pass the table filter.
12. Set **Split By** back to **USUBJID**.
13. **Verify:** the viewer repainted; no errors.

### 3. Reset View

1. Right-click the Timelines viewer and choose **Reset View**.
2. **Verify:** the viewer repainted; no errors.

## Cleanup

Close all.

## Expected results

- The Timelines viewer draws on `ae.csv` without errors; **Color** offers only string columns.
- The legend follows **Legend Visibility**: shown for Auto and Always, hidden for Never.
- A legend click narrows the drawn events, Ctrl+Click adds or removes an item, and clicking the
  only chosen item shows all events again; the viewer never goes blank and the table filter
  stays at 143 rows.
- Changing **Split By** keeps the legend and its filtering working.
- **Reset View** redraws the viewer without errors.

## Automation notes

- The legend is the standard DOM legend, so its items can be counted and clicked by name.
- Ink comparisons are the only reading of the drawn events: the viewer has no widget status for
  its lanes and intervals.

---
{
  "order": 31,
  "datasets": ["System:AppData/Charts/ae.csv"]
}
