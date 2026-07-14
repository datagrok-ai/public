Ручной чеклист. Не входит в автоматизацию PW.

# Statistics viewer — manual tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Statistics

## Statistics for date columns

1. Locate the STARTED row in the Statistics viewer
2. Verify that avg, sum, stdev columns are empty for STARTED
3. Verify that values, nulls, unique columns are populated for STARTED

## Add min and max for date column via submenu

1. Locate the STARTED row in the Statistics viewer
2. Right-click the viewer and hover over Statistics in the context menu
3. Click min in the submenu — a min column shows a date value for STARTED
4. Right-click the viewer and hover over Statistics
5. Click max in the submenu — a max column shows a date value for STARTED

## Column visibility

1. Open the Properties panel for the Statistics viewer
2. In the column manager on the left side of the viewer, uncheck HEIGHT
3. Verify HEIGHT disappears from the viewer rows
4. Re-check HEIGHT — it reappears with correct statistics

## Background color

1. Open the Properties panel for the Statistics viewer
2. Set Back Color to a light grey — viewer background updates immediately
3. Set Back Color back to white

## Reorder stat columns

1. Drag a stat column header (e.g., max) left past min — column order updates
2. Verify no stat values are clipped or overlapping after reorder

## Resize stat columns

1. Drag the edge of a stat column header to widen it
2. Drag the edge of the name column to adjust width
3. Verify no stat values are clipped or overlapping after resize

## Property Pane — general

1. Click the gear icon on the Statistics viewer — Properties panel opens
2. Modify various properties in the Properties panel — all changes apply immediately in the viewer

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
