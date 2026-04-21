# Grid — Manual Test Checklist

Ручной чеклист. Не входит в автоматизацию PW.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog (or SPGI where specified)
3. Add Grid (where specified)

---

## Column and Row Resizing by Drag

*(Drag on canvas border — no DOM target for automation)*

1. Position the mouse on the right border of the AGE column header until a resize cursor appears
2. Drag the border to the right to make the AGE column wider — column width increases
3. Drag the border to the left to make the AGE column very narrow — values adapt (lose precision, eventually show as colored circles)
4. Select the AGE and HEIGHT columns (Shift+click headers), then drag the border of one of them — both selected columns resize simultaneously
5. Drag the border between two row headers to increase row height — all rows become taller and column widths adjust automatically

## Column Reordering by Drag

*(Drag on canvas header — no DOM target for automation)*

1. Drag the RACE column header and drop it before the AGE column — RACE column moves to the new position
2. Drag the HEIGHT column header to the far left — column moves to the beginning
3. Right-click the WEIGHT column header and drag the right border all the way to the left — column is hidden by width

## Block Selection by Mouse Drag

*(Canvas drag — no DOM target for automation)*

1. Click on a data cell in the AGE column and drag across to WEIGHT column and several rows down — both rows and columns get selected (rows in orange, column headers highlighted)
2. Click elsewhere to clear selection

## Column Tooltip Settings

*(Skipped from automation — too many context menu levels, no JS API for tooltip column list)*

*Dataset: demog*

1. Right-click the AGE column header, expand Tooltip > Current Column, select Default — hovering a cell shows the default tooltip
2. Right-click the AGE column header, expand Tooltip > Current Column, select None — hovering a cell shows no tooltip
3. Right-click the AGE column header, expand Tooltip > Current Column, select Custom Columns... — a column selector dialog opens; choose AGE and HEIGHT, click OK
4. Hover a data cell in the AGE column — tooltip shows AGE and HEIGHT values for that row

## Summary Columns — Form Designer

*(Form designer UI not automatable)*

*Dataset: demog*

1. Right-click any data cell, expand Add > Summary Columns, select Design a Form... — a form designer opens
2. Arrange some fields in the form designer and click Close and Apply — a form column is added to the grid
3. Click the form column header — the Context Panel shows form settings

## Column Stats — Visual Verification

*(Stats row values require visual canvas inspection)*

*Prerequisite: add Min and Max stats rows as per grid-tests-pw.md section 16*

1. Verify the Min row shows the minimum AGE, HEIGHT, and WEIGHT values
2. Verify the Max row shows the maximum AGE, HEIGHT, and WEIGHT values

## Context Panel — Column Hamburger Menu Inline Filter

*(SPGI dataset — inline filter in hamburger popup requires complex DOM interaction)*

*Dataset: SPGI*

1. Hover a column header and click the hamburger (three-lines) icon to open the column menu
2. Adjust the inline filter slider — rows should be filtered immediately
3. Click Add filter — the filter is added to the Filter Panel
4. Repeat for numeric, categorical columns

## Context Panel — Permissions

*(Requires multi-user setup — not automatable)*

*Dataset: SPGI*

1. Open the Context Panel and click any column header
2. Go to Advanced > Permissions > Edited by and enter your user name — you should have the ability to edit the column
3. Go to Advanced > Permissions > Edited by and enter only a colleague's user name — you should NOT have the ability to edit the column; a balloon appears when trying to edit

## Column Groups

*(No public JS API — cannot automate via grid-tests-pw.md)*

*Dataset: demog*

1. Right-click a column header and choose Group Columns... (or select several columns and group them)
2. Expand/collapse the group in the header
3. Remove the group via right-click > Ungroup
4. Change group properties (name, color) via the Context Panel

## Pick Up / Apply — Cross-Grid Visual Match

*(Visual verification of color match between two grid instances — not automatable)*

*Dataset: SPGI, SPGI-linked2*

1. For Average Mass set a linear Color Coding and edit the color scheme
2. Right-click the grid and select Pick Up/Apply > Pick Up
3. Open the second SPGI dataset
4. Right-click the grid and select Pick Up/Apply > Apply — color-coding, formatting, and style on both grids should match visually

## Table Switching (Grid Viewer)

*(Complex Context Panel UI — not automatable without chrome-devtools MCP)*

*Dataset: SPGI and SPGI-linked1*

1. Open SPGI and SPGI-linked1 datasets
2. Go to SPGI, click the Add viewer icon, and select Grid
3. With the added Grid viewer selected, open the Context Panel
4. Go to Data > Table and switch to SPGI-linked1 — viewer re-binds to the new table
5. Modify row height and frozen columns — changes are reflected in the viewer



---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/SPGI.csv", "System:DemoFiles/SPGI-linked1.csv"]
}
