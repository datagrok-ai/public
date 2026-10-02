---
feature: boxplot
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  The mouse-over highlight checks are visual judgment and the tooltip walk is
  exploratory — human-eye items not covered by the automated section.
---

# Box plot tests (manual checklist)

> Manual checklist — only items NOT covered by the automated section.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Box plot

## Mouse-over highlight (visual judgment)

1. Set Value to AGE, Category 1 to SEX
2. Hover over a data point -- the point should visually highlight
3. Move to a different point -- the highlight should follow the pointer
4. Move away from all points -- the highlight should disappear
5. Hover over a category region -- the category (row group) should highlight
6. Set **Show Mouse Over Row Group** to false -- the category highlight should
   not appear on hover

## Tooltips testing (exploratory)

1. Right-click the box plot and open the **Tooltip** tab: toggle **Show tooltip** off and on, exercise the
   remaining Tooltip-tab items the same way (e.g. the group tooltip toggle),
   and edit the tooltip columns (e.g. leave only AGE and SEX) — the hover tooltip
   disappears/reappears and its rows follow the edited column list
2. Save to Layout, reload, check
3. Go to Property Pane > Tooltip and repeat the same changes there (Show tooltip off/on,
   tooltip column edits) — the hover tooltip follows each change
4. Save to Layout, reload, check

## Pick Up / Apply

1. Customize the box plot: set Value to HEIGHT, Category to RACE, change the marker color
2. Add a second Box plot with default settings
3. Right-click the customized box plot and select Pick Up / Apply > Pick Up
4. Right-click the second box plot and select Pick Up / Apply > Apply — the second viewer now matches the first (value, category, marker color)

## Cleanup

1. Delete whichever layouts this run saved during the Tooltips steps from the layout
   gallery (View > Layout > Open gallery, or Browse > Platform > Layouts), then Close all.

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
