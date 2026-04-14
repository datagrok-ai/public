### Tree map viewer

1. Open the **demog** dataset.
2. On the **Viewers** tab, click **Tree map**. The viewer opens with an auto-picked split column.
3. In the split row above the canvas, pick a column in the first `<select>` (e.g. `race`).
   * Expected: Canvas splits into rectangles per category.
4. Pick a second column in the trailing empty `<select>` (e.g. `sex`).
   * Expected: A new empty selector is appended; rectangles are now nested by `race` → `sex`.
5. Set the second selector back to the empty option.
   * Expected: That level is removed; the trailing selector chain is truncated.
6. Click the **Color** combo box (to the right of the splits); pick a numerical column (e.g. `age`).
   Change its aggregation `<select>` (e.g. `avg` → `max`).
   * Expected: Rectangles are recolored on a gray→red scale; aggregation change triggers a recolor.
7. Click a rectangle, then Shift+click another, then Ctrl+click a third.
   * Expected: Click replaces selection, Shift adds, Ctrl toggles. Selected rows shown in orange.
8. Hover a rectangle.
   * Expected: Row-group tooltip with the leaf name and rows in that group.
9. Click the **Gear** icon and modify properties in the **Property Pane**:
   `Size Column Name` + `Size Aggr Type`, `Outer Margin *`, `Show Column Selection Panel`,
   `Row Source`, `Filter` (e.g. `${AGE} > 40`).
   * Expected: Each change is reflected immediately without errors; hiding the selection panel
     removes the split + color row above the canvas.

---
{
  "order": 23,
  "datasets": ["System:DemoFiles/demog.csv"]
}