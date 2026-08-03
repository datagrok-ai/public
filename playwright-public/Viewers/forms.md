# Forms tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Forms

## Fields selection

1. Open Context Panel and find the Forms viewer properties
2. Click on **Fields** property — a column selection dialog should open
3. Uncheck all columns, then check only AGE, SEX, and RACE
4. Confirm — the forms should show only AGE, SEX, and RACE fields
5. Drag-and-drop RACE above AGE in the column selection dialog — RACE should appear first in each form card
6. Remove a column by clicking the X icon next to its name in the column headers — that column should disappear from all form cards

## Current row tracking

1. Click on different rows in the grid — a form card with a green stripe on top should update to reflect the current row
2. Verify the current row form is always shown first (topmost position)
3. Set **Show Current Row** to false in the properties — the green-striped card should disappear
4. Set **Show Current Row** back to true — the current row card reappears

## Mouse-over row tracking

1. Hover over different rows in the grid — a form card with a grey stripe should appear showing the hovered row
2. Verify the mouse-over form is positioned right after the current row form
3. Move the mouse away from the grid — the mouse-over card should remain showing the last hovered row
4. Set **Show Mouse Over Row** to false — the grey-striped card should disappear
5. Set **Show Mouse Over Row** back to true — the card reappears on hover

## Selected rows display

1. In the grid, select 3 rows (Ctrl+click) — 3 form cards should appear below the current/mouse-over cards
2. Select additional rows — new form cards should appear
3. Deselect a row (Ctrl+click) — its form card should disappear
4. Set **Show Selected Rows** to false — all selected-row cards should disappear, only current/mouse-over remain
5. Set **Show Selected Rows** back to true — selected-row cards reappear

## Form card click interactions

1. Click on a form card (not the current row one) — that row should become the current row in the grid
2. Ctrl+click a form card — that row's selection should toggle
3. Ctrl+click the same card again — selection should toggle back
4. Click on a specific field value inside a card — that cell should become the current cell in the grid

## Color coding

1. In the grid, apply color coding to the AGE column (e.g., Linear color coding)
2. Verify AGE values in the form cards reflect the color coding (colored backgrounds or text)
3. Set **Color Code** to false in the Forms properties — color coding should disappear from all cards
4. Set **Color Code** back to true — colors reappear

## Grid sort synchronization

1. In the grid, sort by HEIGHT column
2. Verify the selected-row forms follow the grid sort order
3. Set **Use Grid Sort** to false in properties — forms should revert to natural row order
4. Set **Use Grid Sort** back to true — forms should follow grid sort again

## Sort By property

1. Set the **Sort By** property to WEIGHT — selected-row cards should reorder by WEIGHT
2. Change **Sort By** to AGE — cards reorder by AGE
3. Clear the **Sort By** property — sorting falls back to grid sort or natural order

## Renderer size

1. Set **Renderer Size** to "small" — form cards and field values should be compact
2. Set **Renderer Size** to "normal" — cards should be slightly larger
3. Set **Renderer Size** to "large" — cards and rendered content should be noticeably larger
4. Verify column name labels remain readable at each size

## Filtering interaction

1. Open the filter panel and filter SEX to show only "M"
2. Verify the Forms viewer only shows form cards for filtered rows
3. Select rows in the grid — only filtered and selected rows should appear as form cards
4. Remove the filter — previously hidden selected rows should now appear as form cards

## Column removal reaction

1. In the grid, delete a column that is displayed in the Forms viewer
2. Verify the Forms viewer removes that column from all form cards without error
3. The column should also disappear from the column headers on the left

## Layout persistence

1. Configure the Forms viewer: select specific columns, set Sort By to AGE, set Renderer Size to large
2. Save the layout (**View | Layout | Save as**)
3. Close all and reopen demog
4. Apply the saved layout
5. Verify the Forms viewer restores with the same columns, sort, and renderer size

## Molecule rendering (spgi-100)

Setup: Close all, open System:AppData/Chem/tests/spgi-100.csv, add Forms

1. Set **Fields** to Structure, Primary Series Name, Average Mass, TPSA
2. Verify the Structure column renders as a molecule drawing (not raw molblock text) on each form card
3. The molecule should occupy the left portion of the card, with other fields on the right
4. Select 5 rows — verify all 5 form cards render molecules correctly
5. Set **Renderer Size** to "large" — molecule drawings should scale up and be more detailed
6. Set **Renderer Size** to "small" — molecules should still be recognizable at a smaller size
7. Hover over different rows — the mouse-over card should show the correct molecule for the hovered row

## Multiple molecule columns (spgi-100)

Setup: Close all, open System:AppData/Chem/tests/spgi-100.csv, add Forms

1. Set **Fields** to Structure, Core, Primary Series Name
2. Verify both Structure and Core render as molecule drawings
3. Both molecules should appear on the left side of the card, stacked vertically
4. Each molecule should have its column name label displayed above it

## Curves rendering (curves.csv)

Setup: Close all, open System:DemoFiles/curves.csv, add Forms

1. Set **Fields** to smiles, multiple prefit
2. Verify the smiles column renders as a molecule drawing
3. Verify the multiple prefit column renders as a dose-response curve chart (not raw JSON text)
4. Select 3 rows — each form card should show its own molecule and curve
5. Set **Fields** to smiles, multiple styled series, styled proprtional with IC50
6. Verify both curve columns render as separate curve charts with correct styling (colors, line styles)
7. Set **Renderer Size** to "large" — curve charts should be larger and more readable
8. Set **Renderer Size** to "small" — curves should still be visible, axes may be simplified

---
{
  "order": 9,
  "datasets": ["System:DemoFiles/demog.csv,System:AppData/Chem/tests/spgi-100.csv,System:DemoFiles/curves.csv"]
}
