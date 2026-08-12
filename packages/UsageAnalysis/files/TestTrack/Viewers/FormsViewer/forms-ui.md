# Forms viewer tests (manual / exploratory)

> Exploratory checklist — manual steps that are NOT covered by the automated
> section

The first scenario starts on demog:
1. Close all
2. Open demog
3. Add Forms
4. Confirm **Show Selected Rows**, **Show Current Row**, and **Show Mouse Over
   Row** are all ON

The remaining scenarios name their own dataset —
**System:AppData/Chem/tests/spgi-100.csv** (table `spgi-100`) for the molecular
checks and **System:DemoFiles/curves.csv** for the curve checks.

## Mouse-over row binding follows and clears with the pointer (exploratory)

1. Ctrl+click a handful of rows (e.g. rows 2, 5, 10, 15) and click row 5 to make
   it current, so the card area has selected-row cards.
2. Move the pointer onto one of the selected-row cards (not the leading
   current-row card), then onto a DIFFERENT card, then onto a grid row, moving
   slowly. Watch the mouse-over card: it should leave each previous target and
   follow the pointer continuously, without flicker, stalling, or a card that
   keeps a stale highlight after the pointer has left it.

## Molecules occupy the left part of the card (exploratory)

Setup: Close all, open **System:AppData/Chem/tests/spgi-100.csv**, add **Forms**,
set **Renderer Size** to **normal**.

1. Open the **Fields** picker and set the field set to **Structure**, **Primary
   Series Name**, **Average Mass**, **TPSA**. Look at any card: the structure
   drawing takes the left part of the card and the three text fields sit to its
   right, so the card reads as "one picture plus its properties" rather than as a
   column of unrelated boxes.
2. Set the field set to **Structure**, **Core**, **Primary Series Name**. Both
   molecule drawings sit on the left, stacked vertically, each next to its own
   column label in the header pane — the two structures are visually separable
   and neither overlaps the other or the text fields.
3. Ctrl+click five rows in the grid. Every card draws its own structure, and no
   card shows a blank or a leftover drawing from another row.
4. Hover different rows in the grid and watch the mouse-over card: it draws the
   structure of the row under the pointer.

## Renderer sizes stay readable across the whole ladder (exploratory)

Setup: as above, on **spgi-100**, with **Structure** among the fields.

1. Set **Renderer Size** to **large**. The structure drawings are noticeably
   bigger and show more detail — individual atom labels and bond types are easy
   to read.
2. Set **Renderer Size** to **normal**, then **small**. At each step the drawing
   shrinks, but the molecule is still recognizable at small: it degrades
   gracefully rather than turning into an unreadable smudge.
3. At each of the three sizes, read the column name labels in the left header
   pane: every label stays fully readable and stays aligned with the field it
   names — no clipping, no overlap with the field above or below, no truncation
   that hides which field is which.

## Curves draw with the right colours, line styles and axes (exploratory)

Setup: Close all, open **System:DemoFiles/curves.csv**, add **Forms**.

1. Set the field set to **smiles**, **multiple prefit**. Each card shows a
   dose-response chart, not raw text: data points, fitted curve, and axes are all
   present.
2. Set the field set to **smiles**, **multiple styled series**, **styled
   proportional with IC50**. Both curve columns draw as separate charts, and each
   series inside a chart uses its declared colour and line style (solid / dashed
   as configured in the data) — the styling is not collapsed to one default
   colour.
3. Set **Renderer Size** to **large**: the charts get bigger and the axes carry
   more tick labels. Set it to **small**: the charts remain visible and the axes
   simplify — fewer ticks and shorter labels — rather than becoming an unreadable
   tangle.
4. Ctrl+click three rows. Each card draws its own molecule and its own curve, and
   no card reuses the neighbouring row's chart.


## Substructure highlighting and alignment match between grid and form (exploratory)

Setup: Close all, open **System:AppData/Chem/tests/spgi-100.csv**, add **Forms**,
set **Renderer Size** to **normal** and keep **Structure** among the fields.

1. Ctrl+click four rows in the grid so their cards are visible.
2. Open the filter panel, find the **Structure** column's molecular filter, and
   sketch a substructure that matches part of the dataset.
3. Look at a surviving row in the **grid**: its structure cell now draws the
   matched fragment highlighted, and the molecule is aligned on that fragment.
4. Look at the **same row's card** in the Forms viewer: the field shows the same
   highlight on the same fragment and the same alignment — the form and the grid
   agree about what was matched, allowing for the different drawing size.
5. Change the substructure query. Both the grid cell and the form field follow
   the new query together.

---
{
  "order": 900,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv", "System:DemoFiles/curves.csv"]
}
