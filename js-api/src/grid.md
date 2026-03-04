# Grid API Usage Guide

Source: `grid.ts`. Full type reference: see that file or the generated `.d.ts`.
Samples: `../../packages/ApiSamples/scripts/grid`

## Common patterns

```js
const grid = grok.shell.tv.grid;   // accessing the grid in the currently visible table view
```

## Columns (`GridColumnList` / `GridColumn`)

```js
const col = grid.columns.byName('age');  // by name
const col = grid.col('age');             // shortcut
grid.columns.byIndex(1);                 // by index (0 = row header)
grid.columns.rowHeader;                  // row header column
```

Common `GridColumn` properties:

```js
col.selected = true;          // select column (visual highlight); call grid.invalidate() after
col.visible = false;          // hide column
col.width = 150;              // width in pixels
col.editable = false;         // make read-only
col.format = 'date';          // display format
col.backColor = 0xFFFF0000;   // background color (ARGB)
col.column;                   // underlying DataFrame Column (or null for virtual cols)

col.pin();                    // pin column to left
col.unpin();
col.move(2);                  // reorder to position index
col.scrollIntoView();         // scroll grid to show column
col.getVisibleCells();        // Iterable<GridCell> of visible cells
```

Bulk operations:

```js
grid.columns.setOrder(['sex', 'age', 'height']);   // reorder columns
grid.columns.setVisible(['sex', 'age']);            // show only these (hides rest)
grid.columns.add({ cellType: 'html', gridColumnName: 'custom' });  // add virtual column
```

**Selecting columns** — must use grid API; `t.columns.byName().selected` has no visual effect:

```js
grid.columns.byName('age').selected = true;
grid.columns.byName('sex').selected = true;
grid.invalidate();
```

## Sorting

```js
grid.sort(['age', 'sex'], [true, false]);  // ascending age, descending sex
grid.sort(['age']);                         // ascending (default)
grid.sortIndexes((a, b) => ...);           // custom row comparer
grid.setRowOrder(indexes);                 // explicit row order
grid.getRowOrder();                        // Int32Array of current order
```

## Scrolling

```js
grid.scrollToCell('age', 0);   // scroll to column 'age', grid row 0
grid.scrollToPixels(x, y);
grid.vertScroll;               // RangeSlider
grid.horzScroll;               // RangeSlider
```

## Cells (`GridCell`)

```js
const cell = grid.cell('age', 3);  // GridCell at column 'age', grid row 3
cell.value;                         // cell value (may differ from table cell — see onCellPrepare)
cell.tableRowIndex;                 // corresponding table row index
cell.gridRow;                       // grid row index
cell.gridColumn;                    // GridColumn
cell.tableColumn;                   // DataFrame Column
cell.bounds;                        // Rect (position within grid canvas)
cell.documentBounds;                // Rect (relative to document)
cell.style;                         // GridCellStyle (font, color, alignment, etc.)
cell.customText;                    // override display text
cell.setValue(x);                   // set value and fire onCellValueEdited

grid.hitTest(x, y);                 // GridCell at canvas coords (or null)
grid.getVisibleCells();             // Iterable<GridCell> of all visible cells
```

## Repaint

```js
grid.invalidate();                  // force redraw after programmatic changes
grid.runPostponedComputations();    // recalculate layout (e.g., after changing column widths)
grid.autoSize(maxW, maxH);          // resize grid to fit content
```

## Events (rxjs Observable)

```js
grid.onCellRender.subscribe(args => { /* custom rendering */ args.preventDefault(); });
grid.onCellRendered.subscribe(args => ...);
grid.onCellClick.subscribe(cell => ...);
grid.onCellDoubleClick.subscribe(cell => ...);
grid.onCurrentCellChanged.subscribe(cell => ...);
grid.onCellValueEdited.subscribe(cell => ...);
grid.onCellMouseEnter.subscribe(cell => ...);
grid.onColumnResized.subscribe(_ => ...);
grid.onRowsSorted.subscribe(_ => ...);

// Custom cell prepare (set value/style before rendering):
grid.onCellPrepare(cell => { cell.style.backColor = ...; });

// Custom tooltip:
grid.onCellTooltip((cell, x, y) => { /* return true to suppress default */ });
```

## Cell style (`GridCellStyle`)

Accessible via `cell.style`, `gridColumn.contentCellStyle`, `gridColumn.headerCellStyle`:

```js
style.font = 'bold 12px Arial';
style.textColor = 0xFF000000;    // ARGB
style.backColor = 0xFFFFFFFF;
style.horzAlign = 'right';       // 'left' | 'center' | 'right'
style.vertAlign = 'top';         // 'top' | 'center' | 'bottom'
style.textWrap = 'ellipsis';     // 'auto' | 'none' | 'new line' | 'words' | 'characters' | 'ellipsis'
style.margin = 4;
```
