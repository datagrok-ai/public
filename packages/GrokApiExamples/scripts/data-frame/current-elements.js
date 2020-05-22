// Selecting current row, column and cell

let t = grok.data.testData('demog', 100);
grok.shell.addTableView(t);

t.currentRow = 4;
t.currentCol = t.col('age');
t.currentCell = t.cell(5, 'sex');

t.onCurrentCellChanged.subscribe((_) => grok.shell.info(`Cell: ${t.currentCell.rowIndex}, ${t.currentCell.column.name}`));
t.onCurrentColChanged.subscribe((_) => grok.shell.info(`Column: ${t.currentCol.name}`));
t.onCurrentRowChanged.subscribe((_) => grok.shell.info(`Row: ${t.currentRow.idx}`));
