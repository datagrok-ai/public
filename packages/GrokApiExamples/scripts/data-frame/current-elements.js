// Selecting current row, column and cell

let t = grok.testData('demog', 100);
grok.addTableView(t);

t.currentRow = 4;
t.currentCol = t.col('age');
t.currentCell = t.cell(5, 'sex');

t.onCurrentCellChanged(() => grok.balloon.info(`Row: ${t.currentRow.idx}`));
t.onCurrentColChanged(() => grok.balloon.info(`Column: ${t.currentCol.name}`));
t.onCurrentRowChanged(() => grok.balloon.info(`Cell: ${t.currentCell.rowIndex}, ${t.currentCell.column.name}`));
