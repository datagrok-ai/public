let strCol = DG.Column.fromList('string', 'string', ['>5', '4', '<3']);
let qCol = strCol.convertTo(DG.COLUMN_TYPE.QNUM);
qCol.name = 'qnum';

let t = DG.DataFrame.fromColumns([strCol, qCol]);
grok.shell.addTableView(t);