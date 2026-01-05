// Note that bytearray column is empty by default and cant be created so far with Column.fromList method
// double click on cell to insert image

let t = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('images', ['im1', 'im2', 'im3']),
  DG.Column.fromType(DG.COLUMN_TYPE.BYTE_ARRAY, 'byte array col', 3),
]);
grok.shell.addTableView(t);