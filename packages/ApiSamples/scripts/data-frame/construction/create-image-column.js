// A byte array column holds binary values; PowerGrid draws PNG, JPEG, GIF, WebP and BMP bytes as pictures.
// Double-click an image to preview or replace it, or an empty cell to pick an image file.

const photo = await grok.dapi.files.readAsBytes('System:DemoFiles/images/cats/cat1.jpg');
const t = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('name', ['cat', 'empty']),
  DG.Column.fromList(DG.COLUMN_TYPE.BYTE_ARRAY, 'image', [photo, null]),
]);
grok.shell.addTableView(t);
