// Modifying grid column styles

let view = grok.shell.addTableView(grok.data.demo.demog());
let age = view.grid.columns.byName('age');

age.contentCellStyle.font = '10px Roboto, Roboto Local';
age.contentCellStyle.backColor = DG.Color.red;
age.contentCellStyle.horzAlign = 'center';

age.headerCellStyle.backColor = DG.Color.green;