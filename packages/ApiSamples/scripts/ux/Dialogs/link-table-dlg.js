//Example "Link table" dialog layout

let table = grok.data.testData('demog', 10000);
grok.shell.addTableView(table);
let tables = grok.shell.tables;

let linkTablesDlg = ui.dialog('Link Tables');
let linkType = ui.choiceInput('Link Type', 'row to filter', ['row to filter', 'row to row', 'filter to filter']);
let inputs = ui.inputs([
  ui.divH([ui.tableInput('Tables', tables[0], tables), ui.tableInput('', tables[0], tables)]),
  ui.divH([ui.columnInput('Key Columns', tables[0]), ui.columnInput('', tables[0]), ui.button(ui.iconFA('plus', ()=>{
    let newInputRow = ui.divH([ui.columnInput(' ', tables[0]), ui.columnInput('', tables[0]), ui.button(ui.iconFA('trash-alt'),()=>{
      $(newInputRow).detach();
    })])
    inputs.insertBefore(newInputRow, linkType.root);
    $(inputs).find('input, select, .d4-column-selector').css('min-width','120px');
  }))]),
  linkType
])

linkTablesDlg.add(inputs);
linkTablesDlg.show({x:50, y:100});
$(inputs).find('input, select, .d4-column-selector').css('min-width','120px');