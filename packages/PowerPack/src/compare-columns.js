import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function compare(t1, c1, c2) {
  t1.rows.select((row) => row[c1.name] === row[c2.name]);
  grok.shell.add(t1);
}

function addCols(t1, t2) {

  let firstColumnAdded = false,
    secondColumnAdded = false;

  let c1 = ui.choiceInput('Columns', '', grok.shell.tables.find(({name}) => name === t1.value).columns.names(), (chosenColumnName) => {
    firstColumnAdded = true;
    if (secondColumnAdded)
      compare(
        grok.shell.tables.find(({name}) => name === t1.value),
        grok.shell.tables.find(({name}) => name === t1.value).columns.byName(chosenColumnName),
        grok.shell.tables.find(({name}) => name === t2.value).columns.byName(c2.value)
      );
  });

  let c2 = ui.choiceInput('', '', grok.shell.tables.find(({name}) => name === t2.value).columns.names(), (chosenColumnName) => {
    secondColumnAdded = true;
    if (firstColumnAdded)
      compare(
        grok.shell.tables.find(({name}) => name === t1.value),
        grok.shell.tables.find(({name}) => name === t1.value).columns.byName(c1.value),
        grok.shell.tables.find(({name}) => name === t2.value).columns.byName(chosenColumnName)
      );
  });

  return ui.divH([c1, c2]);
}

//name: Compare Columns
//tags: autostart
export function compareColumns() {

  let tables = grok.shell.tables;
  let tablesNames = tables.map((t) => t.name);

  let firstTableAdded = false,
    secondTableAdded = false;

  let t1 = ui.choiceInput('Tables', '', tablesNames, () => {
    firstTableAdded = true;
    if (secondTableAdded) {
      container.innerHTML = '';
      container.append(addCols(t1, t2));
    }
  });

  let t2 = ui.choiceInput('', '', tablesNames, () => {
    secondTableAdded = true;
    if (firstTableAdded) {
      container.innerHTML = '';
      container.append(addCols(t1, t2));
    }
  });

  let container = ui.div();
  let inputSection = ui.form([
    ui.divH([t1, t2], {style: {marginBottom: '10px'}}),
    container
  ], 'ui-form-aligned');

  $(inputSection).css('ui-form-aligned');

  ui.dialog('Compare Columns')
    .add(inputSection)
    .show();
}