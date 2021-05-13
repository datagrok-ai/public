import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function compare(c1, c2) {
  if (c1.name === c2.name) {c2.name = c1.stringValue + ' (1)'};
  let df = DG.DataFrame.fromColumns([c1, c2]);
  df.columns.addNewCalculated(
    c1.name + '==' + c2.name,
    '${' + c1.name + '}==${' + c2.name + '}'
  ).then((_) => {
    grok.shell.addTableView(df);
  });
}

function addCols(t1, t2, container) {

  let firstColumnAdded = false,
    secondColumnAdded = false;

  let c1 = ui.choiceInput('', '', grok.shell.tables.find(({name}) => name === t1.value).columns.names(), (chosenColumnName) => {
    firstColumnAdded = true;
    if (secondColumnAdded)
      compare(
        grok.shell.tables.find(({name}) => name === t1.value).columns.byName(chosenColumnName),
        grok.shell.tables.find(({name}) => name === t2.value).columns.byName(c2.value)
      );
  });

  let c2 = ui.choiceInput('', '', grok.shell.tables.find(({name}) => name === t2.value).columns.names(), (chosenColumnName) => {
    secondColumnAdded = true;
    if (firstColumnAdded)
      compare(
        grok.shell.tables.find(({name}) => name === t1.value).columns.byName(c1.value),
        grok.shell.tables.find(({name}) => name === t2.value).columns.byName(chosenColumnName)
      );
  });

  c1.input.style.width = '150px';
  c2.input.style.width = '150px';

  return ui.div([ui.label('Columns'), ui.divH([c1, c2])]);
}

//name: Compare Columns
//tags: autostart
export function compareColumns() {

  let tables = grok.shell.tables;
  let tablesNames = tables.map((t) => t.name);

  let firstTableAdded = false,
    secondTableAdded = false;

  let t1 = ui.choiceInput('', '', tablesNames, () => {
    firstTableAdded = true;
    if (secondTableAdded) {
      container.innerHTML = '';
      container.append(addCols(t1, t2, container));
    }
  });

  let t2 = ui.choiceInput('', '', tablesNames, () => {
    secondTableAdded = true;
    if (firstTableAdded) {
      container.innerHTML = '';
      container.append(addCols(t1, t2, container));
    }
  });

  t1.input.style.width = '150px';
  t2.input.style.width = '150px';

  let container = ui.div();
  ui.dialog('Compare Columns')
    .add(
      ui.divV([
        ui.label('Tables'),
        ui.divH([t1, t2], {style: {marginBottom: '10px'}}),
        container
      ])
    )
    .show();
}