import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function showOutputs(t1, t2, c1, c2, outputs) {

  let selection = grok.shell.tables.find(({name}) => name === t1.value).selection;
  let col1 = grok.shell.tables.find(({name}) => name === t1.value).columns.byName(c1.value);
  let col2 = grok.shell.tables.find(({name}) => name === t2.value).columns.byName(c2.value);

  outputs.innerHTML = '';

  if (col1.length === col2.length) {

    for (let i = 0; i < col1.length; i++) if (col1.get(i) != col2.get(i)) selection.set(i, true, false);

    let link = ui.element('a');
    link.text = ' Select';
    link.onclick = (e) => {selection.fireChanged();};

    outputs.append(
      ui.buttonsInput([
        ui.divText('Matched: ' + selection.trueCount),
        ui.span(['Mismatched: ' + (selection.length - selection.trueCount), link])
      ])
    );

  } else {

    outputs.append(
      ui.buttonsInput([
        ui.divText('Length mismatch: '),
        ui.divText(col1.name + ': ' + col1.length + ' rows, '),
        ui.divText(col2.name + ': ' + col2.length + ' rows')
      ])
    );
  }
}

function addCols(t1, t2, outputs) {

  let firstColumnAdded = false,
    secondColumnAdded = false;

  let c1 = ui.choiceInput('Columns', '', grok.shell.tables.find(({name}) => name === t1.value).columns.names(), () => {
    firstColumnAdded = true;
    if (secondColumnAdded) showOutputs(t1, t2, c1, c2, outputs);
  });

  let c2 = ui.choiceInput('', '', grok.shell.tables.find(({name}) => name === t2.value).columns.names(), () => {
    secondColumnAdded = true;
    if (firstColumnAdded) showOutputs(t1, t2, c1, c2, outputs);
  });

  return ui.divH([c1, c2]);
}

//name: Compare Columns
//tags: autostart
export function compareColumns() {

  let tablesNames = grok.shell.tables.map((t) => t.name);

  let firstTableAdded = false,
    secondTableAdded = false;

  let t1 = ui.choiceInput('Tables', '', tablesNames, () => {
    firstTableAdded = true;
    if (secondTableAdded) {
      columnsInputs.innerHTML = '';
      columnsInputs.append(addCols(t1, t2, outputs));
    }
  });

  let t2 = ui.choiceInput('', '', tablesNames, () => {
    secondTableAdded = true;
    if (firstTableAdded) {
      columnsInputs.innerHTML = '';
      columnsInputs.append(addCols(t1, t2, outputs));
    }
  });

  let columnsInputs = ui.div();
  let outputs = ui.div();

  let inputSection = ui.form([
    ui.divH([t1, t2]),
    columnsInputs,
    outputs
  ], 'ui-form-aligned');

  $(inputSection).css('ui-form-aligned');

  ui.dialog('Compare Columns')
    .add(inputSection)
    .show();
}