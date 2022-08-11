import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


function showOutputs(t1: { value: string; }, t2: { value: string; }, c1: DG.InputBase, c2: DG.InputBase, outputs: { innerHTML: string; append: (arg0: HTMLDivElement) => void; }) {

  let selection = grok.shell.table(t1.value).selection;
  let col1 = grok.shell.table(t1.value).columns.byName(c1.value);
  let col2 = grok.shell.table(t2.value).columns.byName(c2.value);

  outputs.innerHTML = '';

  if (col1.length === col2.length) {

    let counter = 0;
    for (let i = 0; i < col1.length; i++)
      if (col1.getString(i) == col2.getString(i))
        counter++;

    let selectButton = ui.link(
      'Select',
      () => selection.init((i) => col1.getString(i) != col2.getString(i)));

    let addColumnButton = ui.button('Add Column', () => {
      grok.shell.table(t1.value).columns.addNewCalculated(
        c1.value + ' == ' + c2.value,
        '${' + col1.name + '} == ${' + col2.name + '}',
        DG.TYPE.BOOL,
        false
      ).then((_: any) => {if (t1.value != t2.value) grok.shell.info('Column was added to ' + t1.value)});
    });

    let matched = ui.divText('Matched: ' + counter.toString(), {});
    let mismatched = ui.divH([ui.divText('Mismatched: ' + (selection.length - counter).toString()), selectButton]);

    mismatched.style.marginLeft = '152px';
    matched.style.marginTop = '10px';
    matched.style.marginLeft = '152px';
    selectButton.style.marginLeft = '4px';
    outputs.append(
      ui.divV([
        matched,
        mismatched,
        addColumnButton
      ])
    );

  } else {

    let d = ui.divV([
      ui.divText('Length mismatch: '),
      ui.divText(t1.value + ': ' + col1.length + ' rows, '),
      ui.divText(t2.value + ': ' + col2.length + ' rows')
    ]);
    d.style.marginTop = '10px';
    d.style.marginLeft = '152px';
    d.style.color = 'var(--red-3)';
    outputs.append(d);
  }
}

function addCols(t1: DG.InputBase, t2: DG.InputBase, outputs: HTMLDivElement) {

  let firstColumnAdded = false,
    secondColumnAdded = false;

  let c1 = ui.choiceInput('Columns', '', grok.shell.table(t1.value).columns.names(), () => {
    firstColumnAdded = true;
    if (secondColumnAdded) {
      showOutputs(t1, t2, c1, c2, outputs);
    } else if (t1.value === t2.value && grok.shell.table(t1.value).columns.length === 2) {
      let columnsNames = grok.shell.table(t1.value).columns.names();
      c2.value = (columnsNames.indexOf(c1.value!) === 0) ? columnsNames[1] : columnsNames[0];
    }
  });

  let c2 = ui.choiceInput('', '', grok.shell.table(t2.value).columns.names(), () => {
    secondColumnAdded = true;
    if (firstColumnAdded) {
      showOutputs(t1, t2, c1, c2, outputs);
    } else if (t1.value === t2.value && grok.shell.table(t2.value).columns.length === 2) {
      let columnsNames = grok.shell.table(t2.value).columns.names();
      c1.value = (columnsNames.indexOf(c2.value!) === 0) ? columnsNames[1] : columnsNames[0];
    }
  });

  return ui.divH([c1.root, c2.root]);
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
    } else if (tablesNames.length === 1) {
      t2.value = t1.value;
    }
  });

  let t2 = ui.choiceInput('', '', tablesNames, () => {
    secondTableAdded = true;
    if (firstTableAdded) {
      columnsInputs.innerHTML = '';
      columnsInputs.append(addCols(t1, t2, outputs));
    } else if (tablesNames.length === 1) {
      t1.value = t2.value;
    }
  });

  let columnsInputs = ui.div([]);
  let outputs = ui.div([]);

  let inputSection = ui.div([
    ui.divH([t1.root, t2.root], {}),
    columnsInputs,
    outputs
  ]);

  ui.dialog('Compare Columns')
    .add(inputSection)
    .show();
}