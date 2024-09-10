import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


function showOutputs(t1: { value: string; }, t2: { value: string; },
  c1: DG.InputBase, c2: DG.InputBase, outputs: { innerHTML: string; append: (arg0: HTMLDivElement) => void; }) {
  const selection = grok.shell.table(t1.value).selection;
  const col1 = grok.shell.table(t1.value).columns.byName(c1.value);
  const col2 = grok.shell.table(t2.value).columns.byName(c2.value);

  outputs.innerHTML = '';

  if (col1.length === col2.length) {
    let counter = 0;
    for (let i = 0; i < col1.length; i++) {
      if (col1.getString(i) == col2.getString(i))
        counter++;
    }

    const selectButton = ui.link(
      'Select',
      () => selection.init((i) => col1.getString(i) != col2.getString(i)));

    const addColumnButton = ui.button('Add Column', () => {
      grok.shell.table(t1.value).columns.addNewCalculated(
        c1.value + ' == ' + c2.value,
        '${' + col1.name + '} == ${' + col2.name + '}',
        DG.TYPE.BOOL,
        false,
      ).then((_: any) => {if (t1.value != t2.value) grok.shell.info('Column was added to ' + t1.value);});
    });

    const matched = ui.divText('Matched: ' + counter.toString(), {});
    const mismatched = ui.divH([ui.divText('Mismatched: ' + (selection.length - counter).toString()), selectButton]);

    mismatched.style.marginLeft = '152px';
    matched.style.marginTop = '10px';
    matched.style.marginLeft = '152px';
    selectButton.style.marginLeft = '4px';
    outputs.append(
      ui.divV([
        matched,
        mismatched,
        addColumnButton,
      ]),
    );
  } else {
    const d = ui.divV([
      ui.divText('Length mismatch: '),
      ui.divText(t1.value + ': ' + col1.length + ' rows, '),
      ui.divText(t2.value + ': ' + col2.length + ' rows'),
    ]);
    d.style.marginTop = '10px';
    d.style.marginLeft = '152px';
    d.style.color = 'var(--red-3)';
    outputs.append(d);
  }
}

function addCols(t1: DG.InputBase, t2: DG.InputBase, outputs: HTMLDivElement) {
  let firstColumnAdded = false;
  let secondColumnAdded = false;

  const c1 = ui.input.choice('Columns', {value: '', items: grok.shell.table(t1.value).columns.names(), onValueChanged: (value) => {
    firstColumnAdded = true;
    if (secondColumnAdded)
      showOutputs(t1, t2, c1, c2, outputs);
    else if (t1.value === t2.value && grok.shell.table(t1.value).columns.length === 2) {
      const columnsNames = grok.shell.table(t1.value).columns.names();
      c2.value = (columnsNames.indexOf(value) === 0) ? columnsNames[1] : columnsNames[0];
    }
  }});

  const c2 = ui.input.choice('', {value: '', items: grok.shell.table(t2.value).columns.names(), onValueChanged: (value) => {
    secondColumnAdded = true;
    if (firstColumnAdded)
      showOutputs(t1, t2, c1, c2, outputs);
    else if (t1.value === t2.value && grok.shell.table(t2.value).columns.length === 2) {
      const columnsNames = grok.shell.table(t2.value).columns.names();
      c1.value = (columnsNames.indexOf(value) === 0) ? columnsNames[1] : columnsNames[0];
    }
  }});

  return ui.divH([c1.root, c2.root]);
}

//name: Compare Columns
//tags: autostart
export function compareColumns() {
  const tablesNames = grok.shell.tables.map((t) => t.name);

  let firstTableAdded = false;
  let secondTableAdded = false;

  const t1 = ui.input.choice('Tables', {value: '', items: tablesNames, onValueChanged: (value) => {
    firstTableAdded = true;
    if (secondTableAdded) {
      columnsInputs.innerHTML = '';
      columnsInputs.append(addCols(t1, t2, outputs));
    } else if (tablesNames.length === 1)
      t2.value = value;
  }});

  const t2 = ui.input.choice('', {value: '', items: tablesNames, onValueChanged: (value) => {
    secondTableAdded = true;
    if (firstTableAdded) {
      columnsInputs.innerHTML = '';
      columnsInputs.append(addCols(t1, t2, outputs));
    } else if (tablesNames.length === 1)
      t1.value = value;
  }});

  const columnsInputs = ui.div([]);
  const outputs = ui.div([]);

  const inputSection = ui.div([
    ui.divH([t1.root, t2.root], {}),
    columnsInputs,
    outputs,
  ]);

  ui.dialog('Compare Columns')
    .add(inputSection)
    .show();
}
