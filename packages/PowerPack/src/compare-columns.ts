import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


function showOutputs(t1: { value: string; }, t2: { value: string; }, c1: DG.InputBase, c2: DG.InputBase, outputs: { innerHTML: string; append: (arg0: HTMLDivElement) => void; }) {

  // @ts-ignore
  let selection = grok.shell.tables.find(({name}) => name === t1.value).selection;
  // @ts-ignore
  let col1 = grok.shell.tables.find(({name}) => name === t1.value).columns.byName(c1.value);
  // @ts-ignore
  let col2 = grok.shell.tables.find(({name}) => name === t2.value).columns.byName(c2.value);

  outputs.innerHTML = '';

  if (col1.length === col2.length) {

    let counter = 0;
    for (let i = 0; i < col1.length; i++)
      if (col1.getString(i) == col2.getString(i))
        counter++;

    let link = ui.element('a');
    link.textContent = ' Select';
    link.onclick = (ev) => {selection.init((i) => col1.get(i) != col2.get(i));};

    let str1: string = 'Mismatched: ' + (selection.length - counter).toString();
    // @ts-ignore
    outputs.append(
      ui.div([
        ui.divText('Matched: ' + counter.toString(), {}),
        ui.divText(str1 + link)
        //ui.buttonsInput([ui.span([{[str1]: link}])])
      ])
    );

  } else {

    outputs.append(
        ui.divV([
          ui.divText('Length mismatch: '),
          ui.divText(t1.value + ': ' + col1.length + ' rows, '),
          ui.divText(t2.value + ': ' + col2.length + ' rows')
        ])
    );
  }
}

function addCols(t1: DG.InputBase, t2: DG.InputBase, outputs: HTMLDivElement) {

  let firstColumnAdded = false,
      secondColumnAdded = false;

  // @ts-ignore
  let c1 = ui.choiceInput('Columns', '', grok.shell.tables.find(({name}) => name === t1.value).columns.names(), () => {
    firstColumnAdded = true;
    let df = grok.shell.tables.find(({name}) => name === t1.value);
    if (secondColumnAdded) {
      showOutputs(t1, t2, c1, c2, outputs);
    } else { // @ts-ignore
      if (t1.value === t2.value && df.columns.length === 2) {
            // @ts-ignore
            let columnsNames = df.columns.names();
            c2.value = (columnsNames.indexOf(c1.value) === 0) ? columnsNames[1] : columnsNames[0];
          }
    }
  });

  // @ts-ignore
  let c2 = ui.choiceInput('', '', grok.shell.tables.find(({name}) => name === t2.value).columns.names(), () => {
    secondColumnAdded = true;
    let df = grok.shell.tables.find(({name}) => name === t2.value);
    if (firstColumnAdded) {
      showOutputs(t1, t2, c1, c2, outputs);
    } else { // @ts-ignore
      if (t1.value === t2.value && df.columns.length === 2) {
            // @ts-ignore
            let columnsNames = df.columns.names();
            c1.value = (columnsNames.indexOf(c2.value) === 0) ? columnsNames[1] : columnsNames[0];
          }
    }
  });

  // @ts-ignore
  return ui.divH([c1, c2]);
}

//name: Compare Columns
//tags: autostart
export function compareColumns() {

  grok.shell.info('hello, world! true');

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
  ], 'ui-form-aligned');

  //inputSection.css('ui-form-aligned');

  ui.dialog('Compare Columns')
      .add(inputSection)
      .show();
}