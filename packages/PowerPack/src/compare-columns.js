import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//name: Compare Columns
//tags: autostart
export function compareColumns() {
  let tables = grok.shell.tables;
  let tablesNames = tables.map((t) => t.name);

  let t1 = ui.choiceInput('', '', tablesNames, () => {});
  let t2 = ui.choiceInput('', '', tablesNames, () => {});
  let c1 = ui.choiceInput('', '', tables[0].columns.names(), () => {});
  let c2 = ui.choiceInput('', '', tables[0].columns.names(), () => {});

  t1.input.style.width = '150px';
  t2.input.style.width = '150px';
  c1.input.style.width = '150px';
  c2.input.style.width = '150px';

  ui.dialog('Compare Columns')
    .add(
      ui.divV([
        ui.label('Tables'),
        ui.divH([t1, t2], {style: {marginBottom: '10px'}}),
        ui.label('Columns'),
        ui.divH([c1, c2])
      ])
    )
    .add(ui.button('Compare', () => {
      let tbls = grok.shell.tables[0];
      let a = tbls.columns.byName(c1.value);
      let b = tbls.columns.byName(c2.value);
      tbls.columns.addNewBool('Result').init((i) => (a.get(i) === b.get(i)));
    }))
    .show();
}