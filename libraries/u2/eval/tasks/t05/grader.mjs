import * as grok from 'datagrok-api/grok';
import {Func, FuncParam} from '../../../tests/platform-doubles.mjs';
import {buttons, choose, clickOk, errors, expect, flush, input, need, offered, type} from '../../lib/harness.mjs';
import {frame} from '../../lib/platform.mjs';

function openTables(names) {
  grok.shell.dart.tables = [];
  for (const name of names) {
    const table = grok.shell.addTable(frame([{name: 'smiles', type: 'string', semType: 'Molecule'}], [],
      name ?? 'renamed'));
    if (name === null)
      table.name = undefined;
  }
}

async function open(mod, tables) {
  openTables(tables);
  const runs = [];
  const func = new Func('descriptors', {friendlyName: 'Descriptors', inputs: [new FuncParam('table', 'dataframe'),
    new FuncParam('threshold', 'double', {defaultValue: 0.5, isOptional: true})]});
  const dialog = need(mod, 'settingsDialog')(func, (call) => void runs.push(call));
  await flush();
  const root = dialog?.root ?? dialog;
  expect(root instanceof HTMLElement && root.isConnected, 'settingsDialog() did not show a dialog');
  return {root, runs};
}

export const checks = {
  'form-runs': async (mod) => {
    const {root, runs} = await open(mod, ['demog', 'cars']);
    await choose(input('table', root), 'cars');
    await type(input('threshold', root), '0.8');
    const run = buttons(root, /^(run|ok)$/)[0];
    expect(run !== undefined, 'no RUN button');
    await clickOk(run, 'RUN');
    expect(runs.length === 1, `onRun ran ${runs.length} times`);
    const values = runs[0].inputs;
    expect((values.table?.name ?? values.table) === 'cars', `the call got table ${values.table?.name ?? values.table}`);
    expect(values.threshold === 0.8, `the call got threshold ${values.threshold}`);
  },
  'nameless-table-open': async (mod) => {
    const {root} = await open(mod, ['demog', null, 'cars']);
    expect(errors().length === 0, `building the form logged: ${errors()[0]}`);
    const tables = offered(input('table', root)).map((i) => i.label);
    expect(tables.length === 2 && tables.includes('demog') && tables.includes('cars'),
      `with a table open that has no name, the table choice offers ${JSON.stringify(tables)}`);
  },
};
