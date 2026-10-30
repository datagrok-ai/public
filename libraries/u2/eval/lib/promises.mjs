// node --import ./register.mjs promises.mjs: looks up, in the page, every platform API AGENT-PROMPT.md promises;
// prints the ones missing as JSON (`[]` when the page keeps them all).
import {readFileSync} from 'node:fs';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import './harness.mjs';
import {frame, tableView} from './platform.mjs';

const fn = (x) => typeof x === 'function';
const stream = (x) => fn(x?.subscribe);
const has = (obj, names) => names.every((n) => fn(obj?.[n]));

/** Every `on…` event js-api's DataFrame declares, so the page follows js-api as it grows. */
const declared = readFileSync(new URL('../../node_modules/datagrok-api/dist/src/dataframe/data-frame.d.ts',
  import.meta.url), 'utf8');
const EVENTS = [...new Set(/export declare class DataFrame \{[\s\S]*?\n\}/.exec(declared)[0]
  .match(/get on[A-Z]\w*(?=\(\))/g).map((g) => g.slice(4)))];

const df = frame([{name: 'Molecule', type: 'string', semType: 'Molecule'}, {name: 'ID', type: 'string'}],
  [{Molecule: 'CCO', ID: 'V1'}], 'hits');
const tv = tableView(df);
const view = DG.View.fromRoot(document.createElement('div'));
const func = new DG.Func('f', {friendlyName: 'F', inputs: [new DG.Property('x', 'double', {defaultValue: 1,
  caption: 'X'})]});

const PROMISED = {
  'DG.DataFrame name, rowCount': () => df instanceof DG.DataFrame && df.name === 'hits' && df.rowCount === 1,
  'DataFrame.columns toList, names, byName, bySemType': () => has(df.columns, ['toList', 'names', 'byName', 'bySemType']) &&
    df.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name === 'Molecule',
  'DataFrame.col(name)': () => df.col('ID')?.name === 'ID',
  'DataFrame.rows.addNew(values, notify)': () => {
    df.rows.addNew(['CCN', 'V2'], true);
    return df.rowCount === 2 && df.get('ID', 1) === 'V2';
  },
  ...Object.fromEntries(EVENTS.map((e) => [`DataFrame.${e}`, () => stream(df[e])])),
  'DG.Column name, type, semType, get(i)': () => {
    const c = df.columns.byName('Molecule');
    return c instanceof DG.Column && c.name === 'Molecule' && c.type === 'string' && c.semType === 'Molecule' &&
      c.get(0) === 'CCO';
  },
  'DG.FileInfo fullPath, fileName, fromBytes': () => {
    const f = DG.FileInfo.fromBytes('hits.csv', new Uint8Array());
    return f instanceof DG.FileInfo && f.fileName === 'hits.csv' && 'fullPath' in f;
  },
  'DG.Func name, friendlyName, inputs, prepare': () => func.name === 'f' && func.friendlyName === 'F' &&
    func.inputs.length === 1 && fn(func.prepare),
  'DG.FuncCall inputs, setParamValue': () => {
    const call = func.prepare({x: 2});
    call.setParamValue('x', 3);
    return call instanceof DG.FuncCall && call.inputs.x === 3;
  },
  'DG.Property name, propertyType, defaultValue, caption': () => {
    const p = func.inputs[0];
    return p.name === 'x' && p.propertyType === 'double' && p.defaultValue === 1 && p.caption === 'X';
  },
  'DG.TableView dataFrame, getRibbonPanels, setRibbonPanels, close': () => tv instanceof DG.TableView &&
    tv instanceof DG.View && tv instanceof DG.ViewBase && tv.dataFrame === df &&
    has(tv, ['getRibbonPanels', 'setRibbonPanels', 'close']),
  'DG.View': () => view instanceof DG.View && view instanceof DG.ViewBase && !(view instanceof DG.TableView),
  'DG.TYPE, DG.SEMTYPE': () => DG.TYPE.DATA_FRAME === 'dataframe' && DG.SEMTYPE.MOLECULE === 'Molecule',
  'grok.shell v, views, tables, tableNames, table, addTable, info, warning, error': () => grok.shell.v === tv &&
    grok.shell.views.includes(tv) && Array.isArray(grok.shell.tables) && Array.isArray(grok.shell.tableNames) &&
    has(grok.shell, ['table', 'addTable', 'info', 'warning', 'error']),
  'grok.events onCurrentViewChanged, onViewAdded, onViewRemoved, onTableAdded, onTableRemoved': () =>
    ['onCurrentViewChanged', 'onViewAdded', 'onViewRemoved', 'onTableAdded', 'onTableRemoved']
      .every((e) => stream(grok.events[e])),
  'grok.dapi.files exists, list, readAsText': () => has(grok.dapi.files, ['exists', 'list', 'readAsText']),
  'grok.chem.drawMolecule': () => grok.chem.drawMolecule('CCO', 160, 100) instanceof HTMLElement,
  'ui.pickFileFromFiles': () => fn(ui.pickFileFromFiles),
};

const missing = [];
for (const [name, present] of Object.entries(PROMISED)) {
  try {
    if (!present())
      missing.push(name);
  }
  catch (e) {
    missing.push(`${name} (${e instanceof Error ? e.message : e})`);
  }
}
process.stdout.write(JSON.stringify(missing) + '\n');
process.exit(0);
