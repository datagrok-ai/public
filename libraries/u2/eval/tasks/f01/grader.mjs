import * as grok from 'datagrok-api/grok';
import {FileInfo} from '../../../tests/platform-doubles.mjs';
import {buttons, chooseFiles, clickOk, errors, expect, field, flush, input, inputs, isOff, liveScopes, localFile,
  named, need, offered, press, row, type, until, wait} from '../../lib/harness.mjs';
import {busySession, frame, tableView} from '../../lib/platform.mjs';

const CAMPAIGN = {name: 'KHL-3 hits', filePath: 'System:AppData/HitTriage/KHL-3.csv', status: 'In Progress'};
const SHARE_PATH = 'System:DemoFiles/hits-2.csv';
const mol = (name) => ({name, type: 'string', semType: 'Molecule'});

/** The file read back: a table a campaign saved, with the platform's canonical-SMILES cache column first. */
const incoming = () => frame([mol('~Molecule.canonicalSmiles'), {name: 'ID', type: 'string'}, mol('Molecule'),
  mol('Scaffold')], [{'~Molecule.canonicalSmiles': 'CCO', 'ID': 'V7', 'Molecule': 'OCC', 'Scaffold': 'C1CC1'}],
'hits-2');

async function install(mod, service = {}) {
  busySession();
  const table = frame([mol('Molecule'), {name: 'Stage', type: 'string'}], [{Molecule: 'CCO', Stage: 'Idea'}], 'KHL-3');
  const view = tableView(table);
  const calls = {rename: [], download: 0, submit: 0, read: [], merge: [], refresh: 0, picks: 0, table: null};
  const share = new FileInfo(SHARE_PATH);
  grok.dapi.files.exists = async (path) => path === SHARE_PATH;
  grok.dapi.files.list = async () => [share];
  globalThis.grok_UI_PickFileFromFiles = async () => {
    calls.picks++;
    return share;
  };
  need(mod, 'campaignActions')(view, {...CAMPAIGN}, {
    rename: async (name) => void calls.rename.push(name),
    download: () => void calls.download++,
    submit: async () => void calls.submit++,
    readTable: async (file) => {
      calls.read.push(file);
      await wait(20);
      return calls.table = incoming();
    },
    merge: async (file, read, column) => void calls.merge.push({file, table: read, column}),
    refresh: async () => void calls.refresh++,
    ...service,
  });
  await flush();
  return {table, view, calls};
}

function onRibbon(view, re, what) {
  const el = named(view.dart.panels, re);
  expect(el !== null, `no ${what} on the view's ribbon`);
  return el;
}

async function openDialog(o, re, what, title) {
  await clickOk(onRibbon(o.view, re, what), what);
  await flush();
  const dialog = document.querySelectorAll('[role="dialog"]').reverse().find((d) => title.test(d.textContent));
  expect(dialog !== undefined, `${what} opened no "${title.source}" dialog`);
  return dialog;
}

const submitDialog = (o) => openDialog(o, /submit/i, 'Submit...', /submit campaign/i);
const mergeDialog = (o) => openDialog(o, /merge/i, 'Merge with table', /merge with table/i);

function control(dialog, caption, re, what, except) {
  const el = named(row(caption, dialog).el, re, except) ?? named(dialog, re, except);
  expect(el !== null, `no control to ${what}`);
  return el;
}

async function startEditing(o) {
  const dialog = await submitDialog(o);
  await clickOk(control(dialog, 'Name', /edit|rename|pencil/i, 'edit the name'), 'Edit name');
  return dialog;
}

const FILES_BROWSER = /from (datagrok )?files|datagrok files|files browser|folder-tree/i;

async function fromFiles(o, dialog) {
  await clickOk(control(dialog, 'Source file', /\bfiles\b/i, 'pick from Datagrok Files'), 'the Files picker');
  expect(o.calls.picks === 1, 'the platform\'s Files browser was not opened');
}

async function fromComputer(dialog) {
  const local = control(dialog, 'Source file', /open|choose|computer|upload|local|browse|folder/i,
    'open a local file', FILES_BROWSER);
  await clickOk(local, 'Open file');
  await chooseFiles([localFile('hits-2.csv')]);
}

function mergeButton(dialog) {
  const merge = buttons(dialog, /^(merge|ok)$/)[0];
  expect(merge !== undefined, 'no MERGE button');
  return merge;
}

async function merged(o, dialog) {
  const merge = mergeButton(dialog);
  expect(await until(() => !isOff(merge), 2000), 'MERGE stays off once a file is chosen');
  await clickOk(merge, 'MERGE');
  expect(await until(() => o.calls.merge.length === 1), `service.merge ran ${o.calls.merge.length} times`);
  return o.calls.merge[0];
}

const streams = () => grok.events.onViewRemoved.count + grok.events.onCurrentViewChanged.count +
  grok.events.onTableAdded.count + grok.events.onTableRemoved.count;
const cancel = (dialog) => clickOk(buttons(dialog, /^cancel$/)[0] ?? null, 'CANCEL');

export const checks = {
  'ribbon': async (mod) => {
    const {view} = await install(mod);
    for (const [re, what] of [[/add new row/i, 'Add new row'], [/refresh/i, 'Refresh Campaign'],
      [/merge/i, 'Merge with table'], [/submit/i, 'Submit...']])
      onRibbon(view, re, what);
    view.close();
  },
  'submit-summary': async (mod) => {
    const o = await install(mod);
    const dialog = await submitDialog(o);
    for (const [label, value] of [['Name', CAMPAIGN.name], ['File Path', CAMPAIGN.filePath], ['Status', CAMPAIGN.status]]) {
      const r = row(label, dialog);
      expect(r.text === value, `${label} shows "${r.text}", not "${value}"`);
      expect(r.control === undefined || r.control.getState().readOnly === true, `${label} is not read-only`);
    }
    expect(errors().length === 0, `opening Submit logged: ${errors()[0]}`);
    o.view.close();
  },
  'submit-rename': async (mod) => {
    const o = await install(mod);
    const dialog = await startEditing(o);
    const name = input('Name', dialog);
    expect(name.getState().readOnly === false, 'Name is still read-only after Edit name');
    expect(document.activeElement === field(name), 'the cursor is not in Name');
    await type(name, 'KHL-3 final');
    await press(field(name), 'Enter');
    expect(o.calls.submit === 0, 'Enter in Name submitted the campaign');
    expect(dialog.isConnected, 'Enter in Name closed the dialog');
    expect(o.calls.rename.length === 1 && o.calls.rename[0] === 'KHL-3 final',
      `rename got ${JSON.stringify(o.calls.rename)}`);
    expect(name.getState().readOnly === true, 'Name stays editable after Enter');
    o.view.close();
  },
  'submit-escape': async (mod) => {
    const o = await install(mod);
    const dialog = await startEditing(o);
    const name = input('Name', dialog);
    await type(name, 'zzz');
    await press(field(name), 'Escape');
    expect(dialog.isConnected, 'Escape in Name closed the dialog');
    expect(o.calls.rename.length === 0, 'Escape renamed the campaign');
    const s = name.getState();
    expect(s.text === CAMPAIGN.name, `Name shows "${s.text}" after Escape`);
    expect(s.readOnly === true, 'Name stays editable after Escape');
    o.view.close();
  },
  'submit-download': async (mod) => {
    const o = await install(mod);
    const dialog = await submitDialog(o);
    await clickOk(control(dialog, 'File Path', /download/i, 'download the file'), 'Download');
    expect(o.calls.download === 1, `service.download() ran ${o.calls.download} times`);
    o.view.close();
  },
  'submit-runs': async (mod) => {
    const o = await install(mod);
    const dialog = await submitDialog(o);
    const submit = buttons(dialog, /^submit$/)[0];
    expect(submit !== undefined, 'no SUBMIT button');
    await clickOk(submit, 'SUBMIT');
    expect(await until(() => o.calls.submit === 1), `service.submit() ran ${o.calls.submit} times`);
    o.view.close();
  },
  'merge-needs-file': async (mod) => {
    const o = await install(mod);
    const merge = mergeButton(await mergeDialog(o));
    expect(isOff(merge), 'MERGE is on before a file is chosen');
    merge.click();
    await flush();
    expect(o.calls.merge.length === 0, 'MERGE merged with no file');
    o.view.close();
  },
  'merge-from-files': async (mod) => {
    const o = await install(mod);
    const dialog = await mergeDialog(o);
    await fromFiles(o, dialog);
    const m = await merged(o, dialog);
    expect(m.file?.fullPath === SHARE_PATH, `merged ${m.file?.fullPath ?? m.file}, not ${SHARE_PATH}`);
    expect(m.table === o.calls.table, 'MERGE did not pass the table service.readTable read');
    o.view.close();
  },
  'merge-typed-path': async (mod) => {
    const o = await install(mod);
    const dialog = await mergeDialog(o);
    const source = row('Source file', dialog);
    await type(source.control ?? {root: source.el}, SHARE_PATH);
    const m = await merged(o, dialog);
    expect(m.file?.fullPath === SHARE_PATH, `merged ${m.file?.fullPath ?? m.file}, not ${SHARE_PATH}`);
    o.view.close();
  },
  'merge-local-file': async (mod) => {
    const o = await install(mod);
    const dialog = await mergeDialog(o);
    await fromComputer(dialog);
    const m = await merged(o, dialog);
    const name = m.file?.fileName ?? m.file?.name;
    expect(name === 'hits-2.csv', `merged ${name ?? m.file}, not hits-2.csv`);
    o.view.close();
  },
  'merge-molecules': async (mod) => {
    const o = await install(mod);
    const dialog = await mergeDialog(o);
    await fromComputer(dialog);
    const molecules = () => inputs(dialog).find((i) => /^molecules?$/i.test(String(i.state.label).trim()))?.control;
    expect(await until(() => molecules() !== undefined && offered(molecules()).length > 0, 2000),
      'no Molecules choice over the file\'s columns');
    const items = offered(molecules()).map((i) => i.key);
    for (const name of ['Molecule', 'Scaffold'])
      expect(items.includes(name), `the file's molecule column ${name} is not among ${items.join(', ')}`);
    expect(!items.includes('ID'), 'a column that is not a molecule is offered');
    const hidden = items.filter((n) => n.startsWith('~'));
    expect(hidden.length === 0, `offers ${hidden.join(', ')}, a column the platform keeps out of every list`);
    const m = await merged(o, dialog);
    expect(m.column === 'Molecule', `merged by ${JSON.stringify(m.column)}, not Molecule`);
    o.view.close();
  },
  'add-row': async (mod) => {
    const {table, view} = await install(mod);
    await clickOk(onRibbon(view, /add new row/i, 'Add new row'), 'Add new row');
    expect(table.rowCount === 2, `the table has ${table.rowCount} rows after Add new row`);
    expect(errors().length === 0, `Add new row logged: ${errors()[0]}`);
    view.close();
  },
  'refresh-failure-handled': async (mod) => {
    const o = await install(mod, {refresh: () => Promise.reject(new Error('stand offline'))});
    await clickOk(onRibbon(o.view, /refresh/i, 'Refresh Campaign'), 'Refresh Campaign');
    await flush();
    const dropped = errors().filter((e) => e.startsWith('unhandled rejection'));
    expect(dropped.length === 0,
      'a failed refresh became an unhandled promise rejection: the action dropped its promise');
    o.view.close();
  },
  'gone-with-view': async (mod) => {
    const base = liveScopes();
    const subs = streams();
    const o = await install(mod);
    await cancel(await submitDialog(o));
    const dialog = await mergeDialog(o);
    await fromComputer(dialog);
    await until(() => o.calls.read.length === 1);
    await wait(60);
    await cancel(dialog);
    o.view.close();
    await flush();
    expect(liveScopes() === base, `${liveScopes() - base} scopes outlive the closed view`);
    expect(streams() === subs, `${streams() - subs} platform subscriptions outlive the closed view`);
  },
};
