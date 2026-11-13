import * as grok from 'datagrok-api/grok';
import {FileInfo} from '../../../tests/platform-doubles.mjs';
import {buttons, chooseFiles, clickOk, expect, flush, isOff, localFile, named, need, row, type, until}
  from '../../lib/harness.mjs';

const SHARE_PATH = 'System:DemoFiles/demog.csv';
const SHARE = new FileInfo(SHARE_PATH);

function open(mod) {
  const merged = [];
  const state = {merged, picks: 0};
  grok.dapi.files.exists = async (path) => path === SHARE_PATH;
  grok.dapi.files.list = async () => [SHARE];
  globalThis.grok_UI_PickFileFromFiles = async () => {
    state.picks++;
    return SHARE;
  };
  const dialog = need(mod, 'mergeDialog')((file) => void merged.push(file));
  const root = dialog?.root ?? dialog;
  expect(root instanceof HTMLElement && root.isConnected, 'mergeDialog() did not show a dialog');
  const merge = buttons(root, /^(merge|ok)$/)[0];
  expect(merge !== undefined, 'no MERGE button');
  return {...state, state, root, merge, source: row('Source file', root)};
}

const FILES = /\bfiles\b/i;
const FILES_BROWSER = /from (datagrok )?files|datagrok files|files browser|folder-tree/i;

function control(o, re, what, except) {
  const el = named(o.source.el, re, except) ?? named(o.root, re, except);
  expect(el !== null, `no control to ${what}`);
  return el;
}

async function merge(o) {
  await clickOk(o.merge, 'MERGE');
  expect(o.merged.length === 1, `onMerge ran ${o.merged.length} times`);
  return o.merged[0];
}

export const checks = {
  'merge-needs-file': async (mod) => {
    const o = open(mod);
    expect(isOff(o.merge), 'MERGE is on before a file is chosen');
    o.merge.click();
    await flush();
    expect(o.merged.length === 0, 'MERGE merged with no file');
  },
  'files-browser': async (mod) => {
    const o = open(mod);
    await clickOk(control(o, FILES, 'pick from Datagrok Files'), 'the Files picker');
    expect(o.state.picks === 1, 'the platform\'s Files browser was not opened');
    await until(() => !isOff(o.merge));
    const file = await merge(o);
    expect(file?.fullPath === SHARE_PATH, `merged ${file?.fullPath ?? file}, not ${SHARE_PATH}`);
  },
  'typed-share-path': async (mod) => {
    const o = open(mod);
    await type(o.source.control ?? {root: o.source.el}, SHARE_PATH);
    await until(() => !isOff(o.merge));
    const file = await merge(o);
    expect(file?.fullPath === SHARE_PATH, `merged ${file?.fullPath ?? file}, not ${SHARE_PATH}`);
  },
  'local-file': async (mod) => {
    const o = open(mod);
    const local = control(o, /open|choose|computer|upload|local|browse|folder/i, 'open a local file', FILES_BROWSER);
    await clickOk(local, 'Open file');
    await chooseFiles([localFile('hits.csv')]);
    await until(() => !isOff(o.merge));
    const file = await merge(o);
    const name = file?.fileName ?? file?.name;
    expect(name === 'hits.csv', `merged ${name ?? file}, not hits.csv`);
  },
};
