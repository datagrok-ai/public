import {clickOk, expect, field, input, named, need, press, row, type} from '../../lib/harness.mjs';

const CAMPAIGN = {name: 'KHL-3 hits', filePath: 'System:AppData/HitTriage/KHL-3.csv', status: 'In Progress'};

function open(mod) {
  const calls = {rename: [], download: 0, submit: 0};
  const dialog = need(mod, 'submitDialog')({...CAMPAIGN}, {rename: (n) => void calls.rename.push(n),
    download: () => void calls.download++, submit: () => void calls.submit++});
  const root = dialog?.root ?? dialog;
  expect(root instanceof HTMLElement && root.isConnected, 'submitDialog() did not show a dialog');
  return {calls, root};
}

function action(caption, re, what, dialog) {
  const el = named(row(caption, dialog).el, re) ?? named(dialog, re);
  expect(el !== null, `no ${what} control`);
  return el;
}

async function startEditing(mod) {
  const o = open(mod);
  await clickOk(action('Name', /edit|rename|pencil/i, 'Edit name', o.root), 'Edit name');
  return o;
}

export const checks = {
  'summary-rows': async (mod) => {
    const {root} = open(mod);
    const rows = [['Name', CAMPAIGN.name], ['File Path', CAMPAIGN.filePath], ['Status', CAMPAIGN.status]];
    for (const [label, value] of rows) {
      const r = row(label, root);
      expect(r.text === value, `${label} shows "${r.text}", not "${value}"`);
      expect(r.control === undefined || r.control.getState().readOnly === true, `${label} is not read-only`);
    }
  },
  'download-works': async (mod) => {
    const {calls, root} = open(mod);
    await clickOk(action('File Path', /download/i, 'Download', root), 'Download');
    expect(calls.download === 1, `actions.download() ran ${calls.download} times`);
  },
  'pencil-edits': async (mod) => {
    const {root} = await startEditing(mod);
    const name = input('Name', root);
    expect(name.getState().readOnly === false, 'Name is still read-only after Edit name');
    expect(document.activeElement === field(name), 'the cursor is not in Name');
  },
  'enter-saves': async (mod) => {
    const {calls, root} = await startEditing(mod);
    const name = input('Name', root);
    await type(name, 'KHL-3 final');
    await press(field(name), 'Enter');
    expect(calls.submit === 0, 'Enter in Name submitted the dialog');
    expect(root.isConnected, 'Enter in Name closed the dialog');
    expect(calls.rename.length === 1 && calls.rename[0] === 'KHL-3 final',
      `rename got ${JSON.stringify(calls.rename)}`);
    expect(name.getState().readOnly === true, 'Name stays editable after Enter');
  },
  'escape-cancels': async (mod) => {
    const {calls, root} = await startEditing(mod);
    const name = input('Name', root);
    await type(name, 'zzz');
    await press(field(name), 'Escape');
    expect(root.isConnected, 'Escape in Name closed the dialog');
    expect(calls.rename.length === 0, 'Escape renamed the campaign');
    const s = name.getState();
    expect(s.text === CAMPAIGN.name, `Name shows "${s.text}" after Escape`);
    expect(s.readOnly === true, 'Name stays editable after Escape');
  },
  'submit-runs': async (mod) => {
    const {calls, root} = open(mod);
    const submit = root.querySelectorAll('button').find((b) => b.textContent.trim().toLowerCase() === 'submit');
    expect(submit !== undefined, 'no SUBMIT button');
    await clickOk(submit, 'SUBMIT');
    expect(calls.submit === 1, `actions.submit() ran ${calls.submit} times`);
  },
};
