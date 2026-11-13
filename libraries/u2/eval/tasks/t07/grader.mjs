import {buttons, clickOk, expect, flush, input, isOff, mount, need, told, type} from '../../lib/harness.mjs';

async function build(mod) {
  const saved = [];
  const root = mount(need(mod, 'newTemplateForm')((values) => void saved.push(values)));
  await flush();
  const save = buttons(root, /^save$/)[0];
  expect(save !== undefined, 'no SAVE button');
  return {root, save, saved};
}

async function offFor(save, field, when) {
  expect(isOff(save), `SAVE is on ${when}`);
  const why = await told(save);
  expect(why !== '', `SAVE is off ${when} and says nothing about why`);
  expect(new RegExp(field, 'i').test(why), `SAVE is off ${when} but says "${why}", which does not name ${field}`);
}

export const checks = {
  'save-says-why': async (mod) => {
    const {root, save} = await build(mod);
    await offFor(save, 'name', 'on an empty form');
    await type(input('Name', root), 'Kinase panel');
    await offFor(save, 'key', 'with no key');
    await type(input('Key', root), 'kh1');
    await offFor(save, 'key', 'with the key "kh1"');
  },
  'off-save-does-nothing': async (mod) => {
    const {save, saved} = await build(mod);
    save.click();
    await flush();
    expect(saved.length === 0, 'a press on the off SAVE saved');
  },
  'saves': async (mod) => {
    const {root, save, saved} = await build(mod);
    await type(input('Name', root), 'Kinase panel');
    await type(input('Key', root), 'KIN');
    await clickOk(save, 'SAVE');
    expect(saved.length === 1 && saved[0].name === 'Kinase panel' && saved[0].key === 'KIN',
      `onSave got ${JSON.stringify(saved)}`);
  },
};
