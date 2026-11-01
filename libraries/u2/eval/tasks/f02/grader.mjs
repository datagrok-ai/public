import * as grok from 'datagrok-api/grok';
import {Func, FuncParam} from '../../../tests/platform-doubles.mjs';
import {blocked, buttons, choose, clickOk, css, errors, expect, field, flush, input, inputs, isOff, liveScopes, mount,
  need, offered, sources, told, type, until, wait} from '../../lib/harness.mjs';
import {busySession, frame} from '../../lib/platform.mjs';

const functions = () => [
  new Func('descriptors', {friendlyName: 'Descriptors', inputs: [new FuncParam('table', 'dataframe'),
    new FuncParam('threshold', 'double', {defaultValue: 0.5})]}),
  new Func('similarity', {friendlyName: 'Similarity', inputs: [new FuncParam('table', 'dataframe'),
    new FuncParam('reference', 'dataframe'), new FuncParam('cutoff', 'double', {defaultValue: 0.7})]}),
  new Func('toxicity', {friendlyName: 'Toxicity', inputs: [new FuncParam('table', 'dataframe')]}),
];
const SAVED = {descriptors: {threshold: 0.6}};
const NAMED_TABLES = ['KHL-3', 'demog', 'refs'];

function deferred() {
  let resolve;
  const promise = new Promise((r) => resolve = r);
  return {promise, resolve};
}

const streams = () => grok.events.onTableAdded.count + grok.events.onTableRemoved.count;

async function build(mod) {
  busySession();
  const campaign = grok.shell.addTable(frame([{name: 'Molecule', type: 'string', semType: 'Molecule'}],
    [{Molecule: 'CCO'}], 'KHL-3'));
  const pending = deferred();
  const saves = [];
  const base = {scopes: liveScopes(), streams: streams()};
  const pane = need(mod, 'computeSettings')(campaign, {load: () => pending.promise,
    save: async (settings) => void saves.push(settings)});
  expect(typeof pane?.dispose === 'function', 'computeSettings() did not return a control');
  const root = mount(pane);
  await flush();
  const resolve = () => pending.resolve({functions: functions(), saved: structuredClone(SAVED)});
  return {pane, root, base, resolve, saves};
}

const norm = (s) => String(s ?? '').trim().toLowerCase();
const checkOf = (caption, root) => inputs(root)
  .find((i) => norm(i.state.label) === norm(caption) && typeof i.state.value === 'boolean')?.control;
const saysLoading = (root) => /loading/i.test(root.textContent) || [root, ...root.querySelectorAll('*')]
  .some((el) => el.getAttribute('aria-busy') === 'true' || el.getAttribute('role') === 'status');

async function shown(b) {
  b.resolve();
  expect(await until(() => checkOf('Toxicity', b.root) !== undefined, 2000),
    'the functions never show as checks after the load resolves' + (errors().length ? ` (logged: ${errors()[0]})` : ''));
}

async function loaded(mod) {
  const b = await build(mod);
  await shown(b);
  return b;
}

function saveButton(root) {
  const save = buttons(root, /^save$/)[0];
  expect(save !== undefined, 'no SAVE button');
  return save;
}

const ancestors = (el) => {
  const out = [];
  for (let node = el.parentElement; node; node = node.parentElement)
    out.push(node);
  return out;
};
const words = (node) => node.nodeType === 3 ? node.textContent : [...node.childNodes].map(words).join(' ');
const widest = (el) => Math.min(...['width', 'max-width'].map((p) => css.px(css.computed(el, p)) ?? Infinity));
const WHY = /\b(check|tick|turn on|switch on|enable|select)/i;
const INTERNAL = /\.(u2|d4|ui|grok)-|\[data-u2/;
const INLINE_STYLE = /\.style\.\w+\s*=(?!=)|\.style\.(setProperty|cssText)\b|setAttribute\(\s*['"]style['"]/;

const canon = (x) => x !== null && typeof x === 'object' ?
  Object.fromEntries(Object.keys(x).sort().map((k) => [k, canon(x[k])])) : x;
const same = (a, b) => JSON.stringify(canon(a)) === JSON.stringify(canon(b));

export const checks = {
  'loads': async (mod) => {
    const b = await build(mod);
    expect(await until(() => saysLoading(b.root), 500), 'nothing shows that it is loading');
    await shown(b);
    for (const [caption, on] of [['Descriptors', true], ['Similarity', false], ['Toxicity', false]]) {
      const check = checkOf(caption, b.root);
      expect(check !== undefined, `no check captioned ${caption}`);
      expect(check.getState().value === on, `${caption} is ${on ? 'not ' : ''}checked, and it is ` +
        `${on ? '' : 'not '}in the saved settings`);
    }
  },
  'params': async (mod) => {
    const {root} = await loaded(mod);
    const threshold = input('threshold', root).getState();
    expect(threshold.value === 0.6, `threshold shows ${JSON.stringify(threshold.value)}, not the saved 0.6`);
    const cutoff = input('cutoff', root).getState();
    expect(cutoff.value === 0.7, `cutoff shows ${JSON.stringify(cutoff.value)}, not its default 0.7`);
    input('reference', root);
    expect(!inputs(root).some((i) => norm(i.state.label) === 'table'),
      'the campaign\'s table is asked for, though it is always the campaign\'s');
  },
  'saves': async (mod) => {
    const {root, saves} = await loaded(mod);
    await clickOk(field(checkOf('Similarity', root)), 'the Similarity check');
    await choose(input('reference', root), 'refs');
    await type(input('cutoff', root), '0.8');
    await clickOk(saveButton(root), 'SAVE');
    expect(await until(() => saves.length === 1), `store.save ran ${saves.length} times`);
    const wanted = {descriptors: {threshold: 0.6}, similarity: {reference: 'refs', cutoff: 0.8}};
    expect(same(saves[0], wanted), `store.save got ${JSON.stringify(saves[0])}, not ${JSON.stringify(wanted)}`);
  },
  'locked-unchecked': async (mod) => {
    const {root} = await loaded(mod);
    const cutoff = input('cutoff', root);
    const s = cutoff.getState();
    expect(s.enabled === false || s.readOnly === true || blocked(field(cutoff)) !== null,
      'the cutoff of Similarity, which is not checked, can be changed');
  },
  'unchecked-says-why': async (mod) => {
    const {root} = await loaded(mod);
    const cutoff = input('cutoff', root);
    const check = checkOf('Similarity', root);
    const up = ancestors(cutoff.root);
    const region = up.find((el) => el.contains(check.root)) ?? root;
    const near = [field(cutoff), cutoff.root, ...up.slice(0, up.indexOf(region) + 1)];
    const said = [...await Promise.all(near.map(told)), words(region)].join(' | ');
    expect(WHY.test(said), 'Similarity\'s settings are locked while it is unchecked, and nothing there says why ' +
      `(it reads "${said.replace(/\s+/g, ' ').slice(0, 160)}")`);
  },
  'save-says-why': async (mod) => {
    const {root, saves} = await loaded(mod);
    await clickOk(field(checkOf('Descriptors', root)), 'the Descriptors check');
    const save = saveButton(root);
    expect(isOff(save), 'SAVE is on with no function checked');
    const why = await told(save);
    expect(why !== '', 'SAVE is off with no function checked and says nothing about why');
    expect(/function|check|select|choose|pick|nothing|none/i.test(why),
      `SAVE is off with no function checked but says "${why}"`);
    save.click();
    await flush();
    expect(saves.length === 0, 'a press on the off SAVE saved');
  },
  'nameless-table-open': async (mod) => {
    const {root} = await loaded(mod);
    expect(errors().length === 0, `with a table open that has no name, building the pane logged: ${errors()[0]}`);
    const tables = offered(input('reference', root)).map((i) => i.label).sort();
    expect(JSON.stringify(tables) === JSON.stringify(NAMED_TABLES),
      `with a table open that has no name, the reference offers ${JSON.stringify(tables)}`);
  },
  'dispose-while-loading': async (mod) => {
    const {pane, root, base, resolve} = await build(mod);
    pane.dispose();
    resolve();
    await wait(300);
    expect(liveScopes() === base.scopes, `${liveScopes() - base.scopes} scopes were built for a pane closed while it loaded`);
    expect(root.querySelectorAll('.u2-input-root').length === 0, 'rows were built into the pane after it was closed');
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
  'dispose-releases': async (mod) => {
    const {pane, base} = await loaded(mod);
    pane.dispose();
    await flush();
    expect(liveScopes() === base.scopes, `${liveScopes() - base.scopes} scopes outlive the closed pane`);
    expect(streams() === base.streams, `${streams() - base.streams} platform subscriptions outlive the closed pane`);
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
  'width-cap': async (mod) => {
    const {root} = await loaded(mod);
    const save = saveButton(root);
    const last = checkOf('Toxicity', root).root;
    const capped = [root, ...root.querySelectorAll('*')].filter((el) => el.contains(save) && el.contains(last))
      .some((el) => widest(el) <= 480);
    expect(capped, 'nothing keeps the pane within 480px');
  },
  'save-at-right': async (mod) => {
    const {root} = await loaded(mod);
    const save = saveButton(root);
    const row = save.parentElement;
    const flex = ['flex', 'inline-flex'].includes(css.computed(row, 'display'));
    const end = ['flex-end', 'end', 'right'].includes(css.computed(row, 'justify-content')) ||
      css.computed(save, 'margin-left') === 'auto';
    expect(flex && end, 'SAVE is not at the right');
  },
  'package-css-clean': async (mod) => {
    await loaded(mod);
    for (const rule of css.parse(css.packageCss())) {
      for (const part of rule.parts)
        expect(!INTERNAL.test(part.selector), `package CSS styles u2's insides: ${part.selector}`);
      const forced = rule.decls.find((d) => d.important);
      expect(forced === undefined, `package CSS forces ${forced?.prop} with !important`);
    }
  },
  'no-inline-styles': async (mod, ctx) => {
    const hit = INLINE_STYLE.exec(sources(ctx.file));
    expect(hit === null, `the solution writes inline styles (${hit?.[0]})`);
  },
};
