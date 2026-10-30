// What graders drive a solution with: user gestures, u2's public introspection, logged errors, live scopes.
import {readdirSync, readFileSync} from 'node:fs';
import {basename, dirname, join} from 'node:path';
import {Control, Scope} from 'datagrok-api/u2core';
import {fire, flush, resetDom} from '../../tests/dom-shim.js';
import * as css from './css.mjs';

export {fire, flush, css};

/** The platform's molecule drawing: a `w`×`h` element that says what it draws. */
export function drawMolecule(smiles, w = 150, h = 150) {
  const el = document.createElement('div');
  el.className = 'eval-molecule';
  el.dataset.smiles = smiles;
  el.style.width = `${w}px`;
  el.style.height = `${h}px`;
  return el;
}

(globalThis.__evalChem ??= {}).drawMolecule ??= drawMolecule;
globalThis.grok_UI_PickFileFromFiles ??= async () => null;

export const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

export class Fail extends Error {}

const CHECK_MS = 20_000;

/** Answers whether `condition()` came to hold within `ms`. */
export async function until(condition, ms = 1000) {
  for (const end = Date.now() + ms; Date.now() < end; await wait(5)) {
    if (condition())
      return true;
  }
  return condition();
}

export function expect(condition, message) {
  if (!condition)
    throw new Fail(message);
}

const logged = [];
console.error = (...args) => logged.push(args.map((a) => a instanceof Error ? a.message : String(a)).join(' '));
process.on('unhandledRejection', (e) => logged.push(`unhandled rejection: ${e instanceof Error ? e.message : e}`));

export const errors = () => logged.slice();

export const liveScopes = () => Scope.liveCount;

/** The running check's scope: the grader is the owner of what it asks a solution to build, as a page is. */
let checkOwner;

/** The solution's `name`, called under the check's scope. */
export function need(mod, name) {
  expect(typeof mod?.[name] === 'function', `the solution exports no function \`${name}\``);
  return (...args) => checkOwner.build(() => mod[name](...args));
}

/** The root of what the solution built (an element, a control, a dialog), put on the page. */
export function mount(x) {
  const root = x instanceof HTMLElement ? x : x?.root;
  expect(root instanceof HTMLElement, 'the solution returned neither an element nor a control');
  if (!root.isConnected)
    document.body.append(root);
  return root;
}

export const controlOf = (el) => Control.forElement(el);

export function inputs(root = document.body) {
  return root.querySelectorAll('.u2-input-root').map((el) => Control.forElement(el))
    .filter((c, i, all) => c !== undefined && all.indexOf(c) === i)
    .map((control) => ({control, state: control.getState()}));
}

const norm = (s) => String(s ?? '').replace(/\s+/g, ' ').trim().toLowerCase();

export function input(label, root = document.body) {
  const found = inputs(root).find((i) => norm(i.state.label) === norm(label));
  const labels = () => inputs(root).map((i) => i.state.label).join(', ') || 'none';
  expect(found !== undefined, `no input labelled "${label}" (found: ${labels()})`);
  return found.control;
}

/** The row captioned `caption` — a u2 input, or a table row whose first cell says it — with the text it shows. */
export function row(caption, root = document.body) {
  const found = inputs(root).find((i) => norm(i.state.label) === norm(caption));
  if (found !== undefined)
    return {el: found.control.root, text: String(found.state.text ?? ''), control: found.control};
  const tr = root.querySelectorAll('tr').find((r) => norm(r.children[0]?.textContent) === norm(caption));
  expect(tr !== undefined, `no row captioned "${caption}"`);
  return {el: tr, text: tr.children.at(-1).textContent.trim()};
}

export function field(control) {
  const editor = control.root.querySelector('[data-u2-part="editor"]') ?? control.root;
  return editor.matches('input, textarea, select') ? editor : editor.querySelector('input, textarea, select') ?? editor;
}

export async function type(control, text) {
  const el = field(control);
  el.focus();
  el.value = text;
  fire(el, 'input');
  fire(el, 'change');
  await flush();
}

/** What a choice offers, as its item keys and the labels the user reads. */
export function offered(control) {
  const list = control.getItems?.({all: true})?.list;
  if (list !== undefined)
    return list.map((i) => ({key: String(i.key ?? i.label), label: String(i.label ?? i.key)}));
  return control.root.querySelectorAll('option').filter((o) => o.value !== '')
    .map((o) => ({key: o.value, label: o.textContent.trim()}));
}

/** Picks the option whose value or text is `value`, as a user picks what they read. */
export async function choose(control, value) {
  const el = field(control);
  expect(el.tagName === 'SELECT', `${control.getState().label} is not a list to pick from`);
  const option = el.querySelectorAll('option').find((o) => o.value === value) ??
    el.querySelectorAll('option').find((o) => o.textContent.trim() === value);
  expect(option !== undefined, `${control.getState().label} does not offer ${value}`);
  el.focus();
  el.value = option.value;
  fire(el, 'input');
  fire(el, 'change');
  await flush();
}

/** Presses `key` on `el`; answers whether nothing cancelled it. */
export async function press(el, key) {
  const done = fire(el, 'keydown', {key});
  fire(el, 'keyup', {key});
  await flush();
  return done;
}

/** Why the pointer cannot operate `el` (hidden, disabled, CSS out of reach), or null when it can. */
export function blocked(el) {
  if (!el.isConnected)
    return 'not on the page';
  for (let node = el; node && node.nodeType === 1; node = node.parentNode) {
    if (node.hasAttribute('hidden') || node.hasAttribute('inert') || node.inert === true)
      return 'hidden';
    if (css.computed(node, 'display', {hover: el}) === 'none')
      return 'display: none';
  }
  if (el.hasAttribute('disabled') || el.getAttribute('aria-disabled') === 'true' || el.closest('[disabled]'))
    return 'disabled';
  if (css.computed(el, 'visibility', {hover: el}) === 'hidden')
    return 'visibility: hidden';
  if (css.computed(el, 'pointer-events', {hover: el}) === 'none')
    return 'pointer-events: none';
  return null;
}

/** A user's click; refused, answering why, when `el` is not operable. */
export async function click(el) {
  expect(el instanceof HTMLElement, 'nothing to click');
  const why = blocked(el);
  if (why !== null)
    return why;
  for (const type of ['pointerdown', 'mousedown', 'pointerup', 'mouseup'])
    fire(el, type);
  el.click();
  await flush();
  return null;
}

export async function clickOk(el, what) {
  const why = await click(el);
  expect(why === null, `${what} cannot be clicked: ${why}`);
}

let filePicker = null;
const inputClick = HTMLInputElement.prototype.click;
HTMLInputElement.prototype.click = function() {
  if (this.type === 'file')
    filePicker = this;
  return inputClick.call(this);
};

/** Answers the OS file dialog the last click on an `<input type=file>` opened. */
export async function chooseFiles(files) {
  const picker = filePicker;
  expect(picker !== null, 'no file dialog was opened');
  filePicker = null;
  picker.files = files;
  const event = new Event('change', {bubbles: true});
  picker.dispatchEvent(event);
  if (typeof picker.onchange === 'function')
    picker.onchange(event);
  await wait(20);
}

export function localFile(name, text = 'smiles,id\nCCO,1\n') {
  return new File([text], name, {type: 'text/csv'});
}

/** What `el` is called: its label, title and text, and its Font Awesome glyph's name. */
export function nameOf(el) {
  const glyph = [el, ...el.querySelectorAll('i, span')].map((n) => /(?:^|\s)fa-([\w-]+)/.exec(n.className)?.[1])
    .find((g) => g !== undefined && !['fw', 'lg', 'sm', 'xs'].includes(g));
  return [el.getAttribute('aria-label'), el.getAttribute('title'), el.textContent, glyph].filter((x) => x).join(' ');
}

const CLICKABLE = 'button, a, [role="button"], [role="menuitem"], [tabindex], i, .u2-icon';

/** The outermost clickable under `root` whose name matches `re` and not `except`. */
export function named(root, re, except) {
  const all = root.querySelectorAll(CLICKABLE).filter((el) => re.test(nameOf(el)) && !except?.test(nameOf(el)));
  return all.find((el) => !all.some((other) => other !== el && other.contains(el))) ?? null;
}

export async function hoverText(el) {
  fire(el, 'pointerenter', {bubbles: false});
  await wait(350);
  const tip = document.querySelector('#u2-tooltip');
  const text = tip !== null && tip.style.display !== 'none' ? tip.textContent : '';
  fire(el, 'pointerleave', {bubbles: false});
  await flush();
  return text;
}

/** What `el` tells a user about itself: its description, title and hover text. */
export async function told(el) {
  const said = [el.getAttribute('aria-description'), el.getAttribute('title'), await hoverText(el)].filter((x) => x);
  return [...new Set(said)].join(' | ');
}

export function buttons(root, re) {
  return root.querySelectorAll('button, [role="button"]').filter((b) => re.test(norm(b.textContent)));
}

export function isOff(el) {
  return el.hasAttribute('disabled') || el.getAttribute('aria-disabled') === 'true';
}

async function cleanup() {
  for (const el of [...document.body.children]) {
    try {
      Control.forElement(el)?.dispose?.();
    }
    catch (_) {}
  }
  resetDom();
  await flush();
}

/** A reference solution (`good.ts`, `bad-*.ts`) alone; an agent's entry file with the `.ts` files beside it. */
export function sources(file) {
  const dir = dirname(file);
  const reference = /^(good|bad-.*)\.ts$/.test(basename(file));
  const files = reference ? [basename(file)] : readdirSync(dir).filter((f) => f.endsWith('.ts'));
  return files.map((f) => readFileSync(join(dir, f), 'utf8')).join('\n');
}

/** Runs each check on a clean page and prints the verdicts as the process's last line. */
export async function grade(checks, mod, loadError, ctx) {
  const results = [{id: 'loads', pass: loadError === undefined, ...loadError === undefined ? {} :
    {detail: String(loadError?.stack ?? loadError).split('\n').slice(0, 2).join(' | ')}}];
  for (const [id, check] of Object.entries(checks)) {
    if (loadError !== undefined) {
      results.push({id, pass: false, detail: 'the solution did not load'});
      continue;
    }
    logged.length = 0;
    checkOwner = Scope.root('eval check');
    let timer;
    try {
      await Promise.race([check(mod, ctx), new Promise((_, reject) =>
        timer = setTimeout(() => reject(new Fail(`the check did not finish in ${CHECK_MS / 1000} s`)), CHECK_MS))]);
      results.push({id, pass: true});
    }
    catch (e) {
      const thrown = String(e?.stack ?? e).split('\n').slice(0, 3).join(' | ');
      results.push({id, pass: false, detail: e instanceof Fail ? e.message : `threw ${thrown}`});
    }
    clearTimeout(timer);
    await cleanup();
    checkOwner.dispose();
  }
  process.stdout.write(`\n@@EVAL ${JSON.stringify(results)}\n`);
}
