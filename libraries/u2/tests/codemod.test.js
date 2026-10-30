/* scripts/codemod/u2-v1.mjs over fixtures: each is v0 code typed against a v0-shaped u2, migrated as the codemod
   writes it, and the result run again against a v1-shaped u2 comes back untouched (idempotent in the world it
   produces). The two u2s are stub declarations under a temporary node_modules, at the path the codemod knows u2 by. */

import {test, after} from 'node:test';
import assert from 'node:assert/strict';
import {mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {dirname, join} from 'node:path';
import {run, RULES} from '../scripts/codemod/u2-v1.mjs';

const SIGNALS = `
export declare class ReadonlySignal<T> { get value(): T; peek(): T; }
export declare class Signal<T> extends ReadonlySignal<T> { get value(): T; set value(x: T); }
export declare function signal<T>(x: T): Signal<T>;
export declare function computed<T>(fn: () => T): ReadonlySignal<T>;
export declare function span(text: string | ReadonlySignal<string>): HTMLElement;
export declare function actionsMenu(actions: Action[]): void;
`;

const V0 = {
  'src/index.d.ts': `${SIGNALS}
export declare class Scope {
  static runWith<T>(scope: Scope | undefined, fn: () => T): T;
  own(fn: () => void): void;
  dispose(): void;
}
export declare class Control {
  readonly scope: Scope;
  runInScope<T>(fn: () => T): T;
  own(fn: () => void): this;
  dispose(): void;
}
export interface InputOptions<T> {
  label?: string;
  enabled?: boolean | ReadonlySignal<boolean>;
  onEnter?: (e: KeyboardEvent, input: Input<T>) => void;
}
export declare class Input<T> extends Control {
  readonly value: Signal<T>;
  enabled: boolean;
}
export declare class TextInput extends Input<string> { constructor(options?: InputOptions<string>); }
export declare function subscribe<T>(source: {subscribe(fn: (x: T) => void): unknown}, handler: (x: T) => void,
  owner?: Scope | Control): () => void;
export interface Action { name: string; enabled?: boolean; run: () => void; }
export interface ButtonGroupItem { id: string; label: string; }
export declare class ButtonGroup extends Control { constructor(options: {items: ButtonGroupItem[]}); }
export interface ButtonOptions { disabled?: boolean | string | ReadonlySignal<boolean | string>; }
export declare function button(text: string, onClick: () => void, options?: ButtonOptions): HTMLButtonElement;
export interface TabOptions { id: string; label: string; }
export declare class TabStrip extends Control { addTab(tab: TabOptions): void; setTabLabel(id: string, label: string): void; }
export declare class Dialog extends Control {
  constructor(title: string);
  okEnabled(on: boolean | ReadonlySignal<boolean>): this;
}
`,
  'src/dg/index.d.ts': `import type {Control, Scope} from '../index.js';
export interface RibbonAction { icon: string; tooltip: string; label?: string; run: () => unknown; }
export interface TableViewChromeOptions { actions?: RibbonAction[]; owner?: Scope | Control; }
export declare function tableViewChrome(view: object, options?: TableViewChromeOptions): Control;
export declare function showInView(dialog: Control, view: object): Control;
export declare function viewSignals(view: object, owner?: Scope): {current: boolean};
export declare function ownView(view: object): Control;
export declare function collectionSelection(owner: Control, source: object, rows: object): void;
export declare class ViewHost {
  static of(view: object): ViewHost;
  ribbon(options?: TableViewChromeOptions): Control;
  show(dialog: Control): Control;
  readonly signals: {current: boolean};
}
`,
};

const V1 = {
  'src/index.d.ts': `${SIGNALS}
export type Reason = string | null | ReadonlySignal<string | null>;
export declare class Scope {
  static runWith<T>(scope: Scope | undefined, fn: () => T): T;
  build<T>(fn: () => T): T;
  own(fn: () => void): void;
  dispose(): void;
}
export declare class Control {
  readonly scope: Scope;
  build<T>(fn: () => T): T;
  subscribe<T>(source: {subscribe(fn: (x: T) => void): unknown}, handler: (x: T) => void): () => void;
  own(fn: () => void): this;
  dispose(): void;
}
export interface InputOptions<T> {
  label?: string;
  disabled?: Reason;
  onEnter?: (e: KeyboardEvent, input: Input<T>) => boolean | void | Promise<void>;
}
export declare class Input<T> extends Control {
  readonly signal: Signal<T>;
  get value(): T;
  set value(x: T);
  peek(): T;
  get disabled(): string | null;
  set disabled(reason: string | null);
}
export declare class TextInput extends Input<string> { constructor(options?: InputOptions<string>); }
export declare function subscribe<T>(source: {subscribe(fn: (x: T) => void): unknown}, handler: (x: T) => void):
  () => void;
export declare function bindText(el: HTMLElement, source: ReadonlySignal<string>, options?: {owner?: Control}): void;
export interface Action { text: string; name?: string; disabled?: Reason; run: () => void | Promise<void>; }
export interface ButtonGroupItem { id: string; text?: string; }
export declare class ButtonGroup extends Control { constructor(options: {items: ButtonGroupItem[]}); }
export interface ButtonOptions { disabled?: Reason; }
export declare function button(text: string, onClick: () => void, options?: ButtonOptions): HTMLButtonElement;
export interface TabOptions { id: string; title: string; }
export declare class TabStrip extends Control { addTab(tab: TabOptions): void; setTabTitle(id: string, title: string): void; }
export declare class Dialog extends Control {
  constructor(title: string);
  okDisabled(reason: Reason): this;
}
`,
  'src/dg/index.d.ts': `import type {Control} from '../index.js';
export interface RibbonAction { icon: string; tooltip: string; text?: string; run: () => void | Promise<void>; }
export interface RibbonOptions { actions?: RibbonAction[]; }
export declare function collectionSelection(source: object, rows: object, options: {owner: Control}): void;
export declare class ViewHost {
  static of(view: object): ViewHost;
  ribbon(options?: RibbonOptions): Control;
  show(dialog: Control): Control;
  readonly signals: {current: boolean};
}
`,
};

const SHARED = {
  'package.json': '{"name": "@datagrok-libraries/u2", "types": "src/index.d.ts"}\n',
  'dg/index.d.ts': 'export * from \'../src/dg/index.js\';\n',
  'src/components/inputs/text-input.d.ts': 'export {TextInput} from \'../../index.js\';\n',
  'src/core/elements.d.ts': 'export {subscribe} from \'../index.js\';\n',
  'src/internal.d.ts': 'export declare const secret: number;\n',
};

const root = mkdtempSync(join(tmpdir(), 'u2-codemod-'));
after(() => rmSync(root, {recursive: true, force: true}));

let worlds = 0;

/** A directory whose `@datagrok-libraries/u2` is `stubs`, holding `files`; answers the files' paths. */
function world(stubs, files) {
  const dir = join(root, `w${++worlds}`);
  for (const [path, text] of Object.entries({...SHARED, ...stubs}))
    write(join(dir, 'node_modules/@datagrok-libraries/u2', path), text);
  return Object.fromEntries(Object.entries(files).map(([name, text]) => [name, write(join(dir, name), text)]));
}

function write(path, text) {
  mkdirSync(dirname(path), {recursive: true});
  writeFileSync(path, text);
  return path;
}

/** Migrates `before` in the v0 world; checks it reads `expected` and that the v1 world leaves it alone. */
function migrates(name, files, expected, hands = []) {
  test(name, () => {
    const paths = world(V0, files);
    const report = run({paths: Object.values(paths), write: true, rules: new Set(Object.keys(RULES))});
    for (const [file, text] of Object.entries(expected))
      assert.equal(readFileSync(paths[file], 'utf8'), text, file);
    assert.deepEqual(report.hand.map((h) => `${h.rule}: ${h.message}`).sort(), [...hands].sort());
    const again = world(V1, Object.fromEntries(Object.keys(files).map((f) => [f, readFileSync(paths[f], 'utf8')])));
    const rerun = run({paths: Object.values(again), write: false, rules: new Set(Object.keys(RULES))});
    assert.equal(rerun.edits, 0, `the v1 world takes no edit: ${JSON.stringify(rerun.files)}`);
  });
}

const REASON = 'D1: write the reason the user reads in place of U2_REASON_NEEDED';

migrates('V1–V3: a value is plain, `peek()` reads it untracked, the signal is `.signal`; a cast is reported', {
  'a.ts': `import {TextInput, computed, span} from '@datagrok-libraries/u2';
export function f(found: unknown, maybe?: TextInput): HTMLElement {
  const name = new TextInput({label: 'Name'});
  const upper = computed(() => name.value.value.toUpperCase());
  const now = name.value.peek() + maybe?.value.peek();
  (found as {value: {value: string}}).value.value = 'x';
  const either = name.value;
  return span(name.value);
}
`}, {
  'a.ts': `import {TextInput, computed, span} from '@datagrok-libraries/u2';
export function f(found: unknown, maybe?: TextInput): HTMLElement {
  const name = new TextInput({label: 'Name'});
  const upper = computed(() => name.value.toUpperCase());
  const now = name.peek() + maybe?.peek();
  (found as {value: {value: string}}).value.value = 'x';
  const either = name.value;
  return span(name.signal);
}
`}, ['V1: a `.value.value` through a local type (`{value: {value: …}}`): if it is a u2 control, read `.value` — v1 ' +
  'holds the value itself, and the cast hides it from the type check',
'V3: a bare `.value` of a u2 value class where nothing asks for a signal: `.signal` passes the signal on; in v1 ' +
  '`.value` reads the value']);

test('a migrated file run again in the old world takes no edit: the v0 type check rejecting its reads marks it v1', () => {
  const paths = world(V0, {'a.ts': `import {TextInput, span} from '@datagrok-libraries/u2';
declare function check(actual: unknown, expected: unknown): void;
export function f(): HTMLElement {
  const name = new TextInput({label: 'Name'});
  name.value = 'x';
  check(name.value, 'x');
  return span(name.value);
}
`});
  const report = run({paths: Object.values(paths), write: false, rules: new Set(Object.keys(RULES))});
  assert.equal(report.edits, 0, JSON.stringify(report.files));
  assert.deepEqual(report.hand.map((h) => `${h.rule}:${h.line}`), ['V3:6']);
});

migrates('D1: `enabled` → a `disabled` reason, with the placeholder a person replaces; `okEnabled` → `okDisabled`', {
  'a.ts': `import {TextInput, Dialog, button, signal, computed} from '@datagrok-libraries/u2';
export function g(ready: boolean): void {
  const on = signal(false);
  const input = new TextInput({label: 'Name', enabled: ready});
  input.enabled = false;
  new Dialog('Save').okEnabled(on);
  button('Go', () => {}, {disabled: computed(() => ready || 'Not ready')});
}
`}, {
  'a.ts': `import {TextInput, Dialog, button, signal, computed} from '@datagrok-libraries/u2';
export function g(ready: boolean): void {
  const on = signal(false);
  const input = new TextInput({label: 'Name', disabled: ready ? null : U2_REASON_NEEDED});
  input.disabled = U2_REASON_NEEDED;
  new Dialog('Save').okDisabled(computed(() => on.value ? null : U2_REASON_NEEDED));
  button('Go', () => {}, {disabled: computed(() => ready || 'Not ready')});
}
`}, [REASON, REASON, REASON, 'D1: `disabled` mixes booleans and reasons: answer the reason while off and `null` while ' +
  'on — never `true` or `false`']);

migrates('O1–O3: `runInScope` → `build`, `Scope.runWith(s, …)` → `s.build(…)`, a positional owner → the receiver', {
  'a.ts': `import {Control, Scope, TextInput, subscribe} from '@datagrok-libraries/u2';
export function h(owner: Control, scope: Scope, stream: {subscribe(fn: (x: number) => void): unknown}): void {
  owner.runInScope(() => new TextInput());
  Scope.runWith(scope, () => new TextInput());
  subscribe(stream, (x) => console.log(x), owner);
}
`}, {
  'a.ts': `import {Control, Scope, TextInput} from '@datagrok-libraries/u2';
export function h(owner: Control, scope: Scope, stream: {subscribe(fn: (x: number) => void): unknown}): void {
  owner.build(() => new TextInput());
  scope.build(() => new TextInput());
  owner.subscribe(stream, (x) => console.log(x));
}
`});

const STREAM = 'declare const stream: {subscribe(fn: (x: number) => void): unknown};\n';
const owned = (name, call) => `export function ${name}(owner: Control): void {\n  ${call};\n}\n`;
const unmoved = (name) => owned(name, 'subscribe(stream, (x) => console.log(x), owner)');
const moved = (name) => owned(name, 'owner.subscribe(stream, (x) => console.log(x))');

migrates('O3: the import a moved owner leaves unused goes, that name only, and only while nothing else uses it', {
  'a.ts': `import {\n  Control, subscribe,\n  TextInput,\n} from '@datagrok-libraries/u2';\n${STREAM}` +
    `export const input = new TextInput();\n${unmoved('a')}`,
  'b.ts': `import {subscribe} from '@datagrok-libraries/u2/src/core/elements.js';\n` +
    `import type {Control} from '@datagrok-libraries/u2';\n${STREAM}${unmoved('b')}`,
  'c.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}export const api = {subscribe};\n${unmoved('c')}`,
  'd.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}export {subscribe};\n${unmoved('d')}`,
  'e.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}${unmoved('e')}` +
    `export const off = subscribe(stream, (x) => console.log(x));\n`,
  'f.ts': `import {Control, subscribe as on} from '@datagrok-libraries/u2';\n${STREAM}` +
    owned('f', 'on(stream, (x) => console.log(x), owner)'),
}, {
  'a.ts': `import {\n  Control,\n  TextInput,\n} from '@datagrok-libraries/u2';\n${STREAM}` +
    `export const input = new TextInput();\n${moved('a')}`,
  'b.ts': `import type {Control} from '@datagrok-libraries/u2';\n${STREAM}${moved('b')}`,
  'c.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}export const api = {subscribe};\n${moved('c')}`,
  'd.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}export {subscribe};\n${moved('d')}`,
  'e.ts': `import {Control, subscribe} from '@datagrok-libraries/u2';\n${STREAM}${moved('e')}` +
    `export const off = subscribe(stream, (x) => console.log(x));\n`,
  'f.ts': `import {Control} from '@datagrok-libraries/u2';\n${STREAM}${moved('f')}`,
});

const UNTYPED = {
  'a.ts': `declare const DG: any;
declare const el: Element;
export const rows = DG.U2.control(el).value.peek();
export const other = DG.grid.value.peek();
`,
  'b.ts': `import {TextInput} from '@datagrok-libraries/u2';
export function f(found: any): TextInput {
  found.value.value = 'x';
  return new TextInput();
}
`,
};

migrates('V1/V2 on an untyped (`any`) receiver are reported where it may be a u2 control: DG.U2, or a file using u2',
  UNTYPED, UNTYPED, ['V1: a `.value.value` on an untyped (`any`) receiver: if it is a u2 control, read `.value` — v1 ' +
    'holds the value itself, and the type check cannot see it', 'V2: a `.value.peek()` on an untyped (`any`) ' +
    'receiver: if it is a u2 control, read `.peek()` — v1 holds the value itself, and the type check cannot see it']);

migrates('A1 and R1: captions by role, and a `run` hands back nothing or a promise of nothing', {
  'a.ts': `import {TabStrip, actionsMenu} from '@datagrok-libraries/u2';
import type {Action} from '@datagrok-libraries/u2';
declare const rows: {addNew(): object};
declare function save(): Promise<number>;
export function k(tabs: TabStrip, open: Action): boolean {
  actionsMenu([{name: 'Add', run: () => rows.addNew()}, {name: 'Save', run: () => save()}]);
  tabs.addTab({id: 'a', label: 'A'});
  tabs.setTabLabel('a', 'B');
  return open.name === 'Open';
}
`}, {
  'a.ts': `import {TabStrip, actionsMenu} from '@datagrok-libraries/u2';
import type {Action} from '@datagrok-libraries/u2';
declare const rows: {addNew(): object};
declare function save(): Promise<number>;
export function k(tabs: TabStrip, open: Action): boolean {
  actionsMenu([{text: 'Add', run: () => void rows.addNew()}, {text: 'Save', run: async () => { await save(); }}]);
  tabs.addTab({id: 'a', title: 'A'});
  tabs.setTabTitle('a', 'B');
  return open.text === 'Open';
}
`});

migrates('W1 and I1: the view wrappers through `ViewHost`, and a path inside u2 through its entry point', {
  'a.ts': `import {Control} from '@datagrok-libraries/u2';
import {showInView, tableViewChrome, viewSignals} from '@datagrok-libraries/u2/src/dg/index.js';
import {TextInput} from '@datagrok-libraries/u2/src/components/inputs/text-input.js';
import {secret} from '@datagrok-libraries/u2/src/internal.js';
export function m(view: object, dialog: Control): boolean {
  tableViewChrome(view, {actions: [{icon: 'sync', tooltip: 'Refresh', run: () => {}}]});
  showInView(dialog, view);
  return viewSignals(view).current && new TextInput() !== null && secret > 0;
}
`}, {
  'a.ts': `import {Control} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {TextInput} from '@datagrok-libraries/u2';
import {secret} from '@datagrok-libraries/u2/src/internal.js';
export function m(view: object, dialog: Control): boolean {
  ViewHost.of(view).ribbon({actions: [{icon: 'sync', tooltip: 'Refresh', run: () => {}}]});
  ViewHost.of(view).show(dialog);
  return ViewHost.of(view).signals.current && new TextInput() !== null && secret > 0;
}
`}, ['I1: `@datagrok-libraries/u2/src/internal.js`: `secret` is on no entry point: ask u2 to export it, or do without']);

migrates('S1: a spec document\'s `enabled` → a `disabled` reason; a bound one renamed and reported', {
  'spec.json': `{"tag": "u2-form", "children": [
  {"tag": "u2-text-input", "props": {"label": "Name", "enabled": false}},
  {"tag": "u2-bool-input", "props": {"enabled": true}},
  {"tag": "u2-text-input", "bind": {"enabled": "$.on"}}
]}
`}, {
  'spec.json': `{"tag": "u2-form", "children": [
  {"tag": "u2-text-input", "props": {"label": "Name", "disabled": "U2_REASON_NEEDED"}},
  {"tag": "u2-bool-input", "props": {}},
  {"tag": "u2-text-input", "bind": {"disabled": "$.on"}}
]}
`}, ['S1: write the reason the user reads in place of "U2_REASON_NEEDED"',
  'S1: a bound `disabled` reads a reason (string | null): bind it to one']);

test('untyped JavaScript in the v1 shape is left alone: an `{owner}` option, or too few arguments for a v0 owner', () => {
  const paths = world(V1, {'a.js': `import {bindText} from '@datagrok-libraries/u2';
export function f(el, source, owner) {
  bindText(el, source, {owner});
  bindText(el, source);
}
`});
  assert.equal(run({paths: Object.values(paths), write: false, loose: true, rules: new Set(Object.keys(RULES))}).edits, 0);
});

test('a JavaScript file keeps a `run` that answers a value: R1 is a type contract', () => {
  const paths = world(V0, {'a.js': `import {actionsMenu} from '@datagrok-libraries/u2';
const log = [];
actionsMenu([{text: 'Log', run: () => log.push(1)}]);
`});
  assert.equal(run({paths: Object.values(paths), write: false, rules: new Set(['R1'])}).edits, 0);
});

migrates('D1: an `enabled` that may be undefined is off only when false; one that may be a signal is reported', {
  'a.ts': `import {TextInput, ReadonlySignal} from '@datagrok-libraries/u2';
export function d(on?: boolean, live?: boolean | ReadonlySignal<boolean>, options: {ready?: boolean} = {}): void {
  new TextInput({label: 'A', enabled: on});
  new TextInput({label: 'B', enabled: options.ready});
  new TextInput({label: 'C', enabled: live});
}
`}, {
  'a.ts': `import {TextInput, ReadonlySignal} from '@datagrok-libraries/u2';
export function d(on?: boolean, live?: boolean | ReadonlySignal<boolean>, options: {ready?: boolean} = {}): void {
  new TextInput({label: 'A', disabled: on === false ? U2_REASON_NEEDED : null});
  new TextInput({label: 'B', disabled: options.ready === false ? U2_REASON_NEEDED : null});
  new TextInput({label: 'C', disabled: live});
}
`}, [REASON, REASON, 'D1: a value that may be a signal (a live option): a signal is always truthy, so no test of it ' +
  'says whether it is on — write the reason by hand, a `computed(() => …)` over the signal\'s value']);

migrates('A1: a caption in a literal an unannotated `.map` callback returns, by the array the call is passed as', {
  'a.ts': `import {ButtonGroup} from '@datagrok-libraries/u2';
import type {Action} from '@datagrok-libraries/u2';
export function presets(labels: string[]): Action[] {
  new ButtonGroup({items: labels.map((label, at) => ({id: String(at), label}))});
  return labels.map((label) => ({name: label, run: () => {}}));
}
`}, {
  'a.ts': `import {ButtonGroup} from '@datagrok-libraries/u2';
import type {Action} from '@datagrok-libraries/u2';
export function presets(labels: string[]): Action[] {
  new ButtonGroup({items: labels.map((label, at) => ({id: String(at), text: label}))});
  return labels.map((label) => ({text: label, run: () => {}}));
}
`});

test('A1 in untyped JavaScript: a caption read over a list named for actions is reported, another `.name` is not', () => {
  const paths = world(V0, {'a.js': `export const names = (list) => list.actionsFor().map((a) => a.name);
export const rows = (list) => list.rows().map((r) => r.name);
`});
  const report = run({paths: Object.values(paths), write: false, loose: true, rules: new Set(Object.keys(RULES))});
  assert.equal(report.edits, 0);
  assert.deepEqual(report.hand.map((h) => `${h.rule}:${h.line}`), ['A1:1']);
});

migrates('O3: a positional owner on a function taken off a dynamic import moves too', {
  'a.js': `export async function select(owner, source, rows) {
  const {collectionSelection} = await import('@datagrok-libraries/u2/dg');
  collectionSelection(owner, source, rows);
}
`}, {
  'a.js': `export async function select(owner, source, rows) {
  const {collectionSelection} = await import('@datagrok-libraries/u2/dg');
  collectionSelection(source, rows, {owner: owner});
}
`});

test('I1: an explicit `exports` key of u2 is an entry (its test helpers); a `./src/…` key is a path inside u2', () => {
  const pkg = {name: '@datagrok-libraries/u2', types: 'src/index.d.ts', exports: {'.': {types: './src/index.d.ts'},
    './dg': {types: './dg/index.d.ts'}, './tests/dom-shim.js': './tests/dom-shim.js',
    './src/dg': {types: './src/dg/index.d.ts'}, './*': {types: './*.d.ts'}}};
  const paths = world({...V0, 'package.json': JSON.stringify(pkg),
    'tests/dom-shim.d.ts': 'export declare function fire(el: unknown, type: string): void;\n'}, {'a.ts': `import '@datagrok-libraries/u2/tests/dom-shim.js';
import {fire} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {ViewHost} from '@datagrok-libraries/u2/src/dg';
export const used = [fire, ViewHost];
`});
  const report = run({paths: Object.values(paths), write: true, rules: new Set(['I1'])});
  assert.deepEqual(report.hand, []);
  assert.equal(readFileSync(paths['a.ts'], 'utf8'), `import '@datagrok-libraries/u2/tests/dom-shim.js';
import {fire} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {ViewHost} from '@datagrok-libraries/u2/dg';
export const used = [fire, ViewHost];
`);
});
