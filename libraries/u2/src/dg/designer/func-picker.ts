/* The function picker (Q6): one dialog over what the client registry knows, with the picked
   function's inputs — properties, like every other u2 form — as the arguments pane. It answers a
   name, its literal arguments and its bound ones; what that becomes (an event's `args`, a
   func-source's `params`) belongs to whoever opened it, so the commit stays with the panel and the
   tray, whose funnels are the one authority on what may be written. */
import * as DG from 'datagrok-api/dg';
import {Dialog} from '../../components/containers/dialog.js';
import {VirtualList} from '../../components/collections/list.js';
import {TextInput} from '../../components/inputs/text-input.js';
import {div, divV, span} from '../../core/elements.js';
import {Scope} from '../../core/scope.js';
import type {PropertyLike} from '../../core/property-like.js';
import {propertyForm} from '../forms/object-form.js';
import {handlerRenderer} from '../entities/entity.js';
import type {HandlerRenderer} from '../entities/entity.js';
import type {SpecEventEntry, SpecInstance} from '../../spec/spec.js';
import {bindPicker, bindPickerButton} from './bind-picker.js';

/** What the picker reads off a platform function — `DG.Func` satisfies it. */
export interface FuncLike {
  name: string;
  friendlyName?: string;
  description?: string;
  package?: {name?: string} | null;
  inputs?: PropertyLike[];
}

/** One row of the picker: the platform object it renders, and the plain strings the search matches
 * — read once, so typing never walks the dart boundary. Every one of them is drawn on the row, so
 * a hit is never on something the user cannot see. */
export interface FuncEntry {
  /** What a spec writes for this function. */
  name: string;
  label: string;
  description: string;
  search: string;
  func: FuncLike;
}

/** A picked function and what it is called with: literals by param name, and the params that
 * follow a `$.` path instead. */
export interface FuncPick {
  name: string;
  args: Record<string, unknown>;
  binds: Record<string, string>;
}

export interface FuncPickerOptions {
  title?: string;
  /** What the per-argument bind buttons walk. */
  instance: SpecInstance;
  /** What a row already wired to a function starts with. */
  pick?: FuncPick;
  onPick: (pick: FuncPick) => void;
}

/** What a spec can name: `cmd:<name>` takes an identifier. The client registry also answers
 * registrations keyed by a path — a file viewer's `/molecule.json` — which no spec can call and
 * which the picker therefore has no business offering. */
const CALLABLE = /^[A-Za-z_$][A-Za-z0-9_$]*$/;

/** The rows, by label, with the name each of them writes: the BARE name, because a client-registered
 * function has a namespace rather than a package and `find({package})` never matches it — qualified
 * only where the bare name would name more than one function. */
export function funcEntries(funcs: FuncLike[]): FuncEntry[] {
  const callable = funcs.filter((func) => CALLABLE.test(func.name));
  const counts = new Map<string, number>();
  for (const func of callable)
    counts.set(func.name, (counts.get(func.name) ?? 0) + 1);
  const entries = callable.map((func): FuncEntry => {
    const pack = func.package?.name;
    const name = pack && (counts.get(func.name) ?? 0) > 1 ? `${pack}:${func.name}` : func.name;
    const label = func.friendlyName ?? func.name;
    const description = func.description ?? '';
    return {name, label, description,
      search: `${name} ${label} ${description}`.toLowerCase(), func};
  });
  entries.sort((a, b) => a.label.localeCompare(b.label));
  return entries;
}

/** `Func.find`'s own `name` is an exact match, so the search is one enumeration and a local
 * substring filter over what a user can see: the name, the caption and the description. */
export function filterFuncs(entries: FuncEntry[], query: string): FuncEntry[] {
  const text = query.trim().toLowerCase();
  return text === '' ? entries : entries.filter((entry) => entry.search.includes(text));
}

/** The func's inputs as the form edits them: get/set closures over `values`, and — for a param the
 * binding owns — the path it follows, read-only. An input whose type has no editor renders
 * read-only too, which is what leaves a binding as its only way in. */
export function paramProps(inputs: PropertyLike[], values: Record<string, unknown>,
  bound: Record<string, string> = {}): PropertyLike[] {
  return inputs.map((prop): PropertyLike => {
    const name = prop.name;
    const caption = prop.caption ?? prop.friendlyName ?? name;
    const path = bound[name];
    if (path !== undefined)
      return {name, caption, type: 'string', description: `Bound to ${path}`, get: () => path};
    const type = prop.propertyType ?? prop.type ?? 'string';
    return {name, caption, type: type === 'string_list' ? 'list' : type,
      description: prop.description, choices: prop.choices, nullable: prop.nullable,
      get: (t) => t[name], set: (t, v) => t[name] = v};
  });
}

/** What a call is made with, minus the fields the user never filled — an empty editor is not a
 * value, as it is not one in the property panel. */
export function paramValues(args: Record<string, unknown>): Record<string, unknown> {
  const values: Record<string, unknown> = {};
  for (const [name, value] of Object.entries(args)) {
    if (value !== undefined && value !== null && value !== '')
      values[name] = value;
  }
  return values;
}

/** The event entry a pick becomes: the plain string form while there is nothing to pass, the
 * structured one otherwise. A literal that looks like a path is escaped once — `$$.` is what the
 * renderer unescapes when the event fires. */
export function eventEntry(pick: FuncPick): SpecEventEntry {
  const args = paramValues(pick.args);
  for (const [name, value] of Object.entries(args)) {
    if (typeof value === 'string' && value.startsWith('$.'))
      args[name] = `$${value}`;
  }
  Object.assign(args, pick.binds);
  return Object.keys(args).length === 0 ? `cmd:${pick.name}` : {cmd: `cmd:${pick.name}`, args};
}

/** The pick a wired event describes, or nothing where the picker cannot represent it: code-behind
 * and a component's own function are the other `cmd:` tiers, and neither is a platform function. */
export function eventPick(entry: SpecEventEntry | undefined): FuncPick | undefined {
  const cmd = entry === undefined ? '' : typeof entry === 'string' ? entry : entry.cmd;
  if (!cmd.startsWith('cmd:'))
    return undefined;
  const name = cmd.slice(4);
  if (name === '' || name.startsWith('#') || (name.includes('.') && !name.includes(':')))
    return undefined;
  const pick: FuncPick = {name, args: {}, binds: {}};
  const wired = typeof entry === 'object' ? entry.args ?? {} : {};
  for (const [key, value] of Object.entries(wired)) {
    if (typeof value === 'string' && value.startsWith('$.'))
      pick.binds[key] = value;
    else
      pick.args[key] = typeof value === 'string' && value.startsWith('$$.') ? value.slice(1) : value;
  }
  return pick;
}

/** Opens the picker over every function the client registry knows. `onPick` gets the function and
 * its arguments when OK closes it on a row; closing on nothing commits nothing. */
export function funcPicker(options: FuncPickerOptions): Dialog {
  const dialog = new Dialog(options.title ?? 'Pick a function');
  dialog.root.classList.add('u2-func-picker');
  const entries = funcEntries(DG.Func.find({}) as unknown as FuncLike[]);
  const renderer = handlerRenderer<FuncLike>();
  const search = dialog.run(() =>
    new TextInput({search: true, inline: true, placeholder: 'Search functions'}));
  const empty = span('', 'u2-picker-empty');
  const host = div([], 'u2-func-picker-params');
  const args = {...options.pick?.args};
  const binds = {...options.pick?.binds};
  let shown = entries;
  let picked: FuncEntry | undefined;
  // the arguments pane is rebuilt whenever what it describes changes; its form dies with it
  let pane = new Scope();
  dialog.own(() => pane.dispose());

  const showParams = (): void => {
    pane.dispose();
    pane = new Scope();
    host.textContent = '';
    if (picked === undefined) {
      host.append(span('Pick a function to set its arguments.', 'u2-func-picker-hint'));
      return;
    }
    const inputs = picked.func.inputs ?? [];
    // what another function's arguments left behind is not this one's
    const names = new Set(inputs.map((prop) => prop.name));
    for (const key of Object.keys(args)) {
      if (!names.has(key))
        delete args[key];
    }
    for (const key of Object.keys(binds)) {
      if (!names.has(key))
        delete binds[key];
    }
    if (inputs.length === 0) {
      host.append(span(`${picked.label} takes no arguments.`, 'u2-func-picker-hint'));
      return;
    }
    const form = Scope.runWith(pane, () =>
      propertyForm(paramProps(inputs, args, binds), args, {condensed: true}));
    for (const prop of inputs) {
      const input = form.input(prop.name);
      if (input === undefined)
        continue;
      input.root.setAttribute('data-u2-prop', prop.name);
      // owned by the dialog, not by the pane the pick rebuilds — the bind picker outlives it
      bindPickerButton(input, prop.name, () => dialog.run(() => bindPicker(options.instance, (path) => {
        binds[prop.name] = path;
        showParams();
      })));
    }
    host.append(form.root);
  };

  const list = dialog.run(() => new VirtualList<FuncEntry>({
    itemHeight: 28,
    keyOf: (entry) => entry.name,
    render: (entry, _index, row) => {
      row.dataset.u2Func = entry.name;
      return funcRow(entry, renderer);
    },
  }));
  dialog.effect(() => {
    const query = search.value.value;
    shown = filterFuncs(entries, query);
    list.setItems(shown);
    // a blank white box reads as a dialog that is still loading — two acceptance passes read it
    // exactly that way before waiting it out
    const found = shown.length > 0;
    empty.textContent = found ? '' : query.trim() === '' ?
      'No functions are registered.' : `No functions match "${query.trim()}".`;
    list.root.style.display = found ? '' : 'none';
    empty.style.display = found ? 'none' : '';
  });
  dialog.effect(() => {
    const entry = shown[list.selectedIndex.value];
    if (entry === picked)
      return;
    picked = entry;
    showParams();
  });
  const initial = options.pick === undefined ? -1 :
    shown.findIndex((entry) => entry.name === options.pick!.name);
  if (initial >= 0)
    list.selectedIndex.value = initial;

  dialog.add(divV([search.root, list.root, empty, host], 'u2-func-picker-body'));
  dialog.onCancel(() => dialog.dispose());
  dialog.onOK(() => {
    if (picked !== undefined)
      options.onPick({name: picked.name, args: {...args}, binds: {...binds}});
    dialog.dispose();
  });
  return dialog.show({modal: true, width: 560});
}

/** The row: the platform's icon, what the function is called, and what it says about itself — laid
 * out in one line, because the platform's own entity markup drew its icon over the caption and left
 * the descenders cut off at the row edge. Everything the search matches is on the row. */
function funcRow(entry: FuncEntry, renderer: HandlerRenderer<FuncLike>): HTMLElement {
  const parts: HTMLElement[] = [];
  // the base ObjectHandler answers the caption as text (js-api ui.ts:1745) — that is not an icon,
  // and putting it in a 16px box would just clip the name a second time
  const icon = renderer.icon(entry.func);
  if ((icon.textContent ?? '').trim() === '')
    parts.push(div([icon], 'u2-func-row-icon'));
  parts.push(span(entry.label, 'u2-func-row-label'));
  if (entry.name !== entry.label)
    parts.push(span(entry.name, 'u2-func-row-name'));
  if (entry.description !== '')
    parts.push(span(entry.description, 'u2-func-row-desc'));
  return div(parts, 'u2-func-row');
}
