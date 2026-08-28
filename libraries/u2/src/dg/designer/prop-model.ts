/* The panel property model of a spec node (D-9): what the component declares or a plain HTML tag
   honors, the parameters of the function a source names, its bindings and its events — one model
   read by the writable panel and the read-only one alike. */
import {Component} from '../../core/component.js';
import type {IProperty} from '../../core/property-like.js';
import type {NamedProperty} from '../../core/widget-like.js';
import type {FieldOverride} from '../forms/object-form.js';
import {APPEARANCE_CATEGORY, APPEARANCE_PROPS} from '../../spec/appearance.js';
import {htmlProps, isHtmlTag} from '../../spec/spec.js';
import type {SpecEventEntry, SpecNode} from '../../spec/spec.js';
import {backends} from '../../sources/backends.js';
import type {SpecNodeRef} from './node-ref.js';

/** Heading for the properties a component declares without a category of their own. */
const UNGROUPED = 'Properties';
/** The categories the platform's own panel keeps at the bottom, in this order; the rest first-seen. */
const LAST = [APPEARANCE_CATEGORY, 'Misc', 'Description', 'Events'];
/** The types with an editor of their own — and the ones the document carries from a look-grid
 * edit; anything else (`object`) renders as read-only JSON. */
export const EDITABLE: ReadonlySet<string> = new Set(['string', 'int', 'double', 'bool', 'string_list']);
/** A frame reaches a node through its Bindings row only — never through a field, nor a grid row. */
export const BIND_ONLY = 'dataframe';
const TABLE = 'table';

export const TABLE_HINT = 'Bind `table` to a data source (Bindings › table)';

/** One section of a node's panel property model: what a form renders, the mutable snapshot it
 * edits, and what the document says now — where a refused edit goes back to, and what an undo
 * refreshes to. */
export interface PropSection {
  title: string;
  /** The declared metadata with get/set closures over {@link values}; a bound or structured prop
   * gets no set, so the form routes it to the read-only field. */
  props: NamedProperty[];
  values: Record<string, unknown>;
  read: () => Record<string, unknown>;
}

/** The panel property model of a node (D-9): what the component declares — grouped by category —
 * or what a plain HTML tag honors, plus the props its parent reads off it — a pane title — under
 * the parent's own section. Values come from the live component where one built, from the node
 * itself otherwise. */
export function propsFor(ref: SpecNodeRef): PropSection[] {
  const node = ref.node;
  const meta = ref.meta();
  const declared = meta ? meta.props : isHtmlTag(node.tag) ? htmlProps(node.tag) : [];
  const groups = new Map<string, NamedProperty[]>();
  for (const prop of declared) {
    if (typeOf(prop) === BIND_ONLY)
      continue;
    const title = prop.category ?? UNGROUPED;
    const group = groups.get(title);
    if (group)
      group.push(prop);
    else
      groups.set(title, [prop]);
  }
  const parent = ref.parent();
  const childProps = parent ? ref.instance.registry.get(parent.tag)?.childProps : undefined;
  if (childProps !== undefined && childProps.length > 0)
    groups.set(`Parent (${parent!.tag})`, childProps);

  const sections: PropSection[] = [];
  for (const [title, props] of [...groups].sort((a, b) => rank(a[0]) - rank(b[0]))) {
    const read = (): Record<string, unknown> => {
      const live = liveValues(ref);
      const values: Record<string, unknown> = {};
      for (const prop of props) {
        const name = prop.name;
        const value = live[name] ?? node.props?.[name];
        values[name] = EDITABLE.has(typeOf(prop)) ? value : json(value);
      }
      return values;
    };
    sections.push({title, props: props.map((prop) => editable(prop, node)), values: read(), read});
  }
  return sections;
}

const rank = (title: string): number => LAST.indexOf(title) + 1;

/** The shared appearance-group names this node answers — by identity, so a component-own prop
 * that declares the category or collides by name stays out. */
export function sharedAppearance(x: SpecNodeRef): Set<string> {
  const meta = x.meta();
  return new Set(APPEARANCE_PROPS.filter((p) => meta === undefined || meta.props.includes(p))
    .map((p) => p.name));
}

/** A node whose props are the built object's own (a platform viewer's look): the writable panel
 * shows them in the platform's property grid, the read-only one reads them through the tier. */
export function propertyTier(ref: SpecNodeRef): boolean {
  const built = ref.built();
  return Component.is(built) && built.propertyTier;
}

/** A node that wants a frame — the renderer said so, or its bind-only `table` is unbound (a viewer
 * builds over an empty placeholder without one, so this is true for a built node too). */
export function missingTable(x: SpecNodeRef): boolean {
  return x.error()?.includes('needs a table') === true || (x.node.bind?.[TABLE] === undefined &&
    x.meta()?.props.some((prop) => prop.name === TABLE && typeOf(prop) === BIND_ONLY) === true);
}

/** The declared prop as the form edits it: get/set closures over the section snapshot. A prop the
 * node binds is the context's to change — the Bindings field is where it is edited — and a
 * structured value has no editor yet: neither gets a set, so the form renders them read-only.
 * `string_list` is spoken as the form's `list`, which routes it to the list editor. */
function editable(prop: NamedProperty, node: SpecNode): NamedProperty {
  const name = prop.name;
  const type = typeOf(prop);
  return {...prop, propertyType: undefined, type: type === 'string_list' ? 'list' : type,
    get: (t) => t[name],
    set: node.bind?.[name] === undefined && EDITABLE.has(type) ? (t, v) => t[name] = v : undefined};
}

function liveValues(ref: SpecNodeRef): Record<string, unknown> {
  const built = ref.built();
  const live: Record<string, unknown> = {};
  if (Component.is(built)) {
    for (const prop of built.getProperties())
      live[prop.name] = prop.get?.(built.propertyTarget);
  }
  return live;
}

export function typeOf(prop: IProperty): string {
  return prop.propertyType ?? prop.type ?? 'string';
}

function json(value: unknown): string {
  return value === undefined ? '' : JSON.stringify(value);
}

export function stringProps(values: Record<string, string>, description: string): NamedProperty[] {
  return Object.keys(values).map((name) => ({name, type: 'string', description,
    get: (t: any) => t[name], set: (t: any, v: any) => t[name] = v}));
}

export function commitOnChange(values: Record<string, unknown>): Record<string, FieldOverride> {
  const overrides: Record<string, FieldOverride> = {};
  for (const name of Object.keys(values))
    overrides[name] = {commitOn: 'change'};
  return overrides;
}

/** What a source that names a function is called with: the prop its params live under, and the
 * function's own inputs — properties, so the panel edits them as it edits everything else. */
export function paramsOf(x: SpecNodeRef): {prop: string, inputs: NamedProperty[]} | null {
  const prop = x.instance.registry.get(x.node.tag)?.props.find((p) => p.subBindable);
  const func = x.node.props?.func;
  if (prop === undefined || typeof func !== 'string' || func === '')
    return null;
  const descriptor = backends.funcRunner?.find(func);
  return descriptor && descriptor.inputs.length > 0 ?
    {prop: prop.name, inputs: descriptor.inputs} : null;
}

export function paramValuesOf(node: SpecNode, prop: string): Record<string, unknown> {
  return {...node.props?.[prop] as Record<string, unknown> | undefined};
}

export function paramBinds(node: SpecNode, prop: string): Record<string, string> {
  const binds: Record<string, string> = {};
  for (const [key, path] of Object.entries(node.bind ?? {})) {
    if (key.startsWith(`${prop}.`))
      binds[key.slice(prop.length + 1)] = path;
  }
  return binds;
}

/** Every prop the component declares bindable — except the unbound appearance group, which would
 * bury the section under thirteen empty rows — plus anything the node binds beyond them: a dotted
 * sub-bind key, or a prop of a tag the registry no longer knows. An empty row is where a binding
 * is added, by hand or through the picker. */
export function bindsOf(x: SpecNodeRef): Record<string, string> {
  const bind = x.node.bind ?? {};
  const binds: Record<string, string> = {};
  for (const prop of x.instance.registry.get(x.node.tag)?.props ?? []) {
    if (prop.bindable && !APPEARANCE_PROPS.includes(prop))
      binds[prop.name] = bind[prop.name] ?? '';
  }
  for (const name of Object.keys(bind))
    binds[name] = bind[name];
  return binds;
}

/** The Bindings rows a panel shows: every bindable prop on a u2 control; on a property-tier node —
 * a viewer with forty look props — only what is bound, the rest one "Add binding…" row away. */
export function bindRowsOf(x: SpecNodeRef): Record<string, string> {
  const all = bindsOf(x);
  if (!propertyTier(x))
    return all;
  const bound: Record<string, string> = {};
  for (const name of Object.keys(x.node.bind ?? {}))
    bound[name] = all[name];
  return bound;
}

/** What "Add binding…" offers: the bindable props not bound yet — a property-tier node's look
 * props, a plain node's appearance group. */
export function unboundOf(x: SpecNodeRef): string[] {
  return (x.instance.registry.get(x.node.tag)?.props ?? [])
    .filter((p) => p.bindable && x.node.bind?.[p.name] === undefined).map((p) => p.name);
}

/** Every event the component declares, wired or not, plus anything the node wires beyond them. */
export function eventsOf(x: SpecNodeRef): Record<string, string> {
  const on = x.node.on ?? {};
  const events: Record<string, string> = {};
  for (const name of x.meta()?.events ?? [])
    events[name] = shownCommand(on[name]);
  for (const name of Object.keys(on))
    events[name] = shownCommand(on[name]);
  return events;
}

/** A structured entry shows its cmd; retyping the field replaces the whole entry with the
 * plain string — args editing arrives with the function picker. */
export function shownCommand(entry: SpecEventEntry | undefined): string {
  return entry === undefined ? '' : typeof entry === 'string' ? entry : entry.cmd;
}
