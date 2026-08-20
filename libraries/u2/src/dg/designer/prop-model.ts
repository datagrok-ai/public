/* The panel property model of a spec node (D-9): what the component declares or a plain HTML tag
   honors, the parameters of the function a source names, its bindings and its events — one model
   read by the writable panel and the read-only one alike. */
import {Component} from '../../core/component.js';
import type {PropertyLike} from '../../core/property-like.js';
import type {FieldOverride} from '../object-form.js';
import {htmlProps, isHtmlTag} from '../../spec/spec.js';
import type {SpecEventEntry, SpecNode} from '../../spec/spec.js';
import {backends} from '../../sources/backends.js';
import type {SpecNodeRef} from './node-ref.js';

/** Heading for the properties a component declares without a category of their own. */
const UNGROUPED = 'Properties';
/** The types with an editor of their own; anything else (`object`) renders as read-only JSON. */
const EDITABLE = new Set(['string', 'int', 'double', 'bool', 'string_list']);

/** One section of a node's panel property model: what a form renders, the mutable snapshot it
 * edits, and what the document says now — where a refused edit goes back to, and what an undo
 * refreshes to. */
export interface PropSection {
  title: string;
  /** The declared metadata with get/set closures over {@link values}; a bound or structured prop
   * gets no set, so the form routes it to the read-only field. */
  props: PropertyLike[];
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
  const groups = new Map<string, PropertyLike[]>();
  for (const prop of declared) {
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
  for (const [title, props] of groups) {
    const read = (): Record<string, unknown> => {
      const live = liveValues(ref);
      const values: Record<string, unknown> = {};
      for (const prop of props) {
        const value = live[prop.name] ?? node.props?.[prop.name];
        values[prop.name] = EDITABLE.has(typeOf(prop)) ? value : json(value);
      }
      return values;
    };
    sections.push({title, props: props.map((prop) => editable(prop, node)), values: read(), read});
  }
  return sections;
}

/** The declared prop as the form edits it: get/set closures over the section snapshot. A prop the
 * node binds is the context's to change — the Bindings field is where it is edited — and a
 * structured value has no editor yet: neither gets a set, so the form renders them read-only.
 * `string_list` is spoken as the form's `list`, which routes it to the list editor. */
function editable(prop: PropertyLike, node: SpecNode): PropertyLike {
  const name = prop.name;
  const type = typeOf(prop);
  return {...prop, propertyType: null, type: type === 'string_list' ? 'list' : type,
    get: (t) => t[name],
    set: node.bind?.[name] === undefined && EDITABLE.has(type) ? (t, v) => t[name] = v : undefined};
}

function liveValues(ref: SpecNodeRef): Record<string, unknown> {
  const built = ref.built();
  const live: Record<string, unknown> = {};
  if (built instanceof Component) {
    for (const prop of built.getProperties())
      live[prop.name] = prop.get?.(null);
  }
  return live;
}

function typeOf(prop: PropertyLike): string {
  return prop.propertyType ?? prop.type ?? 'string';
}

function json(value: unknown): string {
  return value === undefined ? '' : JSON.stringify(value);
}

export function stringProps(values: Record<string, string>, description: string): PropertyLike[] {
  return Object.keys(values).map((name) => ({name, type: 'string', description,
    get: (t: any) => t[name], set: (t: any, v: any) => t[name] = v}));
}

export function commitOnChange(values: Record<string, string>): Record<string, FieldOverride> {
  const overrides: Record<string, FieldOverride> = {};
  for (const name of Object.keys(values))
    overrides[name] = {commitOn: 'change'};
  return overrides;
}

/** What a source that names a function is called with: the prop its params live under, and the
 * function's own inputs — properties, so the panel edits them as it edits everything else. */
export function paramsOf(x: SpecNodeRef): {prop: string, inputs: PropertyLike[]} | null {
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

/** Every prop the component declares bindable, bound or not, plus anything the node binds beyond
 * them — a dotted sub-bind key, or a prop of a tag the registry no longer knows. An empty row is
 * where a binding is added, by hand or through the picker. */
export function bindsOf(x: SpecNodeRef): Record<string, string> {
  const bind = x.node.bind ?? {};
  const binds: Record<string, string> = {};
  for (const prop of x.instance.registry.get(x.node.tag)?.props ?? []) {
    if (prop.bindable)
      binds[prop.name] = bind[prop.name] ?? '';
  }
  for (const name of Object.keys(bind))
    binds[name] = bind[name];
  return binds;
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
