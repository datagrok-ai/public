import {clone, find, isFlat, isGroup, optionsEquals, parentOf, valueEquals, walk} from './model.js';
import type {FilterGroup, FilterNode, FilterKind, FilterScalar, FilterValue, FilterProblem, Lock}
  from './model.js';
import type {IProperty} from '../property-like.js';
import type {ObjectRenderer} from '../object-renderer.js';
import type {Input, InputOptions} from '../input-base.js';

export interface FilterProperty extends IProperty {
  name: string;
  /** Target address of a ref: `'<schema>.<table>'` for a domain ref column, the `refType` of
   * an entity property. */
  ref?: string;
  kind?: FilterKind;
}
export interface FilterValueItem { value: FilterScalar; label?: string; count?: number }
export interface FilterSchema {
  properties: FilterProperty[];
  values?: (prop: FilterProperty, query: string, signal: AbortSignal) => Promise<FilterValueItem[]>;
  /** One FK hop. */
  resolveRef?: (prop: FilterProperty) => Promise<FilterSchema>;
  renderer?: (prop: FilterProperty) => ObjectRenderer<any> | undefined;
}
export interface FilterTemplate {
  root: FilterGroup;
  allowedProperties?: string[];
  allowedOperators?: Record<string, string[]>;
  allowAdvanced?: boolean;
  allowAdd?: boolean;
}
export interface FilterPlace { parentId: string; index: number }
export interface FilterDiff {
  changed: {id: string, property?: string, operator?: string, value?: FilterValue, options?: Record<string, unknown>}[];
  added: {parentId: string, index: number, node: FilterNode}[];
  removed: string[];
  /** Another parent, or another order among the siblings the template also has. */
  moved: {id: string, from: FilterPlace, to: FilterPlace}[];
}
/** null → the core default editor. */
export type FilterValueEditorFactory = (prop: FilterProperty, options: InputOptions<any>) => Input<any> | null;

const LOCK_RANK: Record<Lock, number> = {none: 0, value: 1, all: 2};

/** A literal schema for gallery pages and specs; `values` lists per property name are
 * offered filtered by the typed text. */
export function schema(props: IProperty[], values?: Record<string, FilterScalar[]>): FilterSchema {
  const result: FilterSchema = {
    properties: props.filter((p) => p.name).map((p) => ({...p, name: p.name!})),
  };
  if (values) {
    result.values = (prop, query) => {
      const q = query.toLowerCase();
      const items = (values[prop.name] ?? []).filter((v) => String(v).toLowerCase().includes(q));
      return Promise.resolve(items.map((value) => ({value, label: String(value)})));
    };
  }
  return result;
}

/** The property behind a path's head segment: the exact name, else the one property whose
 * name matches it case-insensitively (`AGE` for `age`; two case variants stay unknown). */
export function property(schema: FilterSchema, path: string): FilterProperty | null {
  const head = path.split('.')[0];
  const exact = schema.properties.find((p) => p.name === head);
  if (exact)
    return exact;
  const lower = head.toLowerCase();
  const loose = schema.properties.filter((p) => p.name.toLowerCase() === lower);
  return loose.length === 1 ? loose[0] : null;
}

interface Placed { node: FilterNode; parent: FilterGroup | null; index: number }

function index(root: FilterGroup): Map<string, Placed> {
  const map = new Map<string, Placed>();
  walk(root, (node, parent, index) => map.set(node.id, {node, parent, index}));
  return map;
}

/** What `value` changed against `template`, by node id — top-most additions and removals only;
 * `moved` ignores the index shifts additions and removals cause. */
export function diff(template: FilterGroup, value: FilterGroup): FilterDiff {
  const before = index(template);
  const after = index(value);
  const result: FilterDiff = {changed: [], added: [], removed: [], moved: []};
  // the order of the nodes both trees hold under the same parent, so an insertion or a removal
  // beside a node does not read as its move
  const common = (group: FilterGroup, other: Map<string, Placed>) =>
    group.nodes.filter((n) => other.get(n.id)?.parent?.id === group.id).map((n) => n.id);
  for (const [id, {node, parent, index}] of after) {
    const old = before.get(id);
    if (!old) {
      if (parent && before.has(parent.id))
        result.added.push({parentId: parent.id, index, node});
      continue;
    }
    if (parent && old.parent && (parent.id !== old.parent.id ||
      common(parent, before).indexOf(id) !== common(old.parent, after).indexOf(id))) {
      result.moved.push({id, from: {parentId: old.parent.id, index: old.index},
        to: {parentId: parent.id, index}});
    }
    if (isGroup(node) || isGroup(old.node))
      continue;
    const change: FilterDiff['changed'][number] = {id};
    if (node.property !== old.node.property)
      change.property = node.property;
    if (node.operator !== old.node.operator)
      change.operator = node.operator;
    if (!valueEquals(node.value, old.node.value))
      change.value = node.value;
    if (!optionsEquals(node.options, old.node.options))
      change.options = node.options ?? {};
    if (Object.keys(change).length > 1)
      result.changed.push(change);
  }
  for (const [id, {parent}] of before) {
    if (!after.has(id) && parent && after.has(parent.id))
      result.removed.push(id);
  }
  return result;
}

/** The strongest lock on the node or any of its ancestors. */
export function lockOf(root: FilterGroup, id: string): Lock {
  let lock: Lock = 'none';
  const raise = (n: FilterNode | null) => {
    if (n && LOCK_RANK[n.lock ?? 'none'] > LOCK_RANK[lock])
      lock = n.lock!;
  };
  raise(find(root, id));
  // a literal's nodes share the id `undefined`, so the walk would find the same parent forever
  const seen = new Set<FilterGroup>();
  for (let parent = parentOf(root, id); parent && !seen.has(parent); parent = parentOf(root, parent.id)) {
    seen.add(parent);
    raise(parent);
  }
  return lock;
}

/** `locked` problems for edits a template forbids — a group's lock covers its subtree. */
export function checkLocks(template: FilterTemplate, value: FilterGroup): FilterProblem[] {
  const problems: FilterProblem[] = [];
  const locked = (nodeId: string, message: string) => problems.push({nodeId, code: 'locked', message});
  const allowed = (c: FilterNode): string | null => {
    if (isGroup(c))
      return null;
    const head = c.property.split('.')[0];
    if (template.allowedProperties && !template.allowedProperties.includes(head))
      return `Property "${c.property}" is not offered here`;
    const ops = template.allowedOperators?.[head];
    return ops && !ops.includes(c.operator) ? `Operator "${c.operator}" is not offered for "${c.property}"` : null;
  };
  const {changed, added, removed, moved} = diff(template.root, value);
  for (const change of changed) {
    const lock = lockOf(template.root, change.id);
    if (lock === 'all' || (lock === 'value' && (change.property !== undefined || change.operator !== undefined ||
      change.options !== undefined)))
      locked(change.id, lock === 'all' ? 'This condition is locked' : 'Only the value of this condition may change');
    else {
      const message = allowed(find(value, change.id)!);
      if (message)
        locked(change.id, message);
    }
  }
  for (const {parentId, node} of added) {
    if (template.allowAdd === false || lockOf(template.root, parentId) !== 'none')
      locked(node.id, 'Conditions cannot be added here');
    else {
      const message = allowed(node);
      if (message)
        locked(node.id, message);
    }
  }
  for (const id of removed) {
    if (lockOf(template.root, id) !== 'none')
      locked(id, 'This condition cannot be removed');
  }
  // a reorder among siblings changes nothing a template fixes; another parent does
  for (const {id, from, to} of moved) {
    if ((template.allowAdd === false && from.parentId !== to.parentId) || lockOf(template.root, id) !== 'none' ||
      lockOf(template.root, to.parentId) !== 'none')
      locked(id, 'This condition cannot be moved');
  }
  walk(value, (node) => {
    const old = find(template.root, node.id);
    if (!old || !isGroup(node) || !isGroup(old) || lockOf(template.root, node.id) === 'none')
      return;
    if (node.op !== old.op || (node.not === true) !== (old.not === true))
      locked(node.id, 'This group is locked');
  });
  if (template.allowAdvanced === false && !isFlat(value))
    locked(value.id, 'Groups are not allowed here');
  return problems;
}

/** The template's tree as a starting value, ids kept so `diff` can follow it. */
export function applyTemplate(template: FilterTemplate): FilterGroup {
  return clone(template.root);
}
