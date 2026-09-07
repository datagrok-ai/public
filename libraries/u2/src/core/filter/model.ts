import {TYPE} from 'datagrok-api/u2core';

/** The kinds the filter feature tells apart: the platform types they coincide with, plus u2's own
 * `float` (double, num, qnum) and `ref`. */
export const KIND = {
  STRING: TYPE.STRING, INT: TYPE.INT, FLOAT: 'float', BIG_INT: TYPE.BIG_INT, DATE_TIME: TYPE.DATE_TIME,
  BOOL: TYPE.BOOL, STRING_LIST: TYPE.STRING_LIST, REF: 'ref',
} as const;
export type FilterKind = typeof KIND[keyof typeof KIND];
export type Lock = 'none' | 'value' | 'all';
export interface FilterRef { type: string; id: string; name?: string }
/** A relative date — `'-1w'`, `'2d'`, `'now'` — resolved against `now` only when the tree leaves
 * for a target, so a "last 7 days" template stays live. */
export interface FilterSpan { span: string }
export type FilterScalar = string | number | boolean | null | Date | FilterRef | FilterSpan;
export type FilterValue = FilterScalar | FilterScalar[];

export interface FilterCondition {
  id: string;
  /** Schema path: `'name'`, `'author.login'`, `'category_id.name'`. */
  property: string;
  /** Operator registry id (`Filters.operators`). */
  operator: string;
  value?: FilterValue;
  /** Operator extras: `{threshold}` for fuzzy, `{raw: true}` for a verbatim LIKE pattern. */
  options?: Record<string, unknown>;
  lock?: Lock;
}
export interface FilterGroup { id: string; op: 'and' | 'or'; not?: boolean; nodes: FilterNode[]; lock?: Lock }
export type FilterNode = FilterCondition | FilterGroup;

export interface FilterProblem {
  /** null = the whole string (syntax). */
  nodeId: string | null;
  code: 'syntax' | 'unknown-property' | 'operator-not-applicable' | 'missing-value' | 'invalid-value' |
    'not-expressible' | 'evaluation' | 'locked';
  message: string;
  /** Syntax problems only: offsets into the text. */
  position?: {start: number, end: number};
}

/** The tree as JSON carries it: the node objects as they are, a `Date` value as `{date: ISO}`
 * (a span-marked one as its `{span}`), refs and spans verbatim. */
export interface FilterJsonGroup { id?: string; op: 'and' | 'or'; not?: boolean; nodes: FilterJsonNode[]; lock?: Lock }
export interface FilterJsonCondition {
  id?: string;
  property: string;
  operator: string;
  value?: unknown;
  options?: Record<string, unknown>;
  lock?: Lock;
}
export type FilterJsonNode = FilterJsonCondition | FilterJsonGroup;

export class FilterError extends Error {
  constructor(message: string, readonly problems: FilterProblem[] = []) {
    super(message);
    this.name = 'FilterError';
  }
}

/** The structural twin of `datagrok-api`'s `DomainConditionNode`: conditions and nested trees
 * joined by `'and'`/`'or'` strings. */
export interface DomainCondition { property: string; operator: string; value?: unknown; threshold?: number | null }
export type DomainConditionNode = DomainCondition | 'and' | 'or' | DomainConditionNode[];
export type DomainConditionTree = DomainConditionNode[];

let idPrefix = Math.random().toString(36).slice(2, 6);
let idCount = 0;

export function newId(): string {
  return `f${idPrefix}${++idCount}`;
}

/** Makes ids predictable — tests call `resetIds('')` for `f1`, `f2`, … */
export function resetIds(prefix: string = ''): void {
  idPrefix = prefix;
  idCount = 0;
}

export function group(op: 'and' | 'or' = 'and', nodes: FilterNode[] = [],
  extra?: {not?: boolean, lock?: Lock}): FilterGroup {
  const group: FilterGroup = {id: newId(), op, nodes};
  if (extra?.not !== undefined)
    group.not = extra.not;
  if (extra?.lock !== undefined)
    group.lock = extra.lock;
  return group;
}

export function cond(property: string, operator: string, value?: FilterValue,
  extra?: {options?: Record<string, unknown>, lock?: Lock}): FilterCondition {
  const cond: FilterCondition = {id: newId(), property, operator};
  if (value !== undefined)
    cond.value = value;
  if (extra?.options !== undefined)
    cond.options = extra.options;
  if (extra?.lock !== undefined)
    cond.lock = extra.lock;
  return cond;
}

export function isGroup(n: FilterNode): n is FilterGroup {
  return 'nodes' in n;
}

export function find(root: FilterGroup, id: string): FilterNode | null {
  if (root.id === id)
    return root;
  for (const n of root.nodes) {
    const found = n.id === id ? n : isGroup(n) ? find(n, id) : null;
    if (found)
      return found;
  }
  return null;
}

export function parentOf(root: FilterGroup, id: string): FilterGroup | null {
  for (const n of root.nodes) {
    if (n.id === id)
      return root;
    const parent = isGroup(n) ? parentOf(n, id) : null;
    if (parent)
      return parent;
  }
  return null;
}

export function walk(root: FilterGroup,
  visit: (n: FilterNode, parent: FilterGroup | null, index: number) => void): void {
  visit(root, null, -1);
  walkChildren(root, visit);
}

function walkChildren(group: FilterGroup,
  visit: (n: FilterNode, parent: FilterGroup | null, index: number) => void): void {
  group.nodes.forEach((n, i) => {
    visit(n, group, i);
    if (isGroup(n))
      walkChildren(n, visit);
  });
}

export function count(root: FilterGroup): number {
  let count = 0;
  walk(root, (n) => {
    if (!isGroup(n))
      count++;
  });
  return count;
}

/** No sub-groups and no negation — what simple mode can show. */
export function isFlat(root: FilterGroup): boolean {
  return root.not !== true && root.nodes.every((n) => !isGroup(n));
}

export function update(root: FilterGroup, id: string,
  patch: Partial<FilterCondition> | Partial<FilterGroup>): FilterGroup {
  if (root.id === id)
    return {...root, ...patch} as FilterGroup;
  return edit(root, id, (n) => ({...n, ...patch} as FilterNode));
}

export function insert(root: FilterGroup, parentId: string, index: number, node: FilterNode): FilterGroup {
  const into = (parent: FilterGroup): FilterGroup => {
    const nodes = parent.nodes.slice();
    nodes.splice(Math.min(Math.max(index, 0), nodes.length), 0, node);
    return {...parent, nodes};
  };
  if (root.id === parentId)
    return into(root);
  return edit(root, parentId, (n) => isGroup(n) ? into(n) : n);
}

export function remove(root: FilterGroup, id: string): FilterGroup {
  return edit(root, id, () => null);
}

/** Detaches the node, then inserts it at `index` in the parent as it is after the detach (the
 * designer's move-patch reading); a self/descendant target or a no-op lands on the same root. */
export function move(root: FilterGroup, id: string, parentId: string, index: number): FilterGroup {
  const node = find(root, id);
  const from = parentOf(root, id);
  if (!node || !from || id === parentId || (isGroup(node) && find(node, parentId)))
    return root;
  const detached = remove(root, id);
  const target = find(detached, parentId);
  if (!target || !isGroup(target))
    return root;
  const at = Math.min(Math.max(index, 0), target.nodes.length);
  if (from.id === parentId && from.nodes.indexOf(node) === at)
    return root;
  return insert(detached, parentId, at, node);
}

export function replace(root: FilterGroup, id: string, node: FilterNode): FilterGroup {
  if (root.id === id)
    return isGroup(node) ? node : root;
  return edit(root, id, () => node);
}

/** The root with every sub-group inlined, or null when nesting cannot go: a sub-group with
 * another connector and more than one node, or a negated one. */
export function flatten(root: FilterGroup): FilterGroup | null {
  if (isFlat(root))
    return root;
  if (root.not === true)
    return null;
  const nodes: FilterNode[] = [];
  const inline = (group: FilterGroup): boolean => {
    for (const n of group.nodes) {
      if (!isGroup(n))
        nodes.push(n);
      else if (n.not === true || (n.op !== root.op && n.nodes.length > 1) || !inline(n))
        return false;
    }
    return true;
  };
  return inline(root) ? {...root, nodes} : null;
}

export function clone(root: FilterGroup, freshIds: boolean = false): FilterGroup {
  const copy = (n: FilterNode): FilterNode => {
    const id = freshIds ? newId() : n.id;
    if (isGroup(n))
      return {...n, id, nodes: n.nodes.map(copy)};
    const cond: FilterCondition = {...n, id};
    if (Array.isArray(n.value))
      cond.value = n.value.slice();
    if (n.options !== undefined)
      cond.options = {...n.options};
    return cond;
  };
  return copy(root) as FilterGroup;
}

/** Structural equality: a ref's display name and the nodes' ids (unless asked for) are noise. */
export function equals(a: FilterNode, b: FilterNode, ignoreIds: boolean = true): boolean {
  if (!ignoreIds && a.id !== b.id)
    return false;
  if (isGroup(a) || isGroup(b)) {
    if (!isGroup(a) || !isGroup(b))
      return false;
    return a.op === b.op && (a.not === true) === (b.not === true) && (a.lock ?? 'none') === (b.lock ?? 'none') &&
      a.nodes.length === b.nodes.length && a.nodes.every((n, i) => equals(n, b.nodes[i], ignoreIds));
  }
  return a.property === b.property && a.operator === b.operator && (a.lock ?? 'none') === (b.lock ?? 'none') &&
    valueEquals(a.value, b.value) && optionsEquals(a.options, b.options);
}

export function optionsEquals(a: Record<string, unknown> | undefined, b: Record<string, unknown> | undefined): boolean {
  return plainEquals(a ?? {}, b ?? {});
}

export function valueEquals(a: FilterValue | undefined, b: FilterValue | undefined): boolean {
  if (Array.isArray(a) || Array.isArray(b)) {
    return Array.isArray(a) && Array.isArray(b) && a.length === b.length &&
      a.every((v, i) => valueEquals(v, b[i]));
  }
  if (a instanceof Date || b instanceof Date)
    return a instanceof Date && b instanceof Date && a.getTime() === b.getTime();
  if (isRef(a) || isRef(b))
    return isRef(a) && isRef(b) && a.type === b.type && a.id === b.id;
  if (isSpan(a) || isSpan(b))
    return isSpan(a) && isSpan(b) && a.span === b.span;
  return a === b;
}

/** A lossless plain-object form of the tree — every node key kept, operator ids untouched (a
 * semType operator the grammar cannot spell survives), `Date` values as `{date: ISO}`. */
export function toJson(root: FilterGroup): FilterJsonGroup {
  const scalar = (v: FilterScalar): unknown =>
    isSpan(v) ? {span: v.span} : v instanceof Date ? {date: v.toISOString()} : v;
  const node = (n: FilterNode): FilterJsonNode => {
    if (isGroup(n))
      return {...n, nodes: n.nodes.map(node)};
    const c: FilterJsonCondition = {...n};
    if (n.value !== undefined)
      c.value = Array.isArray(n.value) ? n.value.map(scalar) : scalar(n.value);
    return c;
  };
  return node(root) as FilterJsonGroup;
}

/** The inverse: a plain object with `nodes` (groups) or `property` + `operator` (conditions);
 * missing ids are made up, unknown keys dropped, `{date}` values become `Date`s. Throws
 * `FilterError` on anything else — a state written by hand, a Dart map that did not convert. */
export function fromJson(json: unknown): FilterGroup {
  const plain = (v: unknown): v is Record<string, unknown> => {
    const proto = typeof v === 'object' && v !== null ? Object.getPrototypeOf(v) : undefined;
    return proto === Object.prototype || proto === null;
  };
  const scalar = (v: unknown): FilterScalar => {
    if (plain(v) && typeof v.date === 'string')
      return new Date(v.date);
    if (v === null || ['string', 'number', 'boolean'].includes(typeof v) || isRef(v) || isSpan(v))
      return v as FilterScalar;
    throw new FilterError(`Filter value ${JSON.stringify(v)} is not a scalar`);
  };
  const lock = (v: unknown): Lock => {
    if (v !== 'none' && v !== 'value' && v !== 'all')
      throw new FilterError(`Filter lock "${v}" is not none, value or all`);
    return v;
  };
  const node = (j: unknown): FilterNode => {
    if (!plain(j))
      throw new FilterError('Filter node is not a plain object');
    const id = typeof j.id === 'string' && j.id !== '' ? j.id : newId();
    if (Array.isArray(j.nodes)) {
      if (j.op !== 'and' && j.op !== 'or')
        throw new FilterError(`Filter group connector "${j.op}" is not and/or`);
      const group: FilterGroup = {id, op: j.op, nodes: j.nodes.map(node)};
      if (j.not !== undefined)
        group.not = j.not === true;
      if (j.lock !== undefined)
        group.lock = lock(j.lock);
      return group;
    }
    if (typeof j.property !== 'string' || typeof j.operator !== 'string')
      throw new FilterError('Filter condition has no property or operator');
    const cond: FilterCondition = {id, property: j.property, operator: j.operator};
    if (j.value !== undefined)
      cond.value = Array.isArray(j.value) ? j.value.map(scalar) : scalar(j.value);
    if (plain(j.options))
      cond.options = {...j.options};
    if (j.lock !== undefined)
      cond.lock = lock(j.lock);
    return cond;
  };
  const root = node(json);
  if (!isGroup(root))
    throw new FilterError('Filter root is not a group');
  return root;
}

export function isRef(v: unknown): v is FilterRef {
  return typeof v === 'object' && v !== null && typeof (v as FilterRef).id === 'string' &&
    typeof (v as FilterRef).type === 'string';
}

export function isSpan(v: unknown): v is FilterSpan {
  return typeof v === 'object' && v !== null && typeof (v as FilterSpan).span === 'string';
}

function plainEquals(a: unknown, b: unknown): boolean {
  if (a === b)
    return true;
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null)
    return false;
  if (Array.isArray(a) || Array.isArray(b)) {
    return Array.isArray(a) && Array.isArray(b) && a.length === b.length &&
      a.every((v, i) => plainEquals(v, b[i]));
  }
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length &&
    keys.every((k) => plainEquals((a as Record<string, unknown>)[k], (b as Record<string, unknown>)[k]));
}

/** Path-copies from the root down to the node `id`; `change` answers its replacement (null removes
 * it, the same node keeps it). Untouched siblings and subtrees keep their identity; an unknown
 * id keeps the root's. */
function edit(group: FilterGroup, id: string, change: (n: FilterNode) => FilterNode | null): FilterGroup {
  const nodes = group.nodes;
  for (let i = 0; i < nodes.length; i++) {
    const n = nodes[i];
    const next = n.id === id ? change(n) : isGroup(n) ? edit(n, id, change) : n;
    if (next === n)
      continue;
    const copy = nodes.slice();
    if (next === null)
      copy.splice(i, 1);
    else
      copy[i] = next;
    return {...group, nodes: copy};
  }
  return group;
}
