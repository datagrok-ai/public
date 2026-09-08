import {FilterError, KIND, cond, group, isGroup, isRef, isSpan} from './model.js';
import type {DomainCondition, DomainConditionNode, DomainConditionTree, FilterCondition, FilterGroup, FilterNode,
  FilterProblem, FilterScalar, FilterValue} from './model.js';
import {property as propertyOf} from './schema.js';
import type {FilterProperty, FilterSchema} from './schema.js';
import {validate} from './validate.js';
import type {FilterTarget} from './validate.js';
import {fuzzyPair, isCondition, negate, parseTree} from './grammar.js';
import {spanOf} from '../span.js';
import {escapeLike, operators} from './operators.js';
import {kindOf} from './kinds.js';

const BARE_PATH = /^[A-Za-z]\w*(\.[A-Za-z]\w*)*$/;
const ESCAPED: Record<string, string> = {'\\': '\\\\', '"': '\\"', '\'': '\\\'', '\n': '\\n', '\r': '\\r', '\t': '\\t',
  '\b': '\\b', '\f': '\\f'};
const NAME_ESCAPES: Record<string, string> = {
  '\\': '\\', ']': ']', '/': '/', '"': '"', '\'': '\'', 'b': '\b', 'f': '\f', 'n': '\n', 'r': '\r', 't': '\t',
};

/** Whether the `%` or `_` at `i` is a wildcard (an even run of backslashes before it). */
function wildcardAt(s: string, i: number): boolean {
  let slashes = 0;
  for (let j = i - 1; j >= 0 && s[j] === '\\'; j--)
    slashes++;
  return slashes % 2 === 0;
}

function hasWildcard(s: string): boolean {
  for (let i = 0; i < s.length; i++) {
    if ((s[i] === '%' || s[i] === '_') && wildcardAt(s, i))
      return true;
  }
  return false;
}

const unescapeLike = (s: string): string => s.replace(/\\([\\%_])/g, '$1');

/** `[a\]b].c` → `a]b.c`: bracket segments lose their brackets and escapes. */
function unbracket(path: string): string {
  if (!path.includes('['))
    return path;
  return path.replace(/\[((?:[^\]\\]|\\.)*)\]/g, (_m, inner: string) =>
    inner.replace(/\\(u[0-9A-Fa-f]{4}|.)/g, (_e, c: string) =>
      c.length === 5 ? String.fromCharCode(parseInt(c.slice(1), 16)) : NAME_ESCAPES[c] ?? `\\${c}`));
}

function typed(v: unknown, prop: FilterProperty | null): FilterScalar {
  if (v instanceof Date) {
    const span = spanOf(v);
    return span === undefined ? v : {span};
  }
  if (typeof v !== 'string' || !prop)
    return v as FilterScalar;
  const kind = kindOf(prop);
  if (kind === KIND.REF && prop.ref && v !== '@current')
    return {type: prop.ref, id: v};
  if (kind === KIND.DATE_TIME && !Number.isNaN(Date.parse(v)))
    return new Date(v);
  return v;
}

function conditionFrom(c: DomainCondition, schema?: FilterSchema): FilterCondition {
  const proto = Object.getPrototypeOf(c);
  if (proto !== Object.prototype && proto !== null)
    throw new FilterError('condition is not a plain object');
  if (typeof c.property !== 'string')
    throw new FilterError('condition has no property');
  // the head segment takes the schema's spelling, so `AGE > 1` formats back as `age > 1`
  const segments = unbracket(c.property).split('.');
  const head = schema ? propertyOf(schema, segments[0]) : null;
  if (head)
    segments[0] = head.name;
  const property = segments.join('.');
  const prop = segments.length === 1 ? head : null;
  const value = (v: unknown): FilterValue => Array.isArray(v) ? v.map((x) => typed(x, prop)) : typed(v, prop);
  const op = c.operator;
  const v = c.value;
  if (op === 'fuzzy') {
    return cond(property, 'fuzzy', String(v),
      c.threshold == null ? undefined : {options: {threshold: c.threshold}});
  }
  if ((op === '=' || op === '!=') && Array.isArray(v))
    return cond(property, op === '=' ? 'in' : 'not in', value(v));
  if ((op === '=' || op === '!=') && v == null)
    return cond(property, op === '=' ? 'is null' : 'is not null');
  if (op === 'is' || op === 'is not')
    return cond(property, op === 'is' ? 'is null' : 'is not null');
  if (op === '~*' || op === '!~*')
    return cond(property, op === '~*' ? 'matches' : '!matches', String(v));
  if (op === 'like' || op === 'not like') {
    const s = String(v);
    const leading = s.startsWith('%');
    const trailing = s.length > (leading ? 1 : 0) && s.endsWith('%') && wildcardAt(s, s.length - 1);
    const inner = s.slice(leading ? 1 : 0, trailing ? -1 : undefined);
    const raw = hasWildcard(inner);
    if (leading && trailing && !raw)
      return cond(property, op === 'like' ? 'like' : '!like', unescapeLike(inner));
    if (op === 'like' && !raw && (leading || trailing))
      return cond(property, trailing ? 'starts' : 'ends', unescapeLike(inner));
    return cond(property, op === 'like' ? 'like' : '!like', s, {options: {raw: true}});
  }
  return cond(property, op, value(v));
}

/** A sub-list that is a fold — a fuzzy pair or a `between` pair — as one condition, else null. */
function foldedList(list: DomainConditionTree, schema?: FilterSchema): FilterCondition | null {
  const pair = fuzzyPair(list);
  if (pair)
    return conditionFrom(pair.fuzzy, schema);
  if (list.length !== 3 || list[1] !== 'and' || !isCondition(list[0]) || !isCondition(list[2]))
    return null;
  const [low, , high] = list as [DomainCondition, 'and', DomainCondition];
  if (low.operator !== '>=' || high.operator !== '<=' || low.property !== high.property)
    return null;
  return conditionFrom({property: low.property, operator: 'between', value: [low.value, high.value]}, schema);
}

function groupFrom(list: DomainConditionTree, schema?: FilterSchema): FilterGroup {
  const runs: FilterNode[][] = [[]];
  let placed = false;
  let sticky: 'and' | 'or' | null = null;
  let pending: 'and' | 'or' | null = null;
  for (const item of list) {
    if (item === 'and' || item === 'or') {
      sticky = item;
      if (placed && pending === null)
        pending = item;
      continue;
    }
    let node: FilterNode | null;
    if (isCondition(item))
      node = conditionFrom(item, schema);
    else {
      node = foldedList(item, schema);
      if (!node) {
        const sub = groupFrom(item, schema);
        node = sub.nodes.length === 0 ? null : sub.nodes.length === 1 ? sub.nodes[0] : sub;
      }
    }
    if (!node)
      continue;
    if (placed && (pending ?? sticky ?? 'and') === 'or')
      runs.push([node]);
    else
      runs[runs.length - 1].push(node);
    placed = true;
    pending = null;
  }
  if (runs.length === 1)
    return group('and', runs[0]);
  return group('or', runs.map((r) => r.length === 1 ? r[0] : group('and', r)));
}

/** AND-over-OR grouping with sticky connectors; `between`, `in`, `is null`, `like`
 * shapes and fuzzy pairs fold back; with a schema, ref and datetime values get typed. */
export function fromDomainTree(tree: DomainConditionTree | DomainCondition, schema?: FilterSchema): FilterGroup {
  let root = groupFrom(Array.isArray(tree) ? tree : [tree], schema);
  while (root.nodes.length === 1 && isGroup(root.nodes[0]))
    root = root.nodes[0];
  return root;
}

/** The inverse — spans resolve against `now`, `not` groups push down. Throws
 * `FilterError` for an operator without a domain form. */
export function toDomainTree(root: FilterGroup, options?: {now?: Date, schema?: FilterSchema}): DomainConditionTree {
  const now = options?.now ?? new Date();
  const emit = (group: FilterGroup): DomainConditionTree => {
    const out: DomainConditionTree = [];
    for (const n of group.nodes) {
      let piece: DomainConditionNode;
      if (isGroup(n)) {
        const sub = emit(n);
        if (sub.length === 0)
          continue;
        piece = sub.length === 1 ? sub[0] : sub;
      } else {
        const prop = options?.schema ? propertyOf(options.schema, n.property) : null;
        const op = operators.get(n.operator, prop ?? undefined) ?? operators.get(n.operator);
        if (!op?.domain)
          throw new FilterError(`Operator "${n.operator}" has no domain form`);
        piece = op.domain(n, prop ?? {name: n.property}, {now});
      }
      if (out.length > 0)
        out.push(group.op);
      out.push(piece);
    }
    return group.not === true ? negate(out) : out;
  };
  return emit(root);
}

function quote(s: string): string {
  return `"${s.replace(/[\\"'\n\r\t\b\f]|[\x00-\x1f]/g, (c) =>
    ESCAPED[c] ?? `\\u${c.charCodeAt(0).toString(16).padStart(4, '0')}`)}"`;
}

/** A property path the way `format` spells it: bare when every segment is an identifier,
 * `[…]`-bracketed with escapes otherwise. */
export function formatProperty(path: string): string {
  if (BARE_PATH.test(path))
    return path;
  return `[${path.replace(/[\\\]]|[\x00-\x1f]/g, (c) =>
    c === '\\' || c === ']' ? `\\${c}` : `\\u${c.charCodeAt(0).toString(16).padStart(4, '0')}`)}]`;
}

/** A value the way `format` spells it: strings quoted, `@current` and spans bare, lists
 * parenthesized, `''` for `undefined`. */
export function formatValue(v: FilterValue | undefined): string {
  if (v === undefined)
    return '';
  if (Array.isArray(v))
    return `(${v.map(formatValue).join(', ')})`;
  if (v === null || typeof v === 'boolean' || typeof v === 'number')
    return String(v);
  if (v instanceof Date)
    return spanOf(v) ?? quote(v.toISOString());
  if (isSpan(v))
    return v.span;
  if (isRef(v))
    return quote(v.id);
  return v === '@current' ? v : quote(v);
}

function formatCondition(c: FilterCondition): string {
  const prop = formatProperty(c.property);
  switch (c.operator) {
    case 'is null': return `${prop} = null`;
    case 'is not null': return `${prop} != null`;
    case 'between': {
      const [low, high] = Array.isArray(c.value) ? c.value : [];
      return `${prop} between ${formatValue(low)} and ${formatValue(high)}`;
    }
    case 'fuzzy': {
      const t = c.options?.threshold;
      return `${prop} fuzzy${typeof t === 'number' ? `(${t})` : ''} ${formatValue(c.value)}`;
    }
    case 'like': case '!like': case 'starts': case 'ends': {
      const s = String(c.value ?? '');
      if (c.options?.raw !== true)
        return `${prop} ${c.operator} ${quote(escapeLike(s))}`;
      const outer = (c.operator === 'like' || c.operator === '!like') && s.length > 1 && s.startsWith('%') &&
        s.endsWith('%') && wildcardAt(s, s.length - 1);
      return `${prop} ${c.operator} ${quote(outer ? s.slice(1, -1) : s)}`;
    }
    default: return `${prop} ${c.operator} ${formatValue(c.value)}`.trimEnd();
  }
}

/** The pieces of a group: conditions and parenthesized sub-groups, empties dropped. */
function formatPieces(group: FilterGroup): string[] {
  const pieces: string[] = [];
  for (const n of group.nodes) {
    const piece = isGroup(n) ? formatGroup(n, true) : formatCondition(n);
    if (piece !== '')
      pieces.push(piece);
  }
  return pieces;
}

function formatGroup(group: FilterGroup, nested: boolean): string {
  const pieces = formatPieces(group);
  if (pieces.length === 0)
    return '';
  const joined = pieces.length === 1 ? pieces[0] : pieces.join(` ${group.op} `);
  if (group.not === true)
    return `not (${pieces.length === 1 && joined.startsWith('(') ? joined.slice(1, -1) : joined})`;
  return nested && pieces.length > 1 ? `(${joined})` : joined;
}

/** The canonical string: `"`-quoted, bracketed names where needed, nested groups
 * parenthesized, single-node groups inlined, `''` for an empty tree. */
export function format(node: FilterNode): string {
  return isGroup(node) ? formatGroup(node, false) : formatCondition(node);
}

/** `parseTree` → `fromDomainTree` → `validate` (with a schema); syntax problems come back
 * with an empty root. */
export function parse(text: string, schema?: FilterSchema, target?: FilterTarget, options?: {now?: Date}):
  {root: FilterGroup, problems: FilterProblem[]} {
  const {tree, errors} = parseTree(text, options);
  if (errors.length > 0)
    return {root: group('and'), problems: errors};
  const root = fromDomainTree(tree, schema);
  return {root, problems: schema ? validate(root, schema, target) : []};
}
