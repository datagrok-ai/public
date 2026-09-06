import {FilterError} from './model.js';
import type {FilterCondition, FilterKind, DomainConditionNode} from './model.js';
import type {FilterProperty} from './schema.js';
import {domainValue, kindOf} from './kinds.js';
// cycle with mask.ts is function-scoped only: `Masks` is used inside operator bodies, never at load time
import {Masks} from './mask.js';
import type {Mask, MaskCell, MaskColumnLike} from './mask.js';

export interface FilterOperator {
  id: string;
  label: string;
  arity: 0 | 1 | 2 | 'n';
  /** Applicable kinds (semType sets may list `['string']`). */
  kinds: FilterKind[];
  /** SemType-specific set; wins over core for that semType. */
  semType?: string;
  /** On a semType set: the property gets that set plus `is null`/`is not null` only — a molecule
   * column offers `has substructure`, not string `contains` beside it. */
  exclusive?: boolean;
  editor: 'default' | 'range' | 'list' | 'none';
  domain?: (c: FilterCondition, prop: FilterProperty, ctx: {now: Date}) => DomainConditionNode;
  mask?: (col: MaskColumnLike, c: FilterCondition, ctx: {now: Date}) => Mask;
  bitset?: (col: unknown, c: FilterCondition, signal: AbortSignal) => Promise<Mask>;
}

/** What an operator provider answers (a `meta.role: filterOperators` package function, say) — plain
 * data, so the producer needs no u2 import. `kinds` defaults to `['string']`, `editor` to `'default'`;
 * `bitset` answers the column's LSB-first words (a `Uint32Array` or its `ArrayBuffer`, the layout
 * `Mask.bits` and `BitArray.buffer` share). */
export interface FilterOperatorSet {
  semType: string;
  exclusive?: boolean;
  operators: {
    id: string;
    label: string;
    arity: 0 | 1 | 2 | 'n';
    kinds?: string[];
    editor?: FilterOperator['editor'];
    bitset(col: unknown, c: FilterCondition, signal: AbortSignal):
      Promise<{bits: Uint32Array | ArrayBuffer, length: number}>;
  }[];
}

export class OperatorRegistry {
  private _ops: FilterOperator[] = [];

  /** Same `id` + `semType` replaces in place. Returns the unregister function. */
  register(ops: FilterOperator | FilterOperator[]): () => void {
    const list = Array.isArray(ops) ? ops : [ops];
    for (const op of list) {
      const i = this._ops.findIndex((o) => o.id === op.id && o.semType === op.semType);
      if (i >= 0)
        this._ops[i] = op;
      else
        this._ops.push(op);
    }
    return () => {
      for (const op of list) {
        const i = this._ops.indexOf(op);
        if (i >= 0)
          this._ops.splice(i, 1);
      }
    };
  }

  /** A provider's whole set as operators of its `semType` (`exclusive` on each), the `bitset` words
   * wrapped into a `Mask` once here. Throws `FilterError` listing `checkSet`'s faults on a bad descriptor;
   * returns the unregister function for the whole set. */
  registerSet(set: FilterOperatorSet): () => void {
    const problems = OperatorRegistry.checkSet(set);
    if (problems.length > 0)
      throw new FilterError(`Bad operator set: ${problems.join('; ')}`);
    return this.register(set.operators.map((o): FilterOperator => ({
      id: o.id, label: o.label, arity: o.arity, kinds: (o.kinds ?? ['string']) as FilterKind[],
      semType: set.semType, exclusive: set.exclusive === true, editor: o.editor ?? 'default',
      bitset: async (col, c, signal) => {
        const r = await o.bitset(col, c, signal);
        return Masks.from(r.bits, r.length);
      },
    })));
  }

  /** The faults in a provider's descriptor, empty when it is usable — for a host to log once per provider. */
  static checkSet(set: unknown): string[] {
    const s = set as Partial<FilterOperatorSet> | null;
    if (typeof s !== 'object' || s === null)
      return ['not an object'];
    const problems: string[] = [];
    if (typeof s.semType !== 'string' || s.semType === '')
      problems.push('semType must be a non-empty string');
    if (s.exclusive !== undefined && typeof s.exclusive !== 'boolean')
      problems.push('exclusive must be a boolean');
    if (!Array.isArray(s.operators) || s.operators.length === 0)
      return [...problems, 'operators must be a non-empty array'];
    for (const [i, o] of s.operators.entries()) {
      const at = `operators[${i}]`;
      if (typeof o !== 'object' || o === null) {
        problems.push(`${at}: not an object`);
        continue;
      }
      if (typeof o.id !== 'string' || o.id === '')
        problems.push(`${at}: id must be a non-empty string`);
      if (typeof o.label !== 'string')
        problems.push(`${at}: label must be a string`);
      if (!ARITIES.includes(o.arity))
        problems.push(`${at}: arity must be 0, 1, 2 or 'n'`);
      if (o.kinds !== undefined && !(Array.isArray(o.kinds) && o.kinds.every((k) => KINDS.includes(k as FilterKind))))
        problems.push(`${at}: kinds must list ${KINDS.join(', ')}`);
      if (o.editor !== undefined && !EDITORS.includes(o.editor))
        problems.push(`${at}: editor must be ${EDITORS.join(', ')}`);
      if (typeof o.bitset !== 'function')
        problems.push(`${at}: bitset must be a function`);
    }
    return problems;
  }

  /** The property's semType set first (if any), then the core operators of its kind — only the
   * null tests of them when the set is exclusive. */
  for(prop: FilterProperty): FilterOperator[] {
    const kind = kindOf(prop);
    const applicable = (o: FilterOperator) => o.kinds.includes(kind);
    const sem = prop.semType ? this._ops.filter((o) => o.semType === prop.semType && applicable(o)) : [];
    const exclusive = sem.some((o) => o.exclusive === true);
    return [...sem, ...this._ops.filter((o) => o.semType === undefined && applicable(o) &&
      (!exclusive || NULL_TESTS.includes(o.id)))];
  }

  get(id: string, prop?: FilterProperty): FilterOperator | undefined {
    if (prop)
      return this.for(prop).find((o) => o.id === id);
    return this._ops.find((o) => o.id === id && o.semType === undefined) ?? this._ops.find((o) => o.id === id);
  }

  all(): FilterOperator[] {
    return this._ops.slice();
  }
}

const NULL_TESTS = ['is null', 'is not null'];
const ALL: FilterKind[] = ['string', 'int', 'float', 'bigint', 'datetime', 'bool', 'ref'];
const KINDS: FilterKind[] = [...ALL, 'string_list'];
const ARITIES: unknown[] = [0, 1, 2, 'n'];
const EDITORS: FilterOperator['editor'][] = ['default', 'range', 'list', 'none'];
const ORDERED: FilterKind[] = ['int', 'float', 'bigint', 'datetime'];
const LISTABLE: FilterKind[] = ['string', 'int', 'float', 'bigint', 'ref'];
const TEXT: FilterKind[] = ['string', 'string_list'];

/** Escapes LIKE metacharacters. */
export function escapeLike(s: string): string {
  return s.replace(/\\/g, '\\\\').replace(/%/g, '\\%').replace(/_/g, '\\_');
}

function likeValue(c: FilterCondition, shape: (v: string) => string): string {
  return c.options?.raw === true ? String(c.value) : shape(escapeLike(String(c.value ?? '')));
}

const node = (property: string, operator: string, value?: unknown): DomainConditionNode =>
  value === undefined ? {property, operator} : {property, operator, value};

const where = (test: (cell: MaskCell, value: any) => boolean) =>
  (col: MaskColumnLike, c: FilterCondition, ctx: {now: Date}): Mask => Masks.where(col, c, ctx, test);

/** Case-insensitive text test, the server's ILIKE; a `{raw: true}` LIKE pattern goes through a RegExp. */
const textMask = (test: (s: string, v: string) => boolean, negate: boolean = false) =>
  (col: MaskColumnLike, c: FilterCondition, ctx: {now: Date}): Mask => {
    const raw = c.options?.raw === true ? Masks.likeRegExp(String(c.value)) : null;
    const v = String(c.value ?? '').toLowerCase();
    return Masks.where(col, c, ctx, (cell) =>
      (raw ? raw.test(String(cell)) : test(String(cell).toLowerCase(), v)) !== negate);
  };

const regexMask = (negate: boolean) =>
  (col: MaskColumnLike, c: FilterCondition, ctx: {now: Date}): Mask => {
    const re = new RegExp(String(c.value), 'i');
    return Masks.where(col, c, ctx, (cell) => re.test(String(cell)) !== negate);
  };

const compare = (op: string, label: string): FilterOperator => ({
  id: op, label, arity: 1, kinds: ORDERED, editor: 'default',
  domain: (c, _prop, ctx) => node(c.property, op, domainValue(c.value, ctx)),
  mask: where(op === '>' ? (x, v) => x > v : op === '>=' ? (x, v) => x >= v : op === '<' ? (x, v) => x < v :
    (x, v) => x <= v),
});

const likeShape = (id: string, label: string, kinds: FilterKind[], operator: string,
  shape: (v: string) => string, test: (s: string, v: string) => boolean): FilterOperator => ({
  id, label, arity: 1, kinds, editor: 'default',
  domain: (c) => node(c.property, operator, likeValue(c, shape)),
  mask: textMask(test),
});

export const CORE_OPERATORS: FilterOperator[] = [
  {
    id: '=', label: 'equals', arity: 1, kinds: ALL, editor: 'default',
    domain: (c, _prop, ctx) => node(c.property, '=', domainValue(c.value, ctx)),
    mask: where((x, v) => x === v),
  },
  {
    id: '!=', label: 'not equals', arity: 1, kinds: ALL, editor: 'default',
    domain: (c, _prop, ctx) => node(c.property, '!=', domainValue(c.value, ctx)),
    mask: where((x, v) => x !== v),
  },
  compare('>', 'greater than'),
  compare('>=', 'at least'),
  compare('<', 'less than'),
  compare('<=', 'at most'),
  {
    id: 'between', label: 'between', arity: 2, kinds: ORDERED, editor: 'range',
    domain: (c, _prop, ctx) => {
      const [lo, hi] = domainValue(c.value, ctx) as unknown[];
      return [node(c.property, '>=', lo), 'and', node(c.property, '<=', hi)];
    },
    mask: where((x, [lo, hi]) => x >= lo && x <= hi),
  },
  {
    id: 'in', label: 'is one of', arity: 'n', kinds: LISTABLE, editor: 'list',
    domain: (c, _prop, ctx) => node(c.property, '=', domainValue(c.value, ctx)),
    mask: where((x, list: MaskCell[]) => list.includes(x)),
  },
  {
    id: 'not in', label: 'is not one of', arity: 'n', kinds: LISTABLE, editor: 'list',
    domain: (c, _prop, ctx) => node(c.property, '!=', domainValue(c.value, ctx)),
    mask: where((x, list: MaskCell[]) => !list.includes(x)),
  },
  likeShape('like', 'contains', TEXT, 'like', (v) => `%${v}%`, (s, v) => s.includes(v)),
  {
    id: '!like', label: 'does not contain', arity: 1, kinds: TEXT, editor: 'default',
    domain: (c) => node(c.property, 'not like', likeValue(c, (v) => `%${v}%`)),
    mask: textMask((s, v) => s.includes(v), true),
  },
  likeShape('starts', 'starts with', ['string'], 'like', (v) => `${v}%`, (s, v) => s.startsWith(v)),
  likeShape('ends', 'ends with', ['string'], 'like', (v) => `%${v}`, (s, v) => s.endsWith(v)),
  {
    id: 'matches', label: 'matches regex', arity: 1, kinds: ['string'], editor: 'default',
    domain: (c) => node(c.property, '~*', String(c.value)),
    mask: regexMask(false),
  },
  {
    id: '!matches', label: 'does not match', arity: 1, kinds: ['string'], editor: 'default',
    domain: (c) => node(c.property, '!~*', String(c.value)),
    mask: regexMask(true),
  },
  {
    id: 'fuzzy', label: 'is similar to', arity: 1, kinds: ['string'], editor: 'default',
    domain: (c) => {
      const value = String(c.value);
      const threshold = typeof c.options?.threshold === 'number' ? c.options.threshold : null;
      return [{property: c.property, operator: 'fuzzy', threshold, value}, 'or',
        node(c.property, 'like', `%${value}%`)];
    },
  },
  {
    id: 'is null', label: 'is empty', arity: 0, kinds: [...ALL, 'string_list'], editor: 'none',
    domain: (c) => node(c.property, '=', null),
    mask: (col) => Masks.nulls(col),
  },
  {
    id: 'is not null', label: 'is not empty', arity: 0, kinds: [...ALL, 'string_list'], editor: 'none',
    domain: (c) => node(c.property, '!=', null),
    mask: (col) => Masks.not(Masks.nulls(col)),
  },
];

/** The registry every surface consults — seeded with the core table at import. */
export const operators = new OperatorRegistry();
operators.register(CORE_OPERATORS);
