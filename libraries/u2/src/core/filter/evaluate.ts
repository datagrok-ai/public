import {BitArray} from 'datagrok-api/u2core';
import {FilterError, KIND, isColumnRef, isGroup, isParam, isSpan} from './model.js';
import type {FilterCondition, FilterGroup, FilterKind, FilterScalar} from './model.js';
import {domainValue, kindOf} from './kinds.js';
import {resolveSpan} from '../span.js';
// cycle with operators.ts is function-scoped only: neither module reads the other at load time
import {operators} from './operators.js';

/** What a leaf needs from a column — `DG.Column`'s raw-data contract: `Int32Array` for ints
 * (`INT_NULL`), `Float32Array`/`Float64Array` for floats (`FLOAT_NULL`), `Float64Array` µs for
 * datetimes, `Int32Array` category indexes for strings (`''` is the null category), `Uint32Array`
 * bits for bools; bigints and string lists go through `get(i)`. */
export interface MaskColumnLike {
  name: string;
  type: string;
  semType?: string;
  length: number;
  getRawData(): ArrayLike<number>;
  categories?: string[];
  get?(i: number): unknown;
}

export interface MaskFrameLike { rowCount: number; column(name: string): MaskColumnLike | null }

/** A non-null cell converted to the column's kind: number (int, float, datetime µs), bigint, string, boolean. */
export type MaskCell = number | bigint | string | boolean;

export const INT_NULL = -2147483648;
export const FLOAT_NULL = 2.6789344063684636e-34;
const FLOAT_NULL_32 = Math.fround(FLOAT_NULL);
const NEVER = new AbortController().signal;

const isFloatNull = (v: number): boolean => v === FLOAT_NULL || v === FLOAT_NULL_32 || Number.isNaN(v);
const escapeRegExp = (s: string): string => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const COMPARE: Record<string, (x: MaskCell, y: MaskCell) => boolean> = {
  '=': (x, y) => x === y, '!=': (x, y) => x !== y, '<': (x, y) => x < y, '<=': (x, y) => x <= y,
  '>': (x, y) => x > y, '>=': (x, y) => x >= y,
};

/** The leaf evaluators over a column's raw data; every mask is a `BitArray` of the column's length. */
export class ColumnEvaluator {
  static kindOf(col: MaskColumnLike): FilterKind {
    return kindOf({name: col.name, type: col.type});
  }

  /** The rows holding null — never a bool (the platform has no null slot there). */
  static nulls(col: MaskColumnLike): BitArray {
    const kind = ColumnEvaluator.kindOf(col);
    if (kind === KIND.BOOL)
      return new BitArray(col.length);
    if (kind === KIND.BIG_INT || kind === KIND.STRING_LIST)
      return BitArray.create(col.length, (i) => col.get!(i) == null);
    const raw = col.getRawData();
    if (kind === KIND.INT)
      return BitArray.create(col.length, (i) => raw[i] === INT_NULL);
    if (kind === KIND.FLOAT || kind === KIND.DATE_TIME)
      return BitArray.create(col.length, (i) => isFloatNull(raw[i]));
    const categories = col.categories ?? [];
    return BitArray.create(col.length, (i) => (categories[raw[i]] ?? '') === '');
  }

  /** The rows whose non-null cell passes `test` against the condition's value converted to the
   * column's kind (a list element-wise); nulls never pass — the SQL reading every other target
   * has. Strings are tested once per category, bools once per side; a string list is tested as
   * its elements joined by newlines, so `like` reads "some element contains" and `!like` "none does". */
  static where(col: MaskColumnLike, c: FilterCondition, ctx: {now: Date},
    test: (cell: MaskCell, value: any) => boolean): BitArray {
    const kind = ColumnEvaluator.kindOf(col);
    const raw = kind === KIND.BIG_INT || kind === KIND.STRING_LIST ? null : col.getRawData();
    const f32 = raw instanceof Float32Array;
    const convert = (v: FilterScalar): MaskCell => ColumnEvaluator._cell(v, kind, ctx, f32);
    const value = Array.isArray(c.value) ? c.value.map(convert) : convert(c.value as FilterScalar);
    switch (kind) {
      case KIND.BOOL: {
        const hit = [test(false, value), test(true, value)];
        return BitArray.create(col.length, (i) => hit[(raw![i >>> 5] >>> (i & 31)) & 1]);
      }
      case KIND.BIG_INT:
        return BitArray.create(col.length, (i) => {
          const x = col.get!(i);
          return x != null && test(BigInt(x as bigint | number | string), value);
        });
      case KIND.STRING_LIST:
        return BitArray.create(col.length, (i) => {
          const x = col.get!(i);
          return Array.isArray(x) && test(x.join('\n'), value);
        });
      case KIND.INT:
        return BitArray.create(col.length, (i) => raw![i] !== INT_NULL && test(raw![i], value));
      case KIND.FLOAT: case KIND.DATE_TIME:
        return BitArray.create(col.length, (i) => !isFloatNull(raw![i]) && test(raw![i], value));
      default: {
        const hit = (col.categories ?? []).map((s) => s !== '' && test(s, value));
        return BitArray.create(col.length, (i) => hit[raw![i]] === true);
      }
    }
  }

  /** Column against column, row-wise, with one of the six comparators: a null on either side
   * never passes; dates compare by their time value, strings by their text. */
  static compare(left: MaskColumnLike, right: MaskColumnLike, operator: string): BitArray {
    const test = COMPARE[operator];
    if (test === undefined)
      throw new FilterError(`Operator "${operator}" does not accept a column`);
    const a = ColumnEvaluator.cells(left);
    const b = ColumnEvaluator.cells(right);
    return BitArray.create(left.length, (i) => {
      const x = a(i);
      const y = b(i);
      return x !== null && y !== null && test(x, y);
    });
  }

  /** A row's cell converted to the column's kind, null for a null cell. */
  static cells(col: MaskColumnLike): (i: number) => MaskCell | null {
    const kind = ColumnEvaluator.kindOf(col);
    if (kind === KIND.BIG_INT) {
      return (i) => {
        const x = col.get!(i);
        return x == null ? null : BigInt(x as bigint | number | string);
      };
    }
    if (kind === KIND.STRING_LIST) {
      return (i) => {
        const x = col.get!(i);
        return Array.isArray(x) ? x.join('\n') : null;
      };
    }
    const raw = col.getRawData();
    switch (kind) {
      case KIND.BOOL: return (i) => ((raw[i >>> 5] >>> (i & 31)) & 1) === 1;
      case KIND.INT: return (i) => raw[i] === INT_NULL ? null : raw[i];
      case KIND.FLOAT: case KIND.DATE_TIME: return (i) => isFloatNull(raw[i]) ? null : raw[i];
      default: {
        const categories = col.categories ?? [];
        return (i) => {
          const s = categories[raw[i]] ?? '';
          return s === '' ? null : s;
        };
      }
    }
  }

  /** A raw SQL LIKE pattern (`%`, `_`, backslash escapes) as an anchored case-insensitive RegExp. */
  static likeRegExp(pattern: string): RegExp {
    let source = '';
    for (let i = 0; i < pattern.length; i++) {
      const ch = pattern[i];
      if (ch === '\\' && i + 1 < pattern.length)
        source += escapeRegExp(pattern[++i]);
      else
        source += ch === '%' ? '[^]*' : ch === '_' ? '[^]' : escapeRegExp(ch);
    }
    return new RegExp(`^${source}$`, 'i');
  }

  private static _cell(v: FilterScalar, kind: FilterKind, ctx: {now: Date}, f32: boolean): MaskCell {
    switch (kind) {
      case KIND.INT: case KIND.FLOAT: {
        const n = typeof v === 'number' ? v : Number(v);
        return f32 ? Math.fround(n) : n;
      }
      case KIND.DATE_TIME: {
        const ms = v instanceof Date ? v.getTime() : isSpan(v) ?
          resolveSpan(v.span, ctx.now).getTime() : typeof v === 'number' ? v : Date.parse(String(v));
        return ms * 1000;
      }
      case KIND.BIG_INT:
        return BigInt(v as number | string);
      case KIND.BOOL:
        return v === true || v === 'true';
      default:
        return String(domainValue(v, ctx));
    }
  }
}

/** Evaluates the tree over the frame's columns: sync `mask` operators and async `bitset` ones
 * per leaf, a column reference row-wise against the other column, bitwise and/or/not per
 * group, an empty group all-true. Aborts between leaves. Throws `FilterError` for an unknown
 * column, an unbound `$param`, an operator without a DataFrame form, or a mask of the wrong length. */
export async function toMask(frame: MaskFrameLike, root: FilterGroup,
  options: {signal?: AbortSignal, now?: Date} = {}): Promise<BitArray> {
  const ctx = {now: options.now ?? new Date()};
  const signal = options.signal ?? NEVER;
  const group = async (g: FilterGroup): Promise<BitArray> => {
    let acc: BitArray | null = null;
    for (const n of g.nodes) {
      signal.throwIfAborted();
      const m = isGroup(n) ? await group(n) : await leafMask(frame, n, ctx, signal);
      // the clone keeps a registered operator's own mask untouched
      acc = acc === null ? m.clone() : g.op === 'and' ? acc.and(m) : acc.or(m);
    }
    const result = acc ?? new BitArray(frame.rowCount, true);
    return g.not === true ? result.invert() : result;
  };
  return group(root);
}

/** The rows a SQL CHECK constraint admits — the same tree read with three-valued logic, as
 * Postgres reads it (`domain_manifest_rules.dart`): a comparison against a null cell is UNKNOWN,
 * `and`/`or`/`not` are Kleene's, and a row is refused only where the tree is FALSE. A WHERE keeps
 * the two-valued reading of {@link toMask}, where UNKNOWN does not pass. */
export async function toCheckMask(frame: MaskFrameLike, root: FilterGroup,
  options: {signal?: AbortSignal, now?: Date} = {}): Promise<BitArray> {
  const ctx = {now: options.now ?? new Date()};
  const signal = options.signal ?? NEVER;
  // [true, unknown] per node; a row neither holds is false
  type Verdict = [BitArray, BitArray];
  const unknown = (c: FilterCondition): BitArray => {
    if (c.operator === 'is null' || c.operator === 'is not null')
      return new BitArray(frame.rowCount);
    const nulls = ColumnEvaluator.nulls(columnOf(frame, c.property, c));
    return isColumnRef(c.value) ? nulls.or(ColumnEvaluator.nulls(columnOf(frame, c.value.column, c))) : nulls;
  };
  const group = async (g: FilterGroup): Promise<Verdict> => {
    let acc: Verdict | null = null;
    for (const n of g.nodes) {
      signal.throwIfAborted();
      const v: Verdict = isGroup(n) ? await group(n) :
        [(await leafMask(frame, n, ctx, signal)).clone(), unknown(n)];
      if (acc === null) {
        acc = v;
        continue;
      }
      // and: unknown wherever neither side is false; or: wherever neither side is true
      acc = g.op === 'and' ?
        [acc[0].clone().and(v[0]),
          acc[1].clone().or(v[1]).and(acc[0].clone().or(acc[1])).and(v[0].clone().or(v[1]))] :
        [acc[0].clone().or(v[0]), acc[1].clone().or(v[1]).andNot(acc[0]).andNot(v[0])];
    }
    const result: Verdict = acc ?? [new BitArray(frame.rowCount, true), new BitArray(frame.rowCount)];
    return g.not === true ? [result[0].clone().or(result[1]).invert(), result[1]] : result;
  };
  const [truth, unsure] = await group(root);
  return truth.or(unsure);
}

function columnOf(frame: MaskFrameLike, name: string, c: FilterCondition): MaskColumnLike {
  const col = frame.column(name);
  if (!col) {
    const message = `Unknown column "${name}"`;
    throw new FilterError(message, [{nodeId: c.id, code: 'unknown-property', message}]);
  }
  return col;
}

async function leafMask(frame: MaskFrameLike, c: FilterCondition, ctx: {now: Date},
  signal: AbortSignal): Promise<BitArray> {
  const col = columnOf(frame, c.property, c);
  const param = (Array.isArray(c.value) ? c.value : [c.value]).find(isParam);
  if (param !== undefined) {
    const message = `Unbound parameter "$${param.param}"`;
    throw new FilterError(message, [{nodeId: c.id, code: 'not-expressible', message}]);
  }
  if (isColumnRef(c.value))
    return ColumnEvaluator.compare(col, columnOf(frame, c.value.column, c), c.operator);
  const op = operators.get(c.operator, {name: col.name, type: col.type, semType: col.semType});
  const m = op?.mask ? op.mask(col, c, ctx) : op?.bitset ? await op.bitset(col, c, signal) : null;
  if (m === null) {
    const message = `"${c.property} ${c.operator}" cannot be evaluated on a DataFrame`;
    throw new FilterError(message, [{nodeId: c.id, code: 'not-expressible', message}]);
  }
  if (m.length !== frame.rowCount) {
    const message = `"${c.property} ${c.operator}" returned a mask of ${m.length} bits for ${frame.rowCount} rows`;
    throw new FilterError(message, [{nodeId: c.id, code: 'evaluation', message}]);
  }
  return m;
}
