import {BitArray} from 'datagrok-api/u2core';
import {FilterError, KIND, isGroup, isSpan} from './model.js';
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
 * per leaf, bitwise and/or/not per group, an empty group all-true. Aborts between leaves.
 * Throws `FilterError` for an unknown column, an operator without a DataFrame form, or a mask of the wrong length. */
export async function toMask(frame: MaskFrameLike, root: FilterGroup,
  options: {signal?: AbortSignal, now?: Date} = {}): Promise<BitArray> {
  const ctx = {now: options.now ?? new Date()};
  const signal = options.signal ?? NEVER;
  const leaf = async (c: FilterCondition): Promise<BitArray> => {
    const col = frame.column(c.property);
    if (!col) {
      const message = `Unknown column "${c.property}"`;
      throw new FilterError(message, [{nodeId: c.id, code: 'unknown-property', message}]);
    }
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
  };
  const group = async (g: FilterGroup): Promise<BitArray> => {
    let acc: BitArray | null = null;
    for (const n of g.nodes) {
      signal.throwIfAborted();
      const m = isGroup(n) ? await group(n) : await leaf(n);
      // the clone keeps a registered operator's own mask untouched
      acc = acc === null ? m.clone() : g.op === 'and' ? acc.and(m) : acc.or(m);
    }
    const result = acc ?? new BitArray(frame.rowCount, true);
    return g.not === true ? result.invert() : result;
  };
  return group(root);
}
