import {FilterError, isGroup, isSpan} from './model.js';
import type {FilterCondition, FilterGroup, FilterKind, FilterScalar} from './model.js';
import {domainValue, kindOf} from './kinds.js';
import {resolveSpan} from '../span.js';
// cycle with operators.ts is function-scoped only: neither module reads the other at load time
import {operators} from './operators.js';

/** LSB-first words, bit i at `bits[i >>> 5] & (1 << (i & 31))` — the layout `DG.BitSet.fromBytes` reads. */
export interface Mask { bits: Uint32Array; length: number }

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

export class Masks {
  static create(length: number, fill: boolean = false): Mask {
    const bits = new Uint32Array((length + 31) >>> 5);
    if (fill)
      bits.fill(0xffffffff);
    return Masks._trim({bits, length});
  }

  static and(a: Mask, b: Mask): Mask {
    const bits = new Uint32Array(a.bits.length);
    for (let i = 0; i < bits.length; i++)
      bits[i] = a.bits[i] & b.bits[i];
    return {bits, length: a.length};
  }

  static or(a: Mask, b: Mask): Mask {
    const bits = new Uint32Array(a.bits.length);
    for (let i = 0; i < bits.length; i++)
      bits[i] = a.bits[i] | b.bits[i];
    return {bits, length: a.length};
  }

  /** The complement, with the tail beyond `length` cleared. */
  static not(a: Mask): Mask {
    const bits = new Uint32Array(a.bits.length);
    for (let i = 0; i < bits.length; i++)
      bits[i] = ~a.bits[i];
    return Masks._trim({bits, length: a.length});
  }

  /** A mask over a provider's LSB-first words — a `Uint32Array` or its `ArrayBuffer` — copied to the
   * word count of `length`, the tail beyond it cleared. */
  static from(bits: Uint32Array | ArrayBuffer, length: number): Mask {
    const words = bits instanceof Uint32Array ? bits : new Uint32Array(bits);
    const mask = Masks.create(length);
    mask.bits.set(words.subarray(0, mask.bits.length));
    return Masks._trim(mask);
  }

  static fromPredicate(length: number, test: (i: number) => boolean): Mask {
    const mask = Masks.create(length);
    for (let i = 0; i < length; i++) {
      if (test(i))
        mask.bits[i >>> 5] |= 1 << (i & 31);
    }
    return mask;
  }

  static get(mask: Mask, i: number): boolean {
    return (mask.bits[i >>> 5] & (1 << (i & 31))) !== 0;
  }

  static toIndexes(mask: Mask): number[] {
    const indexes: number[] = [];
    for (let i = 0; i < mask.length; i++) {
      if (Masks.get(mask, i))
        indexes.push(i);
    }
    return indexes;
  }

  static kindOf(col: MaskColumnLike): FilterKind {
    return kindOf({name: col.name, type: col.type});
  }

  /** The rows holding null — never a bool (the platform has no null slot there). */
  static nulls(col: MaskColumnLike): Mask {
    const kind = Masks.kindOf(col);
    if (kind === 'bool')
      return Masks.create(col.length);
    if (kind === 'bigint' || kind === 'string_list')
      return Masks.fromPredicate(col.length, (i) => col.get!(i) == null);
    const raw = col.getRawData();
    if (kind === 'int')
      return Masks.fromPredicate(col.length, (i) => raw[i] === INT_NULL);
    if (kind === 'float' || kind === 'datetime')
      return Masks.fromPredicate(col.length, (i) => isFloatNull(raw[i]));
    const categories = col.categories ?? [];
    return Masks.fromPredicate(col.length, (i) => (categories[raw[i]] ?? '') === '');
  }

  /** The rows whose non-null cell passes `test` against the condition's value converted to the
   * column's kind (a list element-wise); nulls never pass — the SQL reading every other target
   * has. Strings are tested once per category, bools once per side; a string list is tested as
   * its elements joined by newlines, so `like` reads "some element contains" and `!like` "none does". */
  static where(col: MaskColumnLike, c: FilterCondition, ctx: {now: Date},
    test: (cell: MaskCell, value: any) => boolean): Mask {
    const kind = Masks.kindOf(col);
    const raw = kind === 'bigint' || kind === 'string_list' ? null : col.getRawData();
    const f32 = raw instanceof Float32Array;
    const convert = (v: FilterScalar): MaskCell => Masks._cell(v, kind, ctx, f32);
    const value = Array.isArray(c.value) ? c.value.map(convert) : convert(c.value as FilterScalar);
    switch (kind) {
      case 'bool': {
        const hit = [test(false, value), test(true, value)];
        return Masks.fromPredicate(col.length, (i) => hit[(raw![i >>> 5] >>> (i & 31)) & 1]);
      }
      case 'bigint':
        return Masks.fromPredicate(col.length, (i) => {
          const x = col.get!(i);
          return x != null && test(BigInt(x as bigint | number | string), value);
        });
      case 'string_list':
        return Masks.fromPredicate(col.length, (i) => {
          const x = col.get!(i);
          return Array.isArray(x) && test(x.join('\n'), value);
        });
      case 'int':
        return Masks.fromPredicate(col.length, (i) => raw![i] !== INT_NULL && test(raw![i], value));
      case 'float': case 'datetime':
        return Masks.fromPredicate(col.length, (i) => !isFloatNull(raw![i]) && test(raw![i], value));
      default: {
        const hit = (col.categories ?? []).map((s) => s !== '' && test(s, value));
        return Masks.fromPredicate(col.length, (i) => hit[raw![i]] === true);
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
      case 'int': case 'float': {
        const n = typeof v === 'number' ? v : Number(v);
        return f32 ? Math.fround(n) : n;
      }
      case 'datetime': {
        const ms = v instanceof Date ? v.getTime() : isSpan(v) ?
          resolveSpan(v.span, ctx.now).getTime() : typeof v === 'number' ? v : Date.parse(String(v));
        return ms * 1000;
      }
      case 'bigint':
        return BigInt(v as number | string);
      case 'bool':
        return v === true || v === 'true';
      default:
        return String(domainValue(v, ctx));
    }
  }

  private static _trim(mask: Mask): Mask {
    const tail = mask.length & 31;
    if (tail !== 0 && mask.bits.length > 0)
      mask.bits[mask.bits.length - 1] &= (1 << tail) - 1;
    return mask;
  }
}

/** Evaluates the tree over the frame's columns: sync `mask` operators and async `bitset` ones
 * per leaf, bitwise and/or/not per group, an empty group all-true. Aborts between leaves.
 * Throws `FilterError` for an unknown column or an operator without a DataFrame form. */
export async function toMask(frame: MaskFrameLike, root: FilterGroup,
  options: {signal?: AbortSignal, now?: Date} = {}): Promise<Mask> {
  const ctx = {now: options.now ?? new Date()};
  const signal = options.signal ?? NEVER;
  const leaf = async (c: FilterCondition): Promise<Mask> => {
    const col = frame.column(c.property);
    if (!col) {
      const message = `Unknown column "${c.property}"`;
      throw new FilterError(message, [{nodeId: c.id, code: 'unknown-property', message}]);
    }
    const op = operators.get(c.operator, {name: col.name, type: col.type, semType: col.semType});
    if (op?.mask)
      return op.mask(col, c, ctx);
    if (op?.bitset)
      return await op.bitset(col, c, signal);
    const message = `"${c.property} ${c.operator}" cannot be evaluated on a DataFrame`;
    throw new FilterError(message, [{nodeId: c.id, code: 'not-expressible', message}]);
  };
  const group = async (g: FilterGroup): Promise<Mask> => {
    let acc: Mask | null = null;
    for (const n of g.nodes) {
      signal.throwIfAborted();
      const m = isGroup(n) ? await group(n) : await leaf(n);
      acc = acc === null ? m : g.op === 'and' ? Masks.and(acc, m) : Masks.or(acc, m);
    }
    const result = acc ?? Masks.create(frame.rowCount, true);
    return g.not === true ? Masks.not(result) : result;
  };
  return group(root);
}
