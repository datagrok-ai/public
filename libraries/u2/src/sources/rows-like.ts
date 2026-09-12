/* The collection protocol lists, cards and pickers read (GOAL ruling 2): items by key, so a plain
   array and a DataFrame feed a control alike — `arrayRows` here, `FrameRows` in df-rows.ts — and
   the row conventions every layer shares: the editor's service columns and the draft key. */
import {signal, Signal, ReadonlySignal} from '../core/signals.js';

/** One row of a domain table as controls see it: its columns by name, `id` the key. Live over
 * the frame the source holds. A draft carries a temporary key until it is saved. */
export interface RowView {
  readonly id: string;
  [column: string]: unknown;
}

export interface RowsLike<T> {
  readonly items: ReadonlySignal<readonly T[]>;
  keyOf(item: T): string;
  byKey(key: string): T | undefined;
}

/** The one place the `~` conventions live: the editor's row-state column, the service-column
 * rule (read by name, never enumerated — H7) and the key a draft carries until it is saved. */
export class Rows {
  /** The editor's row-state service column: `''`, `'new'`, `'modified'` or `'deleted'`. */
  static readonly STATE = '~state';
  /** A draft has no id yet: it is keyed by its row index, `~row:<index>`. */
  static readonly DRAFT_PREFIX = '~row:';

  static isService(column: string): boolean {
    return column.startsWith('~');
  }

  static isDraft(x: string | RowView): boolean {
    return (typeof x === 'string' ? x : x.id).startsWith(Rows.DRAFT_PREFIX);
  }

  static draftKey(index: number): string {
    return `${Rows.DRAFT_PREFIX}${index}`;
  }

  static draftIndex(key: string): number {
    return Number(key.slice(Rows.DRAFT_PREFIX.length));
  }

  /** The frame row a draft key names (-1 past the frame); null for a key that is not a draft's. */
  static draftRow(key: string, rowCount: number): number | null {
    if (!Rows.isDraft(key))
      return null;
    const at = Rows.draftIndex(key);
    return at < rowCount ? at : -1;
  }
}

/** Rows over an array, or over a signal of one — what a non-EMS list is fed with. */
export function arrayRows<T>(items: readonly T[] | ReadonlySignal<readonly T[]>,
  keyOf: (item: T) => string): RowsLike<T> {
  const source: ReadonlySignal<readonly T[]> = items instanceof Signal ?
    items as ReadonlySignal<readonly T[]> : signal(items as readonly T[]);
  let indexed: readonly T[] | undefined;
  let index = new Map<string, T>();
  return {
    items: source,
    keyOf,
    byKey: (key) => {
      const now = source.peek();
      if (now !== indexed) {
        index = new Map(now.map((item) => [keyOf(item), item]));
        indexed = now;
      }
      return index.get(key);
    },
  };
}
