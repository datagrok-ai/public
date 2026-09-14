/* The collection protocol lists, cards and pickers read (GOAL ruling 2): items by key, so a plain
   array and a DataFrame feed a control alike — `arrayRows` here, `FrameRows` in df-rows.ts — and
   the row conventions every layer shares: the editor's service columns and the draft id. */
import {signal, Signal, ReadonlySignal} from '../core/signals.js';

/** What a typed row type must carry: the key. An app's own row type (`IssueRow` from the
 * generated `db.ts`) satisfies it as an interface, so no index signature is asked for. */
export interface DomainRowLike {
  readonly id: string;
}

/** One row of a domain table as controls see it: its columns by name, `id` the key. Live over
 * the frame the source holds. A draft carries a temporary id until it is saved. `TRow` types the
 * columns; the bare `RowView` reads any column as unknown. */
export type RowView<TRow extends DomainRowLike = DomainRowLike> =
  TRow & {readonly id: string, [column: string]: unknown};

/** Column values a draft starts from: the typed columns by name, any other as unknown. */
export type RowValues<TRow extends DomainRowLike = DomainRowLike> = Partial<RowView<TRow>>;

/** A column name: the keys of `TRow` are offered, any string is taken — a strict `keyof` would
 * make a typed table unassignable where an untyped one is expected. */
export type ColumnOf<TRow extends DomainRowLike> = (keyof TRow & string) | (string & {});

export interface RowsLike<T> {
  readonly items: ReadonlySignal<readonly T[]>;
  keyOf(item: T): string;
  byKey(key: string): T | undefined;
}

/** The one place the `~` conventions live: the editor's row-state column, the service-column
 * rule (read by name, never enumerated — H7) and the id a draft carries until it is saved. */
export class Rows {
  /** The editor's row-state service column: `''`, `'new'`, `'modified'` or `'deleted'`. */
  static readonly STATE = '~state';
  /** A draft's id, stamped into the `id` cell by the writer that adds it (the js-api editor, the
   * memory edit state): `~new:<uuid>`, so a child may reference its parent before either exists. */
  static readonly DRAFT_PREFIX = '~new:';

  static isService(column: string): boolean {
    return column.startsWith('~');
  }

  static isDraft(x: string | RowView): boolean {
    return (typeof x === 'string' ? x : x.id).startsWith(Rows.DRAFT_PREFIX);
  }

  static draftId(): string {
    return `${Rows.DRAFT_PREFIX}${crypto.randomUUID()}`;
  }

  /** The key of a row that has no id cell at all (a non-EMS frame): its index, never a draft. */
  static unkeyed(index: number): string {
    return `~row:${index}`;
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
