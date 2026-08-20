/* The DD5 DataFrame binding surface over a (possibly repointing) frame signal — one
   implementation serves a func-source output and a constant table-source. Platform-free:
   structural DataFrameLike; DG.DataFrame satisfies it. */
import {signal, computed, batch, Signal, ReadonlySignal} from '../core/signals.js';
import type {Scope} from '../core/scope.js';
import type {ObservableLike} from '../core/widget-like.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';

export interface ColumnLike {
  name: string;
  type: string;
  semType?: string | null;
}

export interface DataFrameLike {
  rowCount: number;
  currentRowIdx: number;
  get(column: string, row: number): unknown;
  set(column: string, row: number, value: unknown): void;
  columns: {names(): string[]; byName(name: string): ColumnLike | null};
  /** The real BitSets, passed as-is. */
  selection: unknown;
  filter: unknown;
  onCurrentRowChanged: ObservableLike<unknown>;
  onValuesChanged: ObservableLike<unknown>;
  onSelectionChanged: ObservableLike<unknown>;
  onFilterChanged: ObservableLike<unknown>;
  onColumnsChanged: ObservableLike<unknown>;
}

const VALUE = Object.getOwnPropertyDescriptor(Signal.prototype, 'value')!;
const LIVE_PROTO = Object.create(Signal.prototype, {value: {get: VALUE.get}});

/** A read-only signal holding a live platform object (a BitSet): dependents must re-fire when
 * the object mutates even though its identity never changes, which both the engine's set-time
 * equality guard and its batch set-and-restore reconciliation would suppress. So a refresh
 * silently stomps the vendored engine's `_value` with a throwaway sentinel (no notification)
 * and then performs one real set — exactly one flush, dependents never observe the sentinel —
 * while the interposed getter-only `value` prototype keeps `isWritableSignal` false. */
function liveSignal(): {sig: Signal<unknown>, set(v: unknown): void} {
  const sig = signal<unknown>(undefined);
  Object.setPrototypeOf(sig, LIVE_PROTO);
  return {
    sig,
    set: (v) => {
      (sig as unknown as {_value: unknown})._value = {};
      VALUE.set!.call(sig, v);
    },
  };
}

/** The {@link BindSource} over a DataFrame signal. Subscriptions and signals live on `scope`.
 * Lazy: a `currentRow.<col>` signal is created on first bindStep and cached; a repoint (the
 * frame signal changing identity) resubscribes events and refreshes every cached signal.
 * A column step always answers a signal — an absent frame or column reads `undefined` and
 * refuses writes, so bound inputs show empty instead of breaking (the containment posture). */
export function dfBindings(df: ReadonlySignal<DataFrameLike | undefined>, scope: Scope): BindSource {
  const columnSigs = new Map<string, Signal<unknown>>();
  const rowIdx = signal(-1);
  const dataVersion = signal(0);
  const columnsVersion = signal(0);
  const selection = liveSignal();
  const filter = liveSignal();

  const dfIdentity = computed(() => df.value) as unknown as Signal<unknown>;
  const rowCount = computed(() => {
    dataVersion.value;
    return df.value?.rowCount ?? 0;
  }) as unknown as Signal<unknown>;
  const columns = computed(() => {
    columnsVersion.value;
    return df.value?.columns.names() ?? [];
  }) as unknown as Signal<unknown>;

  /** The row the form is on, or -1 when there is none to read. `currentRowIdx` is a two-way
   * binding, so it holds whatever was typed into the input driving it — an index past the last row
   * is a row that does not exist, whether the frame clamped it, refused it or kept it. */
  const rowAt = (d: DataFrameLike | undefined): number => {
    const at = rowIdx.peek();
    return d !== undefined && at >= 0 && at < d.rowCount ? at : -1;
  };
  const readColumn = (name: string): unknown => {
    const d = df.peek();
    const at = rowAt(d);
    return at >= 0 && d!.columns.byName(name) !== null ? d!.get(name, at) : undefined;
  };
  const refreshColumns = () => {
    for (const [name, sig] of columnSigs)
      sig.value = readColumn(name);
  };
  const bump = (sig: Signal<number>) => sig.value = sig.peek() + 1;

  let subs: {unsubscribe(): void}[] = [];
  const drop = () => {
    for (const s of subs)
      s.unsubscribe();
    subs = [];
  };
  scope.effect(() => {
    const d = df.value;
    drop();
    if (d !== undefined) {
      subs.push(d.onCurrentRowChanged.subscribe(() => batch(() => {
        rowIdx.value = d.currentRowIdx;
        refreshColumns();
      })));
      subs.push(d.onValuesChanged.subscribe(() => batch(() => {
        bump(dataVersion);
        refreshColumns();
      })));
      subs.push(d.onSelectionChanged.subscribe(() => selection.set(d.selection)));
      subs.push(d.onFilterChanged.subscribe(() => filter.set(d.filter)));
      subs.push(d.onColumnsChanged.subscribe(() => bump(columnsVersion)));
    }
    batch(() => {
      rowIdx.value = d?.currentRowIdx ?? -1;
      refreshColumns();
      bump(dataVersion);
      bump(columnsVersion);
      selection.set(d?.selection);
      filter.set(d?.filter);
    });
  });
  scope.own(drop);

  scope.effect(() => {
    const v = rowIdx.value;
    const d = df.peek();
    if (d === undefined)
      return;
    // first, and unconditionally: the fields read the row that was ASKED for, so an index the
    // frame cannot hold empties them instead of keeping the last row that did exist — which reads
    // as real data. Idempotent, since the same value written to a signal notifies nobody
    refreshColumns();
    // and the frame is only moved where it can go: the platform throws on an out-of-bounds row,
    // and the throw used to take the refresh above down with it
    if (d.currentRowIdx !== v && v >= -1 && v < d.rowCount)
      d.currentRowIdx = v;
  });

  const columnSignal = (name: string): Signal<unknown> => {
    let sig = columnSigs.get(name);
    if (sig === undefined) {
      const created = signal(readColumn(name));
      columnSigs.set(name, created);
      scope.effect(() => {
        const v = created.value;
        const d = df.peek();
        const at = rowAt(d);
        if (at >= 0 && d!.columns.byName(name) !== null && d!.get(name, at) !== v)
          d!.set(name, at, v);
      });
      sig = created;
    }
    return sig;
  };

  const currentRow: BindSource = {
    bindStep: (name) => name === '' ? null : columnSignal(name),
    bindProps: () => {
      const d = df.peek();
      if (d === undefined)
        return [];
      return d.columns.names().map((name): BindProp => {
        const col = d.columns.byName(name);
        return {name, type: col?.type ?? null, semType: col?.semType ?? null, writable: true};
      });
    },
  };

  const steps: Record<string, Signal<unknown> | BindSource> = {
    df: dfIdentity,
    currentRowIdx: rowIdx as unknown as Signal<unknown>,
    currentRow,
    selection: selection.sig,
    filter: filter.sig,
    rowCount,
    columns,
  };
  return {
    bindStep: (name) => steps[name === '' ? 'df' : name] ?? null,
    bindProps: (): BindProp[] => [
      {name: 'df', type: 'dataframe'},
      {name: 'currentRowIdx', type: 'int', writable: true},
      {name: 'currentRow', type: 'object', walkable: true},
      {name: 'selection', type: 'object'},
      {name: 'filter', type: 'object'},
      {name: 'rowCount', type: 'int'},
      {name: 'columns', type: 'string_list'},
    ],
  };
}
