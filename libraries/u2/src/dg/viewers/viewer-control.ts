/* The one construction path of a platform viewer as a u2 control (VP-8): the spec's `create` and
   the fluent `viewers.*` factories both end here. */
import * as DG from 'datagrok-api/dg';
import {Signal, signal, isWritableSignal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import type {Control} from '../../core/component.js';
import {adopt} from './adopt.js';

/** A viewer's options, each either a value or a signal the viewer follows — and writes back to,
 * where the signal is writable. */
export type Bindable<T> = {[K in keyof T]?: T[K] | ReadonlySignal<T[K]>};

/** The frame a viewer shows: a DataFrame, or a signal of one that may still read undefined. */
export type TableRef = DG.DataFrame | ReadonlySignal<DG.DataFrame | undefined>;

type Look = Record<string, unknown>;

/** The viewers mid-repoint — `dataFrame = df` then the literal look re-applied — whose property
 * events are the construction-time look coming back, not an edit (VP-23). */
export const REPOINTING = new WeakSet<DG.Viewer>();

/** `make` over the table — or an empty frame while the table signal reads undefined, so the viewer
 * shows its own empty state until the first frame repoints it — adopted, repointed whenever the
 * signal changes frame identity, signal options linked two-way where writable. */
export function viewerOf<V extends DG.Viewer, S>(make: (df: DG.DataFrame, look: Partial<S>) => V,
  table: TableRef, options?: Bindable<S>): V & Control {
  const look: Look = {};
  const literal: Look = {};
  const bound: [string, Signal<unknown>][] = [];
  for (const [name, value] of Object.entries(options ?? {})) {
    look[name] = value instanceof Signal ? value.peek() : value;
    if (value instanceof Signal)
      bound.push([name, value as Signal<unknown>]);
    else
      literal[name] = value;
  }
  return build(make, table, look, literal, bound);
}

/** The registry's `create`: `table` (a resolved signal, or absent), `look` (the object escape
 * hatch), every other key a look prop by `Property.name`; a named key wins over the same key in
 * `look`. A bound look prop arrives as a signal from the renderer, which bridges it itself through
 * `link` (VP-4) — here it is peeked for the initial look, never linked. Signals are peeled at the
 * top level only: a signal nested inside `look` is not supported. */
export function viewerControl(type: string, props: Record<string, unknown>): DG.Viewer & Control {
  const {table, look, ...named} = props;
  const literal: Look = {...(look as Look | undefined)};
  const options: Look = {...literal};
  for (const [name, value] of Object.entries(named)) {
    options[name] = value instanceof Signal ? value.peek() : value;
    if (value instanceof Signal)
      delete literal[name];
    else
      literal[name] = value;
  }
  return build((df, l) => DG.Viewer.fromType(type, df, l), table as TableRef, options, literal, [], type);
}

/** The platform resolves column-naming look props against the frame the viewer is built over and
 * drops what the empty placeholder lacks, so a repoint re-applies the construction-time literal
 * look — never a bound option, whose link keeps it consistent. A type that cannot be built over
 * the empty placeholder at all (Summary's RangeError) fails with what to do about it, the
 * platform's error as the cause. */
function build<V extends DG.Viewer, S>(make: (df: DG.DataFrame, look: Partial<S>) => V, table: TableRef,
  look: Look, literal: Look, bound: [string, Signal<unknown>][], type = 'This viewer'): V & Control {
  const frame = (table instanceof Signal ? table : signal(table)) as ReadonlySignal<DG.DataFrame | undefined>;
  const initial = frame.peek();
  // a `filters` literal over the placeholder is one platform error per entry, and the panes come
  // twice once the frame arrives: the repoint's `setOptions(literal)` applies it over the real frame
  const {filters: _filters, ...placeholderLook} = look;
  let made: V;
  try {
    made = make(initial ?? DG.DataFrame.create(0), (initial == null ? placeholderLook : look) as Partial<S>);
  } catch (e) {
    if (initial != null)
      throw e;
    throw new Error(`${type} needs a table — bind \`table\` to a data source`, {cause: e});
  }
  const viewer = adopt(made);
  viewer.effect(() => {
    const df = frame.value;
    if (df == null || viewer.dataFrame?.dart === df.dart)
      return;
    REPOINTING.add(viewer);
    try {
      viewer.dataFrame = df;
      viewer.setOptions(literal);
    } finally {
      REPOINTING.delete(viewer);
    }
  });
  viewer.run(() => {
    for (const [name, source] of bound)
      viewer.link(name, source, isWritableSignal(source));
  });
  return viewer;
}
