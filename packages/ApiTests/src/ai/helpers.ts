import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {take} from 'rxjs/operators';
import {expect, expectObject} from '@datagrok-libraries/test/src/test';

// Shared helpers for AI tests. Keep the surface narrow — only patterns that
// repeat across at least 3 files. Tests that need bespoke logic should still
// inline it.

export const demog = (n: number = 50): DG.DataFrame => grok.data.demo.demog(n);

// Shortcut for `v.getOptions(true).look`.
export const look = (v: DG.Viewer): {[k: string]: any} => v.getOptions(true).look;

// Assert every key in `partial` is present and equal in `v.getOptions(true).look`.
export function expectLook(v: DG.Viewer, partial: {[k: string]: any}): void {
  expectObject(look(v), partial);
}

// Assert every key in `partial` is present and equal in both `v.props` and look.
export function expectPropAndLook(v: DG.Viewer, partial: {[k: string]: any}): void {
  for (const k of Object.keys(partial))
    expect(v.props[k], partial[k]);
  expectObject(look(v), partial);
}

// setOptions then expect every set key shows up in look.
export function expectRoundTrip(v: DG.Viewer, opts: {[k: string]: any}): void {
  v.setOptions(opts);
  expectObject(look(v), opts);
}

// setOptions then expect every set key shows up in both props and look.
export function expectRoundTripPropAndLook(v: DG.Viewer, opts: {[k: string]: any}): void {
  v.setOptions(opts);
  expectPropAndLook(v, opts);
}

// Round-trip a plain boolean get/set accessor (false -> true -> restore original).
// For widget accessors like `acc.autoHideTabHeader`; use expectBoolToggle for `v.props[name]`.
export function expectBoolGetSet(obj: {[k: string]: any}, name: string): void {
  const original = obj[name];
  obj[name] = false;
  expect(obj[name], false);
  obj[name] = true;
  expect(obj[name], true);
  obj[name] = original;
  expect(obj[name], original);
}

// Toggle a boolean prop via `v.props[name] = ...` and check props+look mirror.
export function expectBoolToggle(v: DG.Viewer, name: string, sequence: boolean[] = [true, false, true]): void {
  for (const b of sequence) {
    (v.props as any)[name] = b;
    expect((v.props as any)[name], b);
    expect(look(v)[name], b);
  }
}

// Returns true iff the action did not throw.
export function noThrow(fn: () => void): boolean {
  try {
    fn();
    return true;
  } catch (_e) {
    return false;
  }
}

// Assert that `fn` does not throw.
export function expectNoThrow(fn: () => void): void {
  expect(noThrow(fn), true);
}

// Find a property descriptor by name. Returns null if missing.
export function findProp(v: DG.Viewer, name: string): DG.Property | null {
  for (const p of v.getProperties() as DG.Property[]) {
    if (p.name === name)
      return p;
  }
  return null;
}

// Assert that `propName` exposes the given choices. If the property descriptor
// is missing (some props aren't surfaced on every build), returns null so the
// caller can apply a fallback.
export function expectChoices(v: DG.Viewer, propName: string, required: string[]): DG.Property | null {
  const p = findProp(v, propName);
  if (p == null)
    return null;
  expect(Array.isArray(p.choices), true);
  for (const c of required)
    expect(p.choices.indexOf(c) >= 0, true);
  return p;
}

// Subscribe to each stream and return a cleanup function. Asserts shape.
export function subscribeAll(streams: Observable<any>[]): () => void {
  const subs: Subscription[] = [];
  for (const s of streams) {
    expect(s != null, true);
    expect(typeof s.subscribe, 'function');
    const sub = s.subscribe(() => {});
    expect(typeof sub.unsubscribe, 'function');
    subs.push(sub);
  }
  return () => {
    for (const sub of subs)
      sub.unsubscribe();
  };
}

// Run a body with `addTableView` and close in finally.
export async function withTableView(df: DG.DataFrame, body: (tv: DG.TableView) => Promise<void> | void): Promise<void> {
  const tv = grok.shell.addTableView(df);
  try {
    await body(tv);
  } finally {
    tv.close();
  }
}

// Attach the given viewer type to a fresh table view and run a body with it.
export async function withAttachedViewer<T extends DG.Viewer>(
  df: DG.DataFrame, type: string, opts: {[k: string]: any},
  body: (v: T, tv: DG.TableView) => Promise<void> | void): Promise<void> {
  await withTableView(df, async (tv) => {
    const v = tv.addViewer(type, opts) as T;
    await body(v, tv);
  });
}

// Race a one-shot Observable against a timeout. Asserts the event fires.
export async function expectFiresWithin(
  stream: Observable<any>, trigger: () => void, ms: number = 1500): Promise<void> {
  let sub: Subscription | undefined;
  try {
    const fired: Promise<boolean> = new Promise((resolve) => {
      sub = stream.pipe(take(1)).subscribe(() => resolve(true));
    });
    const timeout: Promise<boolean> = (async () => {
      await DG.delay(ms);
      return false;
    })();
    trigger();
    expect(await Promise.race([fired, timeout]), true);
  } finally {
    if (sub != null)
      sub.unsubscribe();
  }
}

// Assert a value is null, undefined, or an empty string — the platform's "cleared" shape.
export function expectCleared(value: any): void {
  expect(value == null || value === '', true);
}

// Wrapper around DG.delay so call sites can read as `await wait(300)`.
export const wait = (ms: number = 300): Promise<void> => DG.delay(ms);

// Build a small DataFrame from named lists. Each entry: [name, dartType, values].
export function df(cols: Array<[string, string, any[]]>): DG.DataFrame {
  return DG.DataFrame.fromColumns(cols.map(([n, t, v]) => DG.Column.fromList(t as DG.ColumnType, n, v)));
}
