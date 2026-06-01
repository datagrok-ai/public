import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {take} from 'rxjs/operators';
import {expect, expectObject, awaitCheck} from '@datagrok-libraries/test/src/test';

// Shared helpers for AI tests.

export const demog = (n: number = 50): DG.DataFrame => grok.data.demo.demog(n);

export const look = (v: DG.Viewer): {[k: string]: any} => v.getOptions(true).look;

export function expectLook(v: DG.Viewer, partial: {[k: string]: any}): void {
  expectObject(look(v), partial);
}

export function expectPropAndLook(v: DG.Viewer, partial: {[k: string]: any}): void {
  for (const k of Object.keys(partial))
    expect(v.props[k], partial[k]);
  expectObject(look(v), partial);
}

export function expectRoundTrip(v: DG.Viewer, opts: {[k: string]: any}): void {
  v.setOptions(opts);
  expectObject(look(v), opts);
}

export function expectRoundTripPropAndLook(v: DG.Viewer, opts: {[k: string]: any}): void {
  v.setOptions(opts);
  expectPropAndLook(v, opts);
}

export function expectBoolGetSet(obj: {[k: string]: any}, name: string): void {
  const original = obj[name];
  obj[name] = false;
  expect(obj[name], false);
  obj[name] = true;
  expect(obj[name], true);
  obj[name] = original;
  expect(obj[name], original);
}

export function expectBoolToggle(v: DG.Viewer, name: string, sequence: boolean[] = [true, false, true]): void {
  for (const b of sequence) {
    (v.props as any)[name] = b;
    expect((v.props as any)[name], b);
    expect(look(v)[name], b);
  }
}

export function noThrow(fn: () => void): boolean {
  try {
    fn();
    return true;
  } catch (_e) {
    return false;
  }
}

export function expectNoThrow(fn: () => void): void {
  expect(noThrow(fn), true);
}

export function findProp(v: DG.Viewer, name: string): DG.Property | null {
  for (const p of v.getProperties() as DG.Property[]) {
    if (p.name === name)
      return p;
  }
  return null;
}

// Returns null if the property descriptor is missing, so the caller can fall back.
export function expectChoices(v: DG.Viewer, propName: string, required: string[]): DG.Property | null {
  const p = findProp(v, propName);
  if (p == null)
    return null;
  expect(Array.isArray(p.choices), true);
  for (const c of required)
    expect(p.choices.indexOf(c) >= 0, true);
  return p;
}

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

export async function withTableView(df: DG.DataFrame, body: (tv: DG.TableView) => Promise<void> | void): Promise<void> {
  const tv = grok.shell.addTableView(df);
  try {
    await body(tv);
  } finally {
    tv.close();
  }
}

export async function withAttachedViewer<T extends DG.Viewer>(
  df: DG.DataFrame, type: string, opts: {[k: string]: any},
  body: (v: T, tv: DG.TableView) => Promise<void> | void): Promise<void> {
  await withTableView(df, async (tv) => {
    const v = tv.addViewer(type, opts) as T;
    await body(v, tv);
  });
}

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

export function expectCleared(value: any): void {
  expect(value == null || value === '', true);
}

export const wait = (ms: number = 300): Promise<void> => DG.delay(ms);

// Poll until cond() holds; transient throws from not-yet-rendered getters are swallowed.
export async function until(cond: () => boolean, ms: number = 2000): Promise<void> {
  await awaitCheck(() => {
    try {
      return cond();
    } catch (_e) {
      return false;
    }
  }, `until: condition not met within ${ms}ms`, ms, 25);
}

// Build a DataFrame from named lists. Each entry: [name, dartType, values].
export function df(cols: Array<[string, string, any[]]>): DG.DataFrame {
  return DG.DataFrame.fromColumns(cols.map(([n, t, v]) => DG.Column.fromList(t as DG.ColumnType, n, v)));
}
