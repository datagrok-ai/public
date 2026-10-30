/* Wrapper over the vendored signals engine — the only file that imports it, so the
   engine stays swappable (TC39 Signals later). Raw `effect` is deliberately not part of
   the public surface: every effect must have an owner (Scope), or it leaks. */
import {Signal, batch as coreBatch, computed, untracked} from './signals-core.js';
import type {ReadonlySignal} from './signals-core.js';

export {signal, computed, untracked, Signal} from './signals-core.js';
export type {ReadonlySignal} from './signals-core.js';
export {effect as rawEffect} from './signals-core.js';

let depth = 0;
let after: (() => void)[] = [];

/** Runs `fn` with every signal write inside it as one change: the effects it sets off run once, when the outermost
 * batch ends — and then what {@link afterFlush} queued. */
export function batch<T>(fn: () => T): T {
  depth++;
  try {
    return coreBatch(fn);
  }
  finally {
    if (--depth === 0 && after.length > 0)
      flushAfter();
  }
}

function flushAfter(): void {
  const queue = after;
  after = [];
  let failed = false;
  let error: unknown;
  for (const fn of queue) {
    try {
      batch(fn);
    }
    catch (e) {
      if (!failed)
        error = e;
      failed = true;
    }
  }
  if (failed)
    throw error;
}

/** Runs `fn` once the outermost {@link batch} has ended and the effects it set off have run — what an effect builds
 * in answer to a write (a new row) exists then; outside a batch, now. u2 runs every handler in a batch. A throw is
 * thrown at whoever ended the batch: hand it a {@link Scope.bind} handler to have it reported instead.
 * @example
 * rows.value = [...rows.value, record];
 * DG.U2.afterFlush(() => list.querySelector('input')?.focus());
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-lens} */
export function afterFlush(fn: () => void): void {
  if (depth === 0)
    fn();
  else
    after.push(fn);
}

/** A writable signal over a part of `source`: reading it reads `get(source.value)` (tracked), writing it writes
 * `source.value = set(source.peek(), x)`. Writing the value it already reads changes nothing, so two of them wired
 * both ways settle at once. It is a {@link Signal}: it binds wherever a signal does.
 * @example
 * const person = DG.U2.signal({name: 'Ann', age: 41});
 * const name = DG.U2.lens(person, (p) => p.name, (p, name) => ({...p, name}));
 * name.value = 'Bo';   // person.value is {name: 'Bo', age: 41}
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-lens} */
export function lens<S, T>(source: Signal<S>, get: (s: S) => T, set: (s: S, value: T) => S): Signal<T> {
  return new Lens(source, get, set);
}

class Lens<S, T> extends Signal<T> {
  private readonly _source: Signal<S>;
  private readonly _read: ReadonlySignal<T>;
  private readonly _write: (s: S, value: T) => S;

  constructor(source: Signal<S>, get: (s: S) => T, set: (s: S, value: T) => S) {
    super(undefined as T);
    this._source = source;
    this._read = computed(() => get(source.value));
    this._write = set;
  }

  get value(): T {
    return this._read.value;
  }

  set value(x: T) {
    if (!Object.is(untracked(() => this._read.value), x))
      this._source.value = this._write(this._source.peek(), x);
  }
}

/** Whether `s` is a signal whose `value` can be written — the two-way gate. A prototype-descriptor
 * check: the vendored `Computed` defines a getter-only `value`, so computeds classify read-only
 * without any engine knowledge. */
export function isWritableSignal(s: unknown): s is Signal<unknown> {
  if (!(s instanceof Signal))
    return false;
  for (let proto = Object.getPrototypeOf(s); proto !== null; proto = Object.getPrototypeOf(proto)) {
    const descriptor = Object.getOwnPropertyDescriptor(proto, 'value');
    if (descriptor)
      return descriptor.set !== undefined;
  }
  return false;
}
