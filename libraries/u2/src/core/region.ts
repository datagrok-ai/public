/* Keyed and conditional regions: the host's children follow a list or a condition without a
   rebuild. `each` keeps one element per key, built once in its own scope and patched in place, so
   focus, caret and whatever the element owns survive a change elsewhere in the list; `show` swaps
   between two builds. Both own the host's children and are owned by the ambient scope. */
import {ReadonlySignal, Signal, untracked} from './signals.js';
import {Scope} from './scope.js';

export interface EachOptions<T> {
  /** Patches a kept key whose item changed (by identity); return false to rebuild it instead. */
  update?: (el: HTMLElement, item: T, index: number) => boolean | void;
  /** Shown while the list is empty. */
  empty?: () => HTMLElement;
}

interface Built<T> {
  el: HTMLElement;
  scope: Scope;
  item: T;
}

function owner(fn: string): Scope {
  const scope = Scope.ambient;
  if (!scope)
    throw new Error(`u2: ${fn}() needs an owner — wrap the code in Control.build(...) or component.runInScope(...)`);
  return scope;
}

function reader<T>(source: ReadonlySignal<T> | (() => T)): () => T {
  return source instanceof Signal ? () => (source as ReadonlySignal<T>).value : source as () => T;
}

/** Builds `make` in a scope of its own under `parent`. */
function built<T>(parent: Scope, label: string, item: T, make: () => HTMLElement): Built<T> {
  const scope = new Scope(parent, {label, lazy: true});
  return {el: Scope.runWith(scope, make), scope, item};
}

function drop(entry: Built<unknown>): void {
  entry.scope.dispose();
  entry.el.remove();
}

/** A keyed region: one element per key, built once in its own scope (a child of the ambient owner),
 * kept across changes (focus and caret survive), moved only when out of place, disposed when its key
 * leaves. Every element gets `data-u2-item = key` unless it carries one. A repeated key is warned
 * about once and suffixed `#<n>`. Throws without an owner. */
export function each<T>(host: HTMLElement, items: ReadonlySignal<readonly T[]> | (() => readonly T[]),
  keyOf: (item: T, index: number) => string, build: (item: T, index: number) => HTMLElement,
  options: EachOptions<T> = {}): void {
  const scope = owner('each');
  const read = reader(items);
  let entries = new Map<string, Built<T>>();
  let empty: Built<undefined> | undefined;
  const warned = new Set<string>();

  const reconcile = (list: readonly T[]): void => {
    const next = new Map<string, Built<T>>();
    for (let i = 0; i < list.length; i++) {
      const item = list[i];
      let key = keyOf(item, i);
      if (next.has(key)) {
        if (!warned.has(key)) {
          warned.add(key);
          console.warn(`u2: each(): duplicate key "${key}"`);
        }
        let n = 2;
        while (next.has(`${key}#${n}`))
          n++;
        key = `${key}#${n}`;
      }
      let entry = entries.get(key);
      if (entry !== undefined && entry.item !== item) {
        const kept = entry;
        const patched = options.update !== undefined &&
          Scope.runWith(kept.scope, () => options.update!(kept.el, item, i)) !== false;
        if (patched)
          kept.item = item;
        else {
          drop(kept);
          entry = undefined;
        }
      }
      if (entry === undefined) {
        entry = built(scope, 'item', item, () => build(item, i));
        if (entry.el.dataset.u2Item === undefined)
          entry.el.dataset.u2Item = key;
      }
      next.set(key, entry);
    }
    for (const [key, entry] of entries) {
      if (next.get(key) !== entry)
        drop(entry);
    }
    entries = next;
    if (list.length > 0 && empty !== undefined) {
      drop(empty);
      empty = undefined;
    }
    let at: ChildNode | null = host.firstChild;
    for (const {el} of next.values()) {
      if (el === at)
        at = at.nextSibling;
      else
        host.insertBefore(el, at);
    }
    if (list.length === 0 && empty === undefined && options.empty !== undefined) {
      empty = built(scope, 'empty', undefined, options.empty);
      host.append(empty.el);
    }
  };

  scope.effect(() => {
    const list = read();
    untracked(() => reconcile(list));
  });
}

/** A conditional region: `build` (own scope) while `when` is true, `otherwise` (own scope) while it
 * is false; each is built when its branch is entered and disposed when it is left. */
export function show(host: HTMLElement, when: ReadonlySignal<boolean> | (() => boolean),
  build: () => HTMLElement, otherwise?: () => HTMLElement): void {
  const scope = owner('show');
  const read = reader(when);
  let shown: boolean | undefined;
  let current: Built<undefined> | undefined;
  scope.effect(() => {
    const on = read() === true;
    if (on === shown)
      return;
    shown = on;
    untracked(() => {
      if (current !== undefined)
        current.scope.dispose();
      const make = on ? build : otherwise;
      current = make === undefined ? undefined : built(scope, 'show', undefined, make);
      host.replaceChildren(...current === undefined ? [] : [current.el]);
    });
  });
}
