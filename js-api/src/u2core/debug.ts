/* Leak and ownership tracking for u2. Off in production: only a live-scope counter is
   kept, never a reference. When on (`globalThis.__U2_DEBUG__`, `localStorage['u2.debug'] = '1'`,
   `U2_DEBUG=1` in Node, or `debug.enable()`), every scope, owned timer and global listener is
   recorded with its creation stack. One registry per page, shared by every engine copy. */

/** A live scope, as {@link U2Debug.liveScopes} reports it. */
export interface ScopeRecord {
  /** Unique per page, increasing in creation order. */
  id: number;
  /** What created it: a control kind (`text-input`), `effect`, `row`, `popup`, `dialog`, `scope`. */
  label: string;
  /** `data-u2-name` / tag of the control owning it, when known. */
  name?: string;
  /** `performance.now()` at creation. */
  createdAt: number;
  /** `new Error().stack` trimmed of engine frames; present when enabled with stacks. */
  stack?: string;
  /** The owning scope's id. */
  parent?: number;
  /** The control root is connected to the document (undefined: the scope has no root). */
  connected?: boolean;
}

/** A pending timeout, interval or animation frame created through a scope helper or while a u2
 * scope was current. */
export interface TimerRecord {
  /** Unique per page, increasing in creation order. */
  id: number;
  /** `setTimeout`, `setInterval` or `requestAnimationFrame`. */
  kind: 'timeout' | 'interval' | 'raf';
  /** The delay or period. */
  ms?: number;
  /** The owning scope's id. */
  owner?: number;
  /** `'singleton'` when the owner is kept for the page's life by design ({@link U2Debug.markSingleton}). */
  label?: string;
  /** `performance.now()` at creation. */
  createdAt: number;
  /** Creation stack; present when enabled with stacks. */
  stack?: string;
}

/** A listener on window, document or body added through a scope helper or while a u2 scope was
 * current. */
export interface ListenerRecord {
  /** Unique per page, increasing in creation order. */
  id: number;
  /** Which global target it listens on. */
  target: 'window' | 'document' | 'body';
  /** The event type. */
  type: string;
  /** Registered for the capture phase. */
  capture: boolean;
  /** The owning scope's id. */
  owner?: number;
  /** `'singleton'` when the owner is kept for the page's life by design ({@link U2Debug.markSingleton}). */
  label?: string;
  /** `performance.now()` at creation. */
  createdAt: number;
  /** Creation stack; present when enabled with stacks. */
  stack?: string;
}

/** What was live at one moment — the baseline {@link U2Debug.diff} compares against. */
export interface DebugSnapshot {
  /** `performance.now()` when taken. */
  at: number;
  /** Tracked live scopes. */
  scopes: number;
  /** Tracked pending timers. */
  timers: number;
  /** Tracked global listeners. */
  listeners: number;
  /** Ids live at snapshot time — `diff` reports what is live now and was not then. */
  scopeIds: number[];
  /** Timer ids live at snapshot time. */
  timerIds: number[];
  /** Listener ids live at snapshot time. */
  listenerIds: number[];
}

/** What a {@link U2Debug.diff} found live that was not live at the snapshot. */
export interface LeakReport {
  /** Scopes created since the snapshot and still live. */
  scopes: ScopeRecord[];
  /** Timers created since the snapshot and still pending. */
  timers: TimerRecord[];
  /** Global listeners added since the snapshot and still attached. */
  listeners: ListenerRecord[];
  /** scopes + timers + listeners, not counting those labelled `'singleton'`. 0 ⇔ back to baseline. */
  total: number;
}

/** The in-page leak and ownership inspector, `DG.U2.debug`. */
export interface U2Debug {
  /** True when tracking is on (dev/test). Production: false, and no global strong refs exist. */
  readonly enabled: boolean;
  /** Turns tracking on from now. `stacks` (default true) captures creation stacks.
   * `orphanSeconds` (default 10) is the orphan-sweep threshold; 0 disables the sweep. */
  enable(opts?: {stacks?: boolean, orphanSeconds?: number}): void;
  /** Turns tracking off and forgets every record. */
  disable(): void;
  /** Live scope count — always available (a counter, not a set), tracking on or off. */
  liveCount(): number;
  /** Scopes created since tracking was turned on and not disposed yet. */
  liveScopes(): ScopeRecord[];
  /** Pending timeouts/intervals/rAFs created through the scope helpers or while a u2 scope is current. */
  timers(): TimerRecord[];
  /** Listeners on window/document/body added through the scope helpers or while a u2 scope is current. */
  listeners(): ListenerRecord[];
  /** What is live now — pass it to {@link diff} later. */
  snapshot(): DebugSnapshot;
  /** What is live now and was not at `since`. */
  diff(since: DebugSnapshot): LeakReport;
  /** Live control scopes whose root was connected once and has been disconnected for more than
   * `orphanSeconds`. Each is also `console.warn`ed once. */
  orphans(): ScopeRecord[];
  /** Runs the orphan sweep now and answers {@link orphans}. */
  sweep(): ScopeRecord[];
  /** Marks a scope (or a component's scope) as kept for the page's life by design — an overlay host,
   * a tooltip service: it, its descendants and their timers and listeners report `label: 'singleton'`
   * and never count in {@link LeakReport.total} or as orphans. */
  markSingleton(scopeOrComponent: object): void;
}

interface ScopeEntry extends ScopeRecord {
  component?: WeakRef<object>;
  singleton?: boolean;
  everConnected?: boolean;
  lastConnected?: number;
  warned?: boolean;
}

interface ListenerEntry extends ListenerRecord {
  raw?: {target: EventTarget, handler: unknown};
}

type Fn = (...args: any[]) => any;

interface Originals {
  add?: Fn;
  remove?: Fn;
  setInterval: Fn;
  clearInterval: Fn;
  raf?: Fn;
  caf?: Fn;
}

interface Registry {
  enabled: boolean;
  stacks: boolean;
  orphanSeconds: number;
  nextId: number;
  live: number;
  ids: WeakMap<object, number>;
  scopes: Map<number, ScopeEntry>;
  timers: Map<number, TimerRecord>;
  listeners: Map<number, ListenerEntry>;
  rawTimers: Map<unknown, number>;
  /** > 0 while a scope helper registers — the raw patches stay out of its way. */
  helper: number;
  /** One getter per engine copy: the scope current in that copy, if any. */
  ambients: (() => object | undefined)[];
  originals?: Originals;
  sweeper?: unknown;
}

const REGISTRY = Symbol.for('u2.debug.registry');
const g = globalThis as Record<PropertyKey, any>;

/** @internal The page-wide registry; `Scope` keeps its live counter here. */
export const registry: Registry = g[REGISTRY] ??= {enabled: false, stacks: true, orphanSeconds: 10, nextId: 1,
  live: 0, ids: new WeakMap(), scopes: new Map(), timers: new Map(), listeners: new Map(), rawTimers: new Map(),
  helper: 0, ambients: []};

const ENGINE_FRAME = /u2core[\\/]|\bScope[.\s]|\bnew Scope\b|\bget scope\b|\b(captureStack|trackScope|trackTimer|trackListener|describeScope|recordRaw\w*)\b/;

/** @internal Creation stack minus the engine's own frames; undefined when stacks are off. */
export function captureStack(): string | undefined {
  if (!registry.stacks)
    return undefined;
  const E = Error as {stackTraceLimit?: number};
  const limit = E.stackTraceLimit;
  E.stackTraceLimit = 40;
  const lines = (new Error().stack ?? '').split('\n');
  E.stackTraceLimit = limit;
  let i = 1;
  while (i < lines.length && ENGINE_FRAME.test(lines[i]))
    i++;
  return lines.slice(i).join('\n');
}

function now(): number {
  return performance.now();
}

function currentOwner(): {current: boolean, id?: number} {
  for (const ambient of registry.ambients) {
    const scope = ambient();
    if (scope !== undefined)
      return {current: true, id: registry.ids.get(scope)};
  }
  return {current: false};
}

function globalTarget(target: unknown): ListenerRecord['target'] | undefined {
  const doc = g.document;
  return target === g.window ? 'window' : doc !== undefined && target === doc ? 'document' :
    doc?.body !== undefined && target === doc.body ? 'body' : undefined;
}

function isCapture(options: unknown): boolean {
  return options === true || (typeof options === 'object' && options !== null &&
    (options as AddEventListenerOptions).capture === true);
}

/** @internal Records a scope (when tracking is on); `parent` links it to its owner. */
export function trackScope(scope: object, label: string, parent: object | undefined): void {
  if (!registry.enabled)
    return;
  const id = registry.nextId++;
  registry.ids.set(scope, id);
  registry.scopes.set(id, {id, label, createdAt: now(), stack: captureStack(),
    parent: parent === undefined ? undefined : registry.ids.get(parent)});
}

/** @internal */
export function forgetScope(scope: object): void {
  const id = registry.ids.get(scope);
  if (id !== undefined)
    registry.scopes.delete(id);
}

/** @internal Attaches the component a scope belongs to: its label, name and root come from it. */
export function describeScope(scope: object, component: object): void {
  if (!registry.enabled || typeof WeakRef === 'undefined')
    return;
  const id = registry.ids.get(scope);
  const entry = id === undefined ? undefined : registry.scopes.get(id);
  if (entry !== undefined)
    entry.component = new WeakRef(component);
}

/** @internal Records an owned timer; answers its id for {@link forgetTimer}. */
export function trackTimer(kind: TimerRecord['kind'], ms: number | undefined, owner: object): number | undefined {
  if (!registry.enabled)
    return undefined;
  const id = registry.nextId++;
  registry.timers.set(id, {id, kind, ms, owner: registry.ids.get(owner), createdAt: now(), stack: captureStack()});
  return id;
}

/** @internal */
export function forgetTimer(id: number | undefined): void {
  if (id !== undefined)
    registry.timers.delete(id);
}

/** @internal Records an owned listener when the target is window, document or body. */
export function trackListener(target: EventTarget, type: string, options: unknown, owner: object): number | undefined {
  const kind = registry.enabled ? globalTarget(target) : undefined;
  if (kind === undefined)
    return undefined;
  const id = registry.nextId++;
  registry.listeners.set(id, {id, target: kind, type, capture: isCapture(options), owner: registry.ids.get(owner),
    createdAt: now(), stack: captureStack()});
  return id;
}

/** @internal */
export function forgetListener(id: number | undefined): void {
  if (id !== undefined)
    registry.listeners.delete(id);
}

/** @internal Runs a helper's own registration with the raw patches muted. */
export function asHelper<T>(fn: () => T): T {
  registry.helper++;
  try {
    return fn();
  }
  finally {
    registry.helper--;
  }
}

function recordRawListener(target: EventTarget, type: string, handler: unknown, options: unknown): void {
  try {
    if (registry.helper > 0 || handler == null)
      return;
    const kind = globalTarget(target);
    if (kind === undefined)
      return;
    // the browser removes these on its own, with no removeEventListener to forget them by
    const opts = typeof options === 'object' && options !== null ? options as AddEventListenerOptions : undefined;
    if (opts?.once === true || opts?.signal?.aborted === true)
      return;
    const owner = currentOwner();
    if (!owner.current)
      return;
    const id = registry.nextId++;
    registry.listeners.set(id, {id, target: kind, type, capture: isCapture(options), owner: owner.id,
      createdAt: now(), stack: captureStack(), raw: {target, handler}});
    opts?.signal?.addEventListener('abort', () => registry.listeners.delete(id), {once: true});
  }
  catch (_) {
  }
}

function forgetRawListener(target: EventTarget, type: string, handler: unknown, options: unknown): void {
  try {
    if (registry.listeners.size === 0 || globalTarget(target) === undefined)
      return;
    const capture = isCapture(options);
    for (const [id, l] of registry.listeners) {
      if (l.raw?.target === target && l.raw.handler === handler && l.type === type && l.capture === capture) {
        registry.listeners.delete(id);
        return;
      }
    }
  }
  catch (_) {
  }
}

function recordRawTimer(kind: TimerRecord['kind'], handle: unknown, ms?: number): number | undefined {
  try {
    if (registry.helper > 0)
      return undefined;
    const owner = currentOwner();
    if (!owner.current)
      return undefined;
    const id = registry.nextId++;
    registry.timers.set(id, {id, kind, ms, owner: owner.id, createdAt: now(), stack: captureStack()});
    registry.rawTimers.set(handle, id);
    return id;
  }
  catch (_) {
    return undefined;
  }
}

function forgetRawTimer(handle: unknown): void {
  const id = registry.rawTimers.get(handle);
  if (id === undefined)
    return;
  registry.rawTimers.delete(handle);
  registry.timers.delete(id);
}

function patch(): void {
  if (registry.originals !== undefined)
    return;
  const proto = g.EventTarget?.prototype;
  const o: Originals = registry.originals = {add: proto?.addEventListener, remove: proto?.removeEventListener,
    setInterval: g.setInterval, clearInterval: g.clearInterval, raf: g.requestAnimationFrame,
    caf: g.cancelAnimationFrame};
  if (proto !== undefined && o.add !== undefined && o.remove !== undefined) {
    proto.addEventListener = function(this: EventTarget, type: string, handler: unknown, options?: unknown) {
      recordRawListener(this, type, handler, options);
      return o.add!.call(this, type, handler, options);
    };
    proto.removeEventListener = function(this: EventTarget, type: string, handler: unknown, options?: unknown) {
      forgetRawListener(this, type, handler, options);
      return o.remove!.call(this, type, handler, options);
    };
  }
  g.setInterval = function(handler: unknown, ms?: number, ...args: unknown[]) {
    const handle = o.setInterval.call(g, handler, ms, ...args);
    recordRawTimer('interval', handle, ms);
    return handle;
  };
  g.clearInterval = function(handle: unknown) {
    forgetRawTimer(handle);
    return o.clearInterval.call(g, handle);
  };
  if (o.raf !== undefined && o.caf !== undefined) {
    g.requestAnimationFrame = function(fn: (time: number) => void) {
      let handle: unknown;
      handle = o.raf!.call(g, (time: number) => {
        forgetRawTimer(handle);
        fn(time);
      });
      recordRawTimer('raf', handle);
      return handle;
    };
    g.cancelAnimationFrame = function(handle: unknown) {
      forgetRawTimer(handle);
      return o.caf!.call(g, handle);
    };
  }
}

function unpatch(): void {
  const o = registry.originals;
  if (o === undefined)
    return;
  registry.originals = undefined;
  const proto = g.EventTarget?.prototype;
  if (proto !== undefined && o.add !== undefined && o.remove !== undefined) {
    proto.addEventListener = o.add;
    proto.removeEventListener = o.remove;
  }
  g.setInterval = o.setInterval;
  g.clearInterval = o.clearInterval;
  if (o.raf !== undefined && o.caf !== undefined) {
    g.requestAnimationFrame = o.raf;
    g.cancelAnimationFrame = o.caf;
  }
}

interface Described {
  componentMeta?: {tag?: string};
  root?: {dataset?: Record<string, string | undefined>, isConnected?: boolean};
  name?: string;
}

function rootOf(entry: ScopeEntry): Described['root'] | undefined {
  try {
    const root = (entry.component?.deref() as Described | undefined)?.root;
    return root != null && typeof root === 'object' && 'isConnected' in root ? root : undefined;
  }
  catch (_) {
    return undefined;
  }
}

function isSingleton(id: number | undefined): boolean {
  for (let entry = id === undefined ? undefined : registry.scopes.get(id); entry !== undefined;
    entry = entry.parent === undefined ? undefined : registry.scopes.get(entry.parent)) {
    if (entry.singleton)
      return true;
  }
  return false;
}

function publicRecord(entry: ScopeEntry): ScopeRecord {
  const record: ScopeRecord = {id: entry.id, label: entry.label, createdAt: entry.createdAt};
  try {
    const c = entry.component?.deref() as Described | undefined;
    if (isSingleton(entry.id))
      record.label = 'singleton';
    else if (c !== undefined) {
      const root = rootOf(entry);
      record.label = c.componentMeta?.tag ?? root?.dataset?.u2 ?? c.constructor.name;
      record.name = root?.dataset?.u2Name ?? c.name;
      if (root !== undefined)
        record.connected = root.isConnected === true;
    }
  }
  catch (_) {
  }
  if (entry.name !== undefined)
    record.name = entry.name;
  if (entry.stack !== undefined)
    record.stack = entry.stack;
  if (entry.parent !== undefined)
    record.parent = entry.parent;
  return record;
}

function owned<T extends TimerRecord | ListenerRecord>(record: T): T {
  const {raw: _, ...rest} = record as T & {raw?: unknown};
  return (isSingleton(record.owner) ? {...rest, label: 'singleton'} : rest) as T;
}

function counted(records: {label?: string}[]): number {
  return records.filter((r) => r.label !== 'singleton').length;
}

function isOrphan(entry: ScopeEntry, at: number): boolean {
  return entry.everConnected === true && rootOf(entry)?.isConnected !== true && !isSingleton(entry.id) &&
    at - (entry.lastConnected ?? at) > registry.orphanSeconds * 1000;
}

function sweep(): ScopeRecord[] {
  const at = now();
  const found: ScopeRecord[] = [];
  for (const entry of registry.scopes.values()) {
    const root = rootOf(entry);
    if (root === undefined)
      continue;
    if (root.isConnected) {
      entry.everConnected = true;
      entry.lastConnected = at;
    }
    else if (registry.orphanSeconds > 0 && isOrphan(entry, at)) {
      const record = publicRecord(entry);
      found.push(record);
      if (!entry.warned) {
        entry.warned = true;
        console.warn(`u2: orphan ${record.label} (${record.name ?? ''}) disconnected for ` +
          `${Math.round((at - entry.lastConnected!) / 1000)}s`, entry.stack ?? '');
      }
    }
  }
  return found;
}

function startSweeper(): void {
  const o = registry.originals;
  if (registry.sweeper !== undefined || registry.orphanSeconds <= 0 || o === undefined)
    return;
  const timer = o.setInterval.call(g, sweep, 5000);
  (timer as {unref?: () => void})?.unref?.();
  registry.sweeper = timer;
}

function stopSweeper(): void {
  if (registry.sweeper === undefined)
    return;
  (registry.originals?.clearInterval ?? g.clearInterval).call(g, registry.sweeper);
  registry.sweeper = undefined;
}

/** `DG.U2.debug` — live scopes with creation stacks, owned timers and global listeners, the
 * snapshot/diff the BDD leak check runs, and the orphan sweep. */
export const debug: U2Debug = {
  get enabled(): boolean {
    return registry.enabled;
  },

  enable(opts?: {stacks?: boolean, orphanSeconds?: number}): void {
    registry.stacks = opts?.stacks ?? true;
    registry.orphanSeconds = opts?.orphanSeconds ?? 10;
    registry.enabled = true;
    patch();
    stopSweeper();
    startSweeper();
  },

  disable(): void {
    stopSweeper();
    unpatch();
    registry.enabled = false;
    registry.scopes.clear();
    registry.timers.clear();
    registry.listeners.clear();
    registry.rawTimers.clear();
  },

  liveCount: () => registry.live,
  liveScopes: () => [...registry.scopes.values()].map(publicRecord),
  timers: () => [...registry.timers.values()].map(owned),
  listeners: () => [...registry.listeners.values()].map(owned),

  snapshot(): DebugSnapshot {
    return {at: now(), scopes: registry.scopes.size, timers: registry.timers.size, listeners: registry.listeners.size,
      scopeIds: [...registry.scopes.keys()], timerIds: [...registry.timers.keys()],
      listenerIds: [...registry.listeners.keys()]};
  },

  diff(since: DebugSnapshot): LeakReport {
    const scopeIds = new Set(since.scopeIds);
    const timerIds = new Set(since.timerIds);
    const listenerIds = new Set(since.listenerIds);
    const scopes = [...registry.scopes.values()].filter((s) => !scopeIds.has(s.id)).map(publicRecord);
    const timers = [...registry.timers.values()].filter((t) => !timerIds.has(t.id)).map(owned);
    const listeners = [...registry.listeners.values()].filter((l) => !listenerIds.has(l.id)).map(owned);
    return {scopes, timers, listeners, total: counted(scopes) + counted(timers) + counted(listeners)};
  },

  orphans(): ScopeRecord[] {
    const at = now();
    return [...registry.scopes.values()].filter((e) => isOrphan(e, at)).map(publicRecord);
  },

  sweep,

  markSingleton(scopeOrComponent: object): void {
    const scope = registry.ids.has(scopeOrComponent) ? scopeOrComponent : (scopeOrComponent as {scope?: object}).scope;
    const id = scope === undefined ? undefined : registry.ids.get(scope);
    const entry = id === undefined ? undefined : registry.scopes.get(id);
    if (entry !== undefined)
      entry.singleton = true;
  },
};

/** `__U2_DEBUG__` is `true` or the {@link U2Debug.enable} options. */
function environmentOptions(): {stacks?: boolean, orphanSeconds?: number} | undefined {
  const flag = g.__U2_DEBUG__;
  if (flag === true || (typeof flag === 'object' && flag !== null))
    return flag === true ? {} : flag;
  try {
    if (g.localStorage?.getItem('u2.debug') === '1')
      return {};
  }
  catch (_) {
  }
  return g.process?.env?.U2_DEBUG === '1' ? {} : undefined;
}

const fromEnvironment = registry.enabled ? undefined : environmentOptions();
if (fromEnvironment !== undefined)
  debug.enable(fromEnvironment);
