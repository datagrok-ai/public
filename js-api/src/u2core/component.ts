import {Scope, whenActive} from './scope.js';
import type {TaskContext, TaskOptions, TaskRun} from './scope.js';
import {Signal, batch, computed, isWritableSignal, signal} from './signals.js';
import {describeScope, registry} from './debug.js';
import {Usage} from './usage.js';
import {propertyFields} from './property-fields.js';
import type {IProperty} from './property-fields.js';
import type {BindProp, BindSource, ComponentMetaBase, FuncLike, NamedProperty, ObservableLike,
  IWidgetStatus, IControlState, IItemsState, JsonValue} from './protocol.js';

/** What a widget announces when one of its properties changed: the property's name, or null when
 * several changed at once (or the widget cannot say which). Structural — nothing platform-typed. */
export type PropertyChange = string | null;

/** One property-tier step: `inner` is what the step reads, `known` the value the tier last
 * observed — read in a refresh or written by the effect — so the write effect never compares
 * against a fresh read (a read answering a new wrapper per call would spin otherwise). */
interface PropertyEntry {
  inner: Signal<unknown>;
  step: Signal<unknown>;
  known: unknown;
}

/** Labels of the scopes torn down and built again while their owner lives: rows of `each()` and
 * `show()`, and effect runs. */
const REBUILT: Record<string, string> = {item: 'each() row', row: 'row', show: 'show() branch',
  empty: 'each() empty state', effect: 'effect run'};
/** The component each component scope belongs to — the adoption link a scope's `parent` lacks. */
const scopeOwners = new WeakMap<Scope, Component>();
/** Per registering scope: the rebuilt scope each call site first registered from. */
const rebuiltSites = new WeakMap<Scope, Map<string, WeakRef<Scope>>>();
const warnedSites = new Set<string>();

/** The per-instance u2 state, kept under one own field so that {@link Component.init} can be
 * applied to a platform `DG.Widget` (a Component by inheritance) without clobbering its own fields. */
export interface ComponentState {
  scope?: Scope;
  aiDescription: string | null;
  properties: NamedProperty[] | null;
  propertiesMeta?: ComponentMetaBase;
  functions: FuncLike[];
  listeners: {id: string | null, next: (x: unknown) => void}[];
  propertySignals: Map<string, PropertyEntry>;
  propertyIndex?: Map<string, NamedProperty>;
  warned?: Set<string>;
  changesWired: boolean;
  /** Names a notification marked since the last refresh; 'all' for an anonymous change. */
  stale?: Set<string> | 'all';
  table?: unknown;
  detaching?: boolean;
  platformKilled?: boolean;
  propertyChangeListeners?: Set<(c: PropertyChange) => void>;
  name?: string;
  /** The ambient scope that adopted this component and the release it registered there — dropped
   * on {@link Component.dispose} so a long-lived owner does not pile up dead closures. */
  adopter?: {scope: Scope, release: () => void};
  /** Disposed because its owner went, not by hand. */
  byOwner?: boolean;
  /** With u2 debug on: the components whose constructors were running when this one was built, by class name. */
  builtDuring?: string[];
}

/** The members a field of a subclass must not take: overriding one breaks the identity, the lifecycle or the
 * automation u2 runs on every control. An accessor is taken by any own value, a method by a non-function. */
const RESERVED_ACCESSORS = ['name', 'root', 'scope', 'visible', 'rootId'];
const RESERVED_METHODS = ['effect', 'own', 'on', 'dispose', 'runInScope', 'build', 'task', 'subscribe', 'getState',
  'getItems', 'onEvent', 'fireEvent'];
/** With u2 debug on: the classes of the components whose field check is pending (their constructors may still be
 * running), with how many of each. */
const constructing = new Map<string, number>();

/** Base of every u2 component, visual or not: a name, an effect/cleanup scope, and the widget
 * introspection surface (DD7) generated from {@link componentMeta} — what the shell, the context
 * panel, automation and the copilot ask of a widget, implemented by `DG.Widget` through
 * inheritance. A component constructed inside
 * another's scope is disposed with it. Platform-free; in Datagrok, `DG.Widget extends Control`,
 * so the kill channel and the introspection surface are this class's own — nothing wires or
 * delegates. Standalone hosts (gallery, tests) call {@link dispose} directly. */
export class Component implements BindSource {
  /** A root has no owner by design (a dialog): built with none, it is not reported as `owner-missing`. */
  static isRoot = false;

  /** Minted on first use; assigned before {@link _scopeMinted} fires, so the hook's own
   * `this.scope.own(...)` reads the stored scope instead of re-entering the mint. Lazy: counted and tracked only
   * once it owns something, so a static control costs nothing. */
  get scope(): Scope {
    const u2 = this._u2;
    if (u2.scope === undefined) {
      u2.scope = new Scope(undefined, {lazy: true});
      scopeOwners.set(u2.scope, this);
      describeScope(u2.scope, this);
      if (u2.adopter === undefined && this.ambientAdoption && !(this.constructor as typeof Component).isRoot)
        whenActive(u2.scope, () => Component._reportUnowned(this));
      this._scopeMinted(u2.scope);
    }
    return u2.scope;
  }

  /** First-mint hook: a platform descendant wires its lifecycle here. Never fired for a
   * scope provided to {@link init} — that path runs inside the constructor chain. */
  protected _scopeMinted(_scope: Scope): void {}
  /** Unique within a spec: the anchor for selection, patches and automation. An accessor pair
   * (not a field), so a platform descendant's own `name` accessors override it legally (TS2611). */
  get name(): string | undefined { return this._u2.name; }
  set name(x: string | undefined) { this._u2.name = x; }

  /** Sets the introspection/form key without claiming an automation id — the label-fallback path:
   * a caption is not a stable identity, so it never reaches `data-u2-name`. */
  protected deriveName(x: string | undefined): void {
    this._u2.name = x;
  }
  /** The registry metadata of the tag that built this component, stamped by the spec renderer —
   * the single source the introspection surface below is generated from. */
  componentMeta?: ComponentMetaBase;
  /** The props the registry's `create` was handed, stamped alongside {@link componentMeta}: the last resort
   * of the property read below, for a component that keeps a prop under another name or not at all
   * (a `create` wrapping a bare element). Construction-time values — a live one always wins. */
  specProps?: Record<string, unknown>;
  readonly _u2!: ComponentState;

  constructor(scope?: Scope) {
    Component.init(this, scope);
  }

  /** Per-instance state, for a constructed component and for an adopted platform object alike. */
  static init(self: Component, scope?: Scope): void {
    const it = self as {_u2: ComponentState};
    it._u2 = {aiDescription: null, properties: null, functions: [], listeners: [],
      propertySignals: new Map(), changesWired: false};
    if (scope !== undefined) {
      it._u2.scope = scope;
      scopeOwners.set(scope, self);
      describeScope(scope, self);
    }
    const owner = self.ambientAdoption ? Scope.ambient : undefined;
    if (owner && owner !== it._u2.scope) {
      if (owner.isDisposed) {
        Usage.report('owner-disposed', `${self.constructor.name} was built under a scope already disposed ` +
          `(${owner.label}), so it was disposed at once — a handler that closed its own owner and then built? Build ` +
          'it under a live owner; one that must outlive it under `Scope.root(\'<what>\')`.');
      }
      const release = () => {
        it._u2.byOwner = true;
        self.dispose();
      };
      it._u2.adopter = {scope: owner, release};
      owner.own(release);
    }
    if (registry.enabled)
      Component._watchFields(self);
    else if (Usage.strict)
      queueMicrotask(() => Component._checkReserved(self));
  }

  private static _reportUnowned(self: Component): void {
    const type = self.constructor.name || 'A control';
    Usage.report('owner-missing', `${type} was built with no owner and now owns something nobody will dispose. ` +
      `Build it under an owner (\`owner.build(() => new ${type}(…))\`, \`component((self) => …)\`); after an ` +
      '`await`, `owner.task(async (t) => … t.build(…))`; a new root on purpose, under `Scope.root(\'<what>\')`.');
  }

  /** With u2 debug on: notes which components' constructors are running (their `new X` frames), and checks this
   * one's fields a microtask later, once its constructor has run. A best-effort check: the `new X (` frame is V8's
   * stack format, a minified build's class names can match the wrong frame or none, and each child built during
   * another's constructor captures a 40-frame stack. */
  private static _watchFields(self: Component): void {
    const type = self.constructor.name;
    if (constructing.size > 0) {
      const E = Error as {stackTraceLimit?: number};
      const limit = E.stackTraceLimit;
      E.stackTraceLimit = 40;
      const stack = new Error().stack ?? '';
      E.stackTraceLimit = limit;
      const during = [...constructing.keys()].filter((name) => stack.includes(`new ${name} (`));
      if (during.length > 0)
        self._u2.builtDuring = during;
    }
    constructing.set(type, (constructing.get(type) ?? 0) + 1);
    queueMicrotask(() => {
      const left = constructing.get(type)! - 1;
      if (left === 0)
        constructing.delete(type);
      else
        constructing.set(type, left);
      Component._checkFields(self);
    });
  }

  private static _checkReserved(self: Component): Set<string> {
    const taken = new Set<string>();
    const fields = self as unknown as Record<string, unknown>;
    for (const key of Object.keys(self)) {
      const method = RESERVED_METHODS.includes(key);
      if (RESERVED_ACCESSORS.includes(key) || method && typeof fields[key] !== 'function') {
        taken.add(key);
        Usage.report('reserved-key', `${self.constructor.name}.${key} overrides Control's own ` +
          `\`${key}${method ? '()' : ''}\`, which u2 runs on every control — name the field for what it holds ` +
          `(\`${key}Input\`, \`${key}Text\`).`);
      }
    }
    return taken;
  }

  /** A field that takes one of the control's own members (`reserved-key`), and a child control its constructor
   * built outside its scope — owned by whoever was building this one, so it outlives this one (a warning). */
  private static _checkFields(self: Component): void {
    const type = self.constructor.name;
    const taken = Component._checkReserved(self);
    for (const key of Object.keys(self)) {
      const value = (self as unknown as Record<string, unknown>)[key];
      if (taken.has(key) || value === self || !(value instanceof Component) || !value._u2.builtDuring?.includes(type))
        continue;
      const adopter = value._u2.adopter?.scope;
      if (adopter !== undefined && adopter === self._u2.adopter?.scope) {
        console.warn(`u2: ${type}.${key} holds a control built outside ${type}'s scope: it belongs to whoever built ` +
          `${type} and outlives it unless disposed by hand. Build it under \`this.build(() => new …)\`.`);
      }
    }
  }

  /** Structural probe by methods only — never `.scope`: the accessor mints. */
  static is(x: unknown): x is Component {
    if (x instanceof Component)
      return true;
    const c = x as Partial<Component> | null | undefined;
    return typeof c?.bindStep === 'function' && typeof c?.bindProps === 'function' &&
      typeof c?.own === 'function' && typeof c?.dispose === 'function';
  }

  /** The equality the property tier and {@link link} suppress echoes by: `Object.is`, arrays
   * element-wise, plain objects key-wise (recursive), everything else by identity. */
  static sameValue(a: unknown, b: unknown): boolean {
    if (Object.is(a, b))
      return true;
    if (Array.isArray(a) && Array.isArray(b))
      return a.length === b.length && a.every((x, i) => Component.sameValue(x, b[i]));
    if (!Component._isPlain(a) || !Component._isPlain(b))
      return false;
    const keys = Object.keys(a);
    return keys.length === Object.keys(b).length &&
      keys.every((k) => k in b && Component.sameValue(a[k], b[k]));
  }

  protected static _isPlain(x: unknown): x is Record<string, unknown> {
    if (x === null || typeof x !== 'object')
      return false;
    const proto = Object.getPrototypeOf(x);
    return proto === Object.prototype || proto === null;
  }

  /** Runs `fn` with this component's scope ambient: what it builds is owned by this component.
   * @deprecated Use the owner verb `build(fn)` ({@link OwnerVerbs.build}): `this.build(fn)`, `dialog.build(fn)`. */
  runInScope<T>(fn: () => T): T {
    return this.scope.build(fn);
  }

  /** Runs `fn` now and on every change of the signals it reads, until this component is disposed.
   * Called while a row of `each()` or an effect run is being built, it belongs on a control built
   * there (`input.effect(...)`) — on a longer-lived one it piles up with every rebuild, and u2 warns. */
  effect(fn: () => void): this {
    Component._checkRebuilt(this, 'effect');
    this.scope.effect(fn);
    return this;
  }

  /** Registers a cleanup that runs when this component is disposed. */
  own(dispose: () => void): this {
    Component._checkRebuilt(this, 'own');
    this.scope.own(dispose);
    return this;
  }

  /** Warns when one call site registers on `self` from two different rebuilt scopes (rows of an
   * `each()`, runs of an effect) that `self` does not live inside: what it registers outlives each of
   * them and is added again by every rebuild. One warning per call site; with u2 debug tracking on (dev,
   * tests, the bdd leak audit), as it captures a stack. */
  private static _checkRebuilt(self: Component, kind: 'effect' | 'own'): void {
    if (!registry.enabled)
      return;
    const ambient = Scope.ambient;
    const target = self._u2.scope;
    if (ambient === undefined || target === undefined || ambient === target)
      return;
    let rebuilt: Scope | undefined;
    for (let s: Scope | undefined = ambient, hops = 0; s !== undefined && hops < 64; s = Component._ownerOf(s), hops++) {
      if (s === target)
        return;
      if (REBUILT[s.label] !== undefined) {
        rebuilt = s;
        break;
      }
    }
    if (rebuilt === undefined)
      return;
    for (let s: Scope | undefined = target, hops = 0; s !== undefined && hops < 64; s = Component._ownerOf(s), hops++) {
      if (s === rebuilt)
        return;
    }
    const site = (new Error().stack ?? '').split('\n').slice(3, 6).join('\n');
    let sites = rebuiltSites.get(target);
    if (sites === undefined)
      rebuiltSites.set(target, sites = new Map());
    const first = sites.get(site);
    if (first === undefined) {
      sites.set(site, new WeakRef(rebuilt));
      return;
    }
    if (first.deref() === rebuilt || warnedSites.has(site))
      return;
    warnedSites.add(site);
    const root = (self as {root?: HTMLElement}).root;
    const who = root?.dataset.u2Name ?? root?.dataset.u2 ?? self.constructor.name;
    const what = REBUILT[rebuilt.label];
    console.warn(`u2: ${kind}() on ${who} is called from inside a rebuilt ${what} and lives until ${who} is ` +
      `disposed — one more per rebuild. Register it on a control built in the ${what} (\`input.${kind}(...)\`) ` +
      `or on its scope (\`Scope.ambient.${kind}(...)\`).\n${site}`);
  }

  /** The scope `scope` is disposed with: its parent, or the scope that adopted its component. */
  private static _ownerOf(scope: Scope): Scope | undefined {
    return scope.parent ?? scopeOwners.get(scope)?._u2.adopter?.scope;
  }

  /** On a never-engaged wrapper the getter mints here and {@link _scopeMinted} wires the
   * lifecycle (`own(kill)` + the platform cleanup for Dart-owned), so dispose still kills it;
   * re-entry from the kill cleanup is a no-op via the scope's disposed flag. */
  dispose(): void {
    const adopter = this._u2.adopter;
    if (adopter !== undefined) {
      this._u2.adopter = undefined;
      adopter.scope.disown(adopter.release);
    }
    this.scope.dispose();
  }

  /** A short AI-facing briefing, seeded from the registry's usage or description. */
  get aiDescription(): string | null {
    return this._u2.aiDescription ?? this.componentMeta?.usage ?? this.componentMeta?.description ?? null;
  }

  set aiDescription(x: string | null) {
    this._u2.aiDescription = x;
  }

  /** The props {@link componentMeta} declares, each closed over the component's own state: a same-named
   * signal reads and writes live, a same-named plain member, an `options` entry or — last —
   * a {@link specProps} entry reads read-only. {@link name} — the component's spec identity —
   * is never a property; an input's
   * declared `name` prop is its form key, a different thing that IS a property like any other.
   * A component never hand-writes its property list. The list is generated once per {@link componentMeta}:
   * the accessors read live state, so only a re-stamp invalidates it. */
  getProperties(): NamedProperty[] {
    const u2 = this._u2;
    if (u2.properties && u2.propertiesMeta === this.componentMeta)
      return u2.properties;
    const self = this as unknown as Record<string, unknown>;
    u2.propertiesMeta = this.componentMeta;
    u2.propertyIndex = undefined;
    u2.properties = (this.componentMeta?.props ?? []).map((p) => {
      const name = p.name;
      const prop: NamedProperty = {...p, get: () => Component._read(self, name)};
      if (self[name] instanceof Signal) {
        prop.set = (_source: unknown, value: unknown) => {
          (self[name] as Signal<unknown>).value = value;
        };
      } else if (Component._settable(self, name)) {
        prop.set = (_source: unknown, value: unknown) => {
          self[name] = value;
        };
      }
      return prop;
    });
    return u2.properties;
  }

  /** A declared prop backed by an accessor with a setter is writable through it — the setter
   * carries the live behavior. A plain data field stays read-only: a bare assignment would change
   * nothing the component watches. */
  private static _settable(self: object, name: string): boolean {
    for (let at: object | null = self; at !== null; at = Object.getPrototypeOf(at)) {
      const descriptor = Object.getOwnPropertyDescriptor(at, name);
      if (descriptor !== undefined)
        return descriptor.set !== undefined;
    }
    return false;
  }

  getFunctions(): FuncLike[] {
    return this._u2.functions;
  }

  /** Events fired by this component, and — on a {@link Control} — the DOM events of that name
   * raised on its root. Without an id, every component event. The two streams are merged as they
   * are, so a subscriber sees either payload and must discriminate (see {@link ObservableLike}). */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    return {
      subscribe: (next: (x: unknown) => void) => {
        const listeners = this._u2.listeners;
        const listener = {id: eventId, next};
        listeners.push(listener);
        const off = () => {
          const i = listeners.indexOf(listener);
          if (i >= 0)
            listeners.splice(i, 1);
        };
        this.scope.own(off);
        return {unsubscribe: () => {
          off();
          this.scope.disown(off);
        }};
      },
    };
  }

  getWidgetStatus(): IWidgetStatus {
    return {
      parts: {},
      hitAreas: {},
      shortcuts: {},
      events: (this.componentMeta?.events ?? []).map((name) => ({name, eventName: name, description: ''})),
      description: this.componentMeta?.description ?? null,
      error: null,
    };
  }

  /** The second binding tier is on: every property {@link getProperties} declares is a bind step,
   * read and written on {@link propertyTarget}. Off for a plain u2 component, whose signal-less
   * meta props are construction-time options. */
  get propertyTier(): boolean {
    return false;
  }

  /** What this component's properties are read from and written to — `prop.get(target)` /
   * `prop.set(target, v)`: the component itself, unless an adoption says otherwise. */
  get propertyTarget(): unknown {
    return this;
  }

  /** Whether a component constructed inside an ambient scope is disposed with it. */
  protected get ambientAdoption(): boolean {
    return true;
  }

  /** What refreshes the property tier's cached signals from outside: null where nothing announces
   * a change (one-way). */
  protected propertyChanges(): ObservableLike<PropertyChange> | null {
    return null;
  }

  /** Keeps this component's `name` step and `source` equal, echo-suppressed by value; the reverse
   * direction only where `twoWay`. A component that took `source` as its own signal needs no link
   * and gets none; a read-only step cannot follow anything and gets none either. */
  link(name: string, source: Signal<unknown>, twoWay: boolean): void {
    const own = this.bindStep(name);
    if (!(own instanceof Signal) || own === source)
      return;
    if (!isWritableSignal(own)) {
      const owner = this.name ?? this.componentMeta?.tag ?? 'component';
      this._warnOnce(name, `${owner}.${name} is read-only — edits will not flow back`);
      return;
    }
    const scope = Scope.ambient ?? this.scope;
    scope.effect(() => {
      const value = source.value;
      if (!Component.sameValue(own.peek(), value))
        own.value = value;
    });
    if (!twoWay)
      return;
    scope.effect(() => {
      const value = own.value;
      if (!Component.sameValue(source.peek(), value))
        source.value = value;
    });
  }

  /** One binding step (DD5): a meta-declared prop's same-named Signal member — a named node is a
   * bind source, generated off {@link componentMeta} the way {@link getProperties} is — then, under
   * {@link propertyTier}, a cached signal over the declared property. '' answers the default
   * binding: the component's own `value` signal, where one exists. */
  bindStep(name: string): Signal<unknown> | BindSource | null {
    const self = this as unknown as Record<string, unknown>;
    const member = name === '' ? self.value :
      this.componentMeta?.props?.some((p) => p.name === name) ? self[name] : undefined;
    if (member instanceof Signal)
      return member as Signal<unknown>;
    return this.propertyTier && name !== '' ? this._propertySignal(name) : null;
  }

  /** The signal-backed subset of the meta props — what a binding picker offers on a named node —
   * plus, under {@link propertyTier}, every declared property. Allocates no signal. */
  bindProps(): BindProp[] {
    const self = this as unknown as Record<string, unknown>;
    const props: BindProp[] = (this.componentMeta?.props ?? [])
      .filter((p) => self[p.name] instanceof Signal)
      .map((p) => ({...p, writable: isWritableSignal(self[p.name]),
        ...(p.name === 'value' ? {default: true} : {})}));
    if (!this.propertyTier)
      return props;
    const own = new Set(props.map((p) => p.name));
    for (const p of this._propertyIndex().values()) {
      if (!own.has(p.name))
        props.push({...propertyFields(p), name: p.name, writable: p.set != null});
    }
    return props;
  }

  private _property(name: string): NamedProperty | undefined {
    return this.propertyTier ? this._propertyIndex().get(name) :
      this.getProperties().find((p) => p.name === name);
  }

  private _propertyIndex(): Map<string, NamedProperty> {
    const u2 = this._u2;
    if (u2.propertyIndex === undefined)
      u2.propertyIndex = new Map(this.getProperties().map((p): [string, NamedProperty] => [p.name, p]));
    return u2.propertyIndex;
  }

  private _propertySignal(name: string): Signal<unknown> | null {
    const u2 = this._u2;
    const cached = u2.propertySignals.get(name);
    if (cached)
      return cached.step;
    const prop = this._property(name);
    if (prop === undefined)
      return null;
    const entry = {known: prop.get?.(this.propertyTarget)} as PropertyEntry;
    entry.inner = signal(entry.known);
    if (prop.set == null)
      entry.step = computed(() => entry.inner.value) as unknown as Signal<unknown>;
    else {
      entry.step = entry.inner;
      this.scope.effect(() => {
        const value = entry.inner.value;
        if (Component.sameValue(entry.known, value))
          return;
        entry.known = value;
        prop.set!(this.propertyTarget, value);
      });
    }
    u2.propertySignals.set(name, entry);
    this._wireChanges();
    return entry.step;
  }

  private _wireChanges(): void {
    const u2 = this._u2;
    if (u2.changesWired)
      return;
    u2.changesWired = true;
    const changes = this.propertyChanges();
    if (changes === null)
      return;
    const subscription = changes.subscribe((name) => this._stale(name));
    this.scope.own(() => subscription.unsubscribe());
  }

  /** A notification never refreshes inside the event: the names pile up and one microtask reads
   * them back — the platform fires on every write, the effect's own included. */
  private _stale(name: PropertyChange): void {
    const u2 = this._u2;
    const first = u2.stale === undefined;
    if (name === null)
      u2.stale = 'all';
    else if (u2.stale === undefined)
      u2.stale = new Set([name]);
    else if (u2.stale !== 'all')
      u2.stale.add(name);
    if (first)
      queueMicrotask(() => this._refreshStale());
  }

  private _refreshStale(): void {
    const u2 = this._u2;
    const stale = u2.stale;
    u2.stale = undefined;
    const scope = u2.scope;
    if (stale === undefined || scope === undefined || scope.isDisposed)
      return;
    for (const name of stale === 'all' ? u2.propertySignals.keys() : stale)
      this._refreshProperty(name);
  }

  private _refreshProperty(name: string): void {
    const entry = this._u2.propertySignals.get(name);
    if (entry === undefined)
      return;
    const value = this._property(name)?.get?.(this.propertyTarget);
    entry.known = value;
    if (!Component.sameValue(entry.inner.peek(), value))
      entry.inner.value = value;
  }

  private _warnOnce(key: string, message: string): void {
    const warned = this._u2.warned ??= new Set();
    if (warned.has(key))
      return;
    warned.add(key);
    console.warn(`u2: ${message}`);
  }

  private static _read(self: Record<string, unknown>, name: string): unknown {
    const member = self[name];
    if (member instanceof Signal)
      return member.peek();
    if (member !== undefined)
      return member;
    const option = (self.options as Record<string, unknown> | undefined)?.[name];
    if (option !== undefined)
      return option;
    const given = (self.specProps as Record<string, unknown> | undefined)?.[name];
    return given instanceof Signal ? given.peek() : given;
  }

  protected registerFunction(f: FuncLike): void {
    this._u2.functions.push(f);
  }

  /** Subscribes `handler` to {@link onEvent}`(eventId)`; it gets the payload and this component.
   * Owned by the ambient scope when there is one (the subscriber's), else by this component; runs
   * as {@link Scope.bind} runs a handler. Answers the unsubscribe. */
  on(eventId: string, handler: (arg: unknown, source: this) => void): () => void {
    const owner = Scope.ambient ?? this.scope;
    const subscription = this.onEvent(eventId).subscribe(owner.bind((x: unknown) => handler(x, this)));
    const off = () => {
      subscription.unsubscribe();
      owner.disown(off);
    };
    owner.own(off);
    return off;
  }

  /** Every matching listener runs inside one batch; a throwing one is reported to this component's
   * scope and does not stop the rest. */
  protected fireEvent(eventId: string, args?: unknown): void {
    const listeners = this._u2.listeners.filter((l) => l.id === null || l.id === eventId);
    if (listeners.length === 0)
      return;
    try {
      batch(() => {
        for (const listener of listeners) {
          try {
            listener.next(args);
          }
          catch (e) {
            this.scope.reportError(e);
          }
        }
      });
    }
    catch (e) {
      this.scope.reportError(e);
    }
  }
}

/** The element → control slot. A registered symbol, so every engine copy on the page reads the same one. */
const SLOT = Symbol.for('u2.control');

/** The `u2-<n>` id counter, shared by every engine copy on the page. */
const ids: {next: number} = (globalThis as Record<PropertyKey, any>)[Symbol.for('u2.ids')] ??= {next: 0};

/** A component that renders: everything visual in u2 is one. */
export class Control extends Component {
  protected _root!: HTMLElement;

  constructor(root?: HTMLElement, scope?: Scope) {
    super(scope);
    this._root = root ?? document.createElement('div');
    this._root.classList.add('u2-component');
    (this._root as unknown as Record<symbol, Control>)[SLOT] = this;
  }

  get root(): HTMLElement {
    return this._root;
  }

  set root(x: HTMLElement) {
    this._root = x;
    (x as unknown as Record<symbol, Control>)[SLOT] = this;
    this._stampName();
  }

  /** Disposes the control; disposed by hand, it also leaves the page — a dead control left on show looks alive. One
   * its owner disposes leaves the DOM to whoever holds it: a region swapping its content, a spec re-rendering. */
  dispose(): void {
    super.dispose();
    if (this._u2.byOwner !== true)
      this.root.remove();
  }

  /** Hidden ⇔ the root carries the `hidden` attribute. */
  get visible(): boolean {
    return !this.root.hidden;
  }

  set visible(x: boolean) {
    this.root.hidden = !x;
  }

  /** Stable element id, minted (`u2-<n>`) on first ask unless the root has one: the anchor of
   * `for`, `aria-*` and `data-u2-owner-id` references. */
  get rootId(): string {
    const root = this.root;
    return root.id || (root.id = `u2-${++ids.next}`);
  }

  /** The state block generic tests read (`DG.U2.state`). This base reads the DOM — kind, name,
   * accessible name, `aria-disabled`, `hidden`, focus, `aria-busy` — plus {@link getItems}; a
   * subclass adds what the DOM cannot say and spreads `super.getState()`. Synchronous, side-effect
   * free, never cached. */
  getState(): IControlState {
    const root = this.root;
    const busy = root.getAttribute('aria-busy');
    const items = this.getItems();
    const enabled = root.getAttribute('aria-disabled') !== 'true' && !root.hasAttribute('disabled');
    return Control._defined({
      kind: root.dataset.u2 ?? this.componentMeta?.tag?.replace(/^u2-/, '') ?? 'control',
      tag: this.componentMeta?.tag ?? root.dataset.u2Tag,
      name: root.dataset.u2Name,
      label: this.accessibleLabel(),
      enabled,
      disabledReason: enabled ? undefined : root.getAttribute('aria-description') || undefined,
      visible: root.isConnected && root.closest('[hidden]') === null &&
        (typeof root.getClientRects !== 'function' || root.getClientRects().length > 0),
      focused: this._hasFocus(),
      busy: busy === null ? undefined : busy === 'true',
      items,
      selection: items?.list.some((item) => item.selected !== undefined) ?
        items.list.filter((item) => item.selected).map((item) => item.key) : undefined,
    });
  }

  /** Items of a collection-like control; undefined for other kinds. `list` holds the first
   * {@link ITEMS_CAP} unless `all`. */
  getItems(_opts?: {all?: boolean}): IItemsState | undefined {
    return undefined;
  }

  /** Brings the item into view and resolves to its element, or null when none matches: the exact
   * key (`data-u2-item`) first, then the label, trimmed and case-insensitive. Several label matches
   * reject, naming the candidates. This base searches the rendered `[data-u2-item]` elements under
   * the root and in the popups it owns; a virtual or collapsible control overrides it to scroll or
   * expand first. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const elements = this._itemElements();
    const byKey = elements.find((el) => el.dataset.u2Item === keyOrLabel);
    if (byKey !== undefined)
      return byKey;
    const wanted = keyOrLabel.trim().toLowerCase();
    const matches = elements.filter((el) => Control._itemLabel(el).toLowerCase() === wanted);
    if (matches.length > 1) {
      const keys = matches.map((el) => el.dataset.u2Item).join(', ');
      throw new Error(`u2: ${matches.length} items match "${keyOrLabel}": ${keys}`);
    }
    return matches[0] ?? null;
  }

  /** The accessible name: `aria-label`, then the `aria-labelledby` texts, then this control's own
   * `data-u2-part="label"` element, then — for a root that is a button — its own text. */
  protected accessibleLabel(): string | undefined {
    const root = this.root;
    const own = Control.labelOf(root);
    if (own !== undefined)
      return own;
    for (const el of root.querySelectorAll('[data-u2-part="label"]')) {
      if (el.closest('.u2-component') === root)
        return Control.text(el) || undefined;
    }
    return root.matches('button, [role="button"]') ? Control.text(root) || undefined : undefined;
  }

  /** `aria-label`, else the joined texts of the `aria-labelledby` elements; undefined when neither names `el`. */
  protected static labelOf(el: Element): string | undefined {
    const aria = el.getAttribute('aria-label');
    if (aria)
      return aria;
    const by = el.getAttribute('aria-labelledby');
    const doc = el.ownerDocument ?? globalThis.document;
    const text = by === null ? '' : by.split(/\s+/)
      .map((id) => doc?.getElementById(id))
      .map((labelEl) => labelEl == null ? '' : Control.text(labelEl))
      .filter((t) => t !== '').join(' ');
    return text || undefined;
  }

  /** The text the user reads: trimmed, single-spaced. */
  protected static text(el: Element): string {
    return (el.textContent ?? '').replace(/\s+/g, ' ').trim();
  }

  private _hasFocus(): boolean {
    const active = globalThis.document?.activeElement;
    if (active == null)
      return false;
    if (this.root.contains(active))
      return true;
    const id = this.root.id;
    return id !== '' && active.closest(Control._ownedBy(id)) !== null;
  }

  private _itemElements(): HTMLElement[] {
    const elements = [...this.root.querySelectorAll<HTMLElement>('[data-u2-item]')];
    const id = this.root.id;
    if (id !== '') {
      for (const popup of globalThis.document.querySelectorAll(Control._ownedBy(id)))
        elements.push(...popup.querySelectorAll<HTMLElement>('[data-u2-item]'));
    }
    return elements;
  }

  /** The selector of the popups a control with this root id owns; no `CSS` in Node. */
  private static _ownedBy(id: string): string {
    return `[data-u2-owner-id="${typeof CSS === 'undefined' ? id : CSS.escape(id)}"]`;
  }

  private static _itemLabel(el: HTMLElement): string {
    return Control.text(el.querySelector('[data-u2-part="item-label"]') ?? el);
  }

  /** The JSON-safe projection behind {@link IControlState.value}: a Date → its ISO string, an object
   * with a string `name` (an entity, a column) → the name, arrays and sets element-wise, plain objects
   * key-wise (depth 4), a signal → its current value; functions and elements → undefined; anything
   * else → `String(x)`. */
  static jsonValue(x: unknown, depth: number = 0): JsonValue | undefined {
    if (x === null || x === undefined)
      return x as null | undefined;
    if (x instanceof Signal)
      return Control.jsonValue(x.peek(), depth);
    switch (typeof x) {
    case 'boolean':
    case 'string':
      return x;
    case 'number':
      return Number.isFinite(x) ? x : String(x);
    case 'function':
    case 'symbol':
      return undefined;
    }
    if (typeof x !== 'object')
      return String(x);
    if (x instanceof Date)
      return Number.isNaN(x.getTime()) ? String(x) : x.toISOString();
    if (typeof Element !== 'undefined' && x instanceof Element)
      return undefined;
    if (depth >= 4)
      return String(x);
    if (Array.isArray(x) || x instanceof Set)
      return [...x].map((item) => Control.jsonValue(item, depth + 1) ?? null);
    const name = (x as {name?: unknown}).name;
    if (typeof name === 'string')
      return name;
    if (!Component._isPlain(x))
      return String(x);
    const out: {[key: string]: JsonValue} = {};
    for (const [key, value] of Object.entries(x)) {
      const projected = Control.jsonValue(value, depth + 1);
      if (projected !== undefined)
        out[key] = projected;
    }
    return out;
  }

  private static _defined<T extends object>(x: T): T {
    for (const key of Object.keys(x) as (keyof T)[]) {
      if (x[key] === undefined)
        delete x[key];
    }
    return x;
  }

  /** The base status plus {@link getState}. */
  getWidgetStatus(): IWidgetStatus {
    return {...super.getWidgetStatus(), state: this.getState()};
  }

  /** A named control carries its name in the DOM: `data-u2-name` is the stable automation id
   * selectors, tutorials and agents address the control by — the same attribute the spec renderer
   * stamps on named nodes, so hand-built and spec-built UI share one selector vocabulary.
   * Derived names (an input's label fallback) never stamp — see {@link Component.deriveName}. */
  get name(): string | undefined {
    return this._u2.name;
  }

  set name(x: string | undefined) {
    this._u2.name = x;
    this._stampName();
  }

  private _stampName(): void {
    const name = this._u2.name;
    if (name === undefined)
      delete this._root.dataset.u2Name;
    else
      this._root.dataset.u2Name = name;
  }

  /** The nearest control owning `el` — the DOM → component door for tooling, automation and
   * inspectors (the u2 counterpart of the platform's `Widget.find`): the nearest root at or above
   * `el`, else, inside a popup portaled to the body, the popup's owner (`data-u2-owner-id`).
   * Best-effort: the last control constructed on (or assigned) a root wins it. The root carries
   * the control under `Symbol.for('u2.control')`, so every engine copy on the page answers alike. */
  static forElement(el: Element | null): Control | undefined {
    for (let node = el; node; node = node.parentElement) {
      const control = (node as unknown as Record<symbol, Control | undefined>)[SLOT];
      if (control !== undefined)
        return control;
    }
    const host = el?.closest<HTMLElement>('[data-u2-owner-id]');
    const owner = host == null ? null : document.getElementById(host.dataset.u2OwnerId!);
    return owner === null ? undefined : Control.forElement(owner);
  }

  static is(x: unknown): x is Control {
    return (x as Partial<Control> | null | undefined)?.root instanceof HTMLElement && Component.is(x);
  }

  /** Component events plus the DOM events of that name raised on {@link root} — the listener is
   * wired per subscription and dropped with it, or with the scope if the subscriber never does. */
  onEvent(eventId: string | null = null): ObservableLike<unknown> {
    const events = super.onEvent(eventId);
    if (eventId === null)
      return events;
    return {
      subscribe: (next: (x: unknown) => void) => {
        const subscription = events.subscribe(next);
        const handler = (e: Event) => next(e);
        this.root.addEventListener(eventId, handler);
        const off = () => this.root.removeEventListener(eventId, handler);
        this.scope.own(off);
        return {unsubscribe: () => {
          subscription.unsubscribe();
          off();
          this.scope.disown(off);
        }};
      },
    };
  }

  /** The function component {@link component} makes, under its old name.
   * @deprecated Use {@link component}: one spelling, apart from an owner's `build(fn)`. */
  static build(builder: (self: Scope) => BuiltRoot, options?: ControlBuildOptions): Control;
  static build<R extends BuiltHandle>(builder: (self: Scope) => R, options?: ControlBuildOptions):
    Control & Omit<R, 'root'>;
  static build(builder: (self: Scope) => unknown, options?: ControlBuildOptions): Control {
    const scope = new Scope(undefined, {lazy: true});
    let built: unknown;
    try {
      built = scope.build(() => builder(scope));
    }
    catch (e) {
      scope.dispose();
      throw e;
    }
    const handle = Control._isHandle(built) ? built : undefined;
    const content = handle === undefined ? built : handle.root;
    const root = content == null ? undefined : content instanceof HTMLElement ? content :
      Control._wrap(content as BuiltRoot);
    const control = new Control(root, scope);
    if (handle !== undefined)
      Control._attach(control, handle);
    if (options?.name !== undefined)
      control.name = options.name;
    return control;
  }

  private static _isHandle(x: unknown): x is {root: unknown} {
    return typeof x === 'object' && x !== null && !Array.isArray(x) && !(x instanceof HTMLElement) &&
      !Control.is(x) && 'root' in x;
  }

  private static _wrap(content: BuiltRoot): HTMLElement {
    const root = document.createElement('div');
    for (const child of Array.isArray(content) ? content : [content])
      root.append(Control.is(child) ? child.root : child);
    return root;
  }

  private static _attach(control: Control, handle: object): void {
    for (const [key, descriptor] of Object.entries(Object.getOwnPropertyDescriptors(handle))) {
      if (key === 'root')
        continue;
      if (key in control) {
        Usage.report('reserved-key', `component()'s handle has \`${key}\`, a member every control has — name it ` +
          `for what it holds (\`${key}Text\`, \`${key}Input\`); it was not copied.`);
        continue;
      }
      Object.defineProperty(control, key, descriptor);
    }
  }
}

/** The content a {@link component} builder answers: an element, a control, a text, or several (in a `div`). */
export type BuiltRoot = HTMLElement | Control | string | (HTMLElement | Control | string)[];

/** A handle a {@link component} builder answers: the content under `root`, beside the API the control carries;
 * no other key may be a Control member. */
export type BuiltHandle = {root: BuiltRoot} & {[K in Exclude<keyof Control, 'root'> | 'fireEvent' | keyof OwnerVerbs]?:
  never};

/** Options of {@link component}. */
export interface ControlBuildOptions {
  /** The automation id, stamped as `data-u2-name`. */
  name?: string;
}

/** Builds a control — the function component: `builder` runs now with the new control's scope ambient and gets
 * it as `self`, so what it builds is owned by the result; after an `await` it gets back in with
 * `self.task(async (t) => … t.build(…))`. It answers the content (an element, a control, a text, or several, which
 * go into a `div`), or a handle `{root, ...api}`: the result is then `Control & api`, the api's properties copied
 * onto it as they are (getters stay live). A handle key that is a Control member is a type error, and reported
 * (`reserved-key`) and skipped at run time. `name` is the automation id (`data-u2-name`). A built control is a
 * root: built with no owner, it is not reported.
 * @example
 * const counter = DG.U2.component((self) => {
 *   const count = DG.U2.signal(0);
 *   const text = ui.divText('0');
 *   self.effect(() => text.textContent = String(count.value));
 *   return {root: text, increment: () => count.value++};
 * }, {name: 'counter'});
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-build} */
export function component(builder: (self: Scope) => BuiltRoot, options?: ControlBuildOptions): Control;
export function component<R extends BuiltHandle>(builder: (self: Scope) => R, options?: ControlBuildOptions):
  Control & Omit<R, 'root'>;
export function component(builder: (self: Scope) => unknown, options?: ControlBuildOptions): Control {
  return Control.build(builder as (self: Scope) => BuiltRoot, options);
}

/** The owner verbs every component has, its {@link Component.scope}'s: `build(fn)` runs `fn` with it ambient, `task`
 * runs async work for it, `subscribe` owns a stream. The engine puts them on every component; their types come with
 * u2 (`@datagrok-libraries/u2` declares them on `Component`), so a platform class that has a member of one of these
 * names of its own (a function view's `build()`) keeps it and still compiles. */
export interface OwnerVerbs {
  /** Runs `fn` now with this owner ambient: what it builds is owned here. See {@link Scope.build}. */
  build<T>(fn: () => T): T;
  /** Async work owned here. See {@link Scope.task}. */
  task<R>(fn: (t: TaskContext) => R | PromiseLike<R>, options?: TaskOptions): TaskRun;
  /** Subscribes `handler` to `source` until this owner is disposed. See {@link Scope.subscribe}. */
  subscribe<T>(source: ObservableLike<T>, handler: (value: T) => void): () => void;
}

/** `task` and `subscribe` are typed on every component, `DG.Widget`, viewers and views included; `build` only with
 * u2, as compute-utils' function views (`DG.ViewBase` subclasses) declare a `build(): void` of their own. */
export interface Component extends Pick<OwnerVerbs, 'task' | 'subscribe'> {}

const OWNER_VERBS: {[K in keyof OwnerVerbs]: (this: Component, ...args: Parameters<OwnerVerbs[K]>) => unknown} = {
  build: function(fn) {
    return this.scope.build(fn);
  },
  task: function(fn, options) {
    return this.scope.task(fn, options);
  },
  subscribe: function(source, handler) {
    return this.scope.subscribe(source, handler);
  },
};
for (const [name, value] of Object.entries(OWNER_VERBS))
  Object.defineProperty(Component.prototype, name, {value, writable: true, configurable: true});
